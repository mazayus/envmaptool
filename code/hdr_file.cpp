/*
 * Copyright 2017 Milan Izai <milan.izai@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "hdr_file.h"

#define HDR_MAX_HEADER_LINE 2048

#define HDR_MIN_RLE_SCANLINE 8
#define HDR_MAX_RLE_SCANLINE 0x7FFF

extern "C" double ldexp(double x, int exp);
extern "C" double frexp(double x, int* exp);

struct rgbe_t
{
    unsigned char r, g, b, e;
};

static FORCE_INLINE void rgbe_pack(rgbe_t* rgbe, float r, float g, float b)
{
    float f = r;
    if (f < g) f = g;
    if (f < b) f = b;

    if (f <= 1e-32)
    {
        rgbe->r = 0;
        rgbe->g = 0;
        rgbe->b = 0;
        rgbe->e = 0;
    }
    else
    {
        int e = 0;
        f = frexp(f, &e) / f;

        rgbe->r = (r > 0) ? (r * 255 * f + 0.5) : 0;
        rgbe->g = (g > 0) ? (g * 255 * f + 0.5) : 0;
        rgbe->b = (b > 0) ? (b * 255 * f + 0.5) : 0;
        rgbe->e = e + 128;
    }
}

static FORCE_INLINE void rgbe_unpack(rgbe_t rgbe, float* r, float* g, float* b)
{
    if (rgbe.e == 0)
    {
        *r = 0.0f;
        *g = 0.0f;
        *b = 0.0f;
    }
    else
    {
        float e = ldexp(1.0, (int) rgbe.e - 128);

        *r = (rgbe.r + 0.5) / 256.0 * e;
        *g = (rgbe.g + 0.5) / 256.0 * e;
        *b = (rgbe.b + 0.5) / 256.0 * e;
    }
}

static bool HDRFileReadHeader(FILE* fp, int* _width, int* _height)
{
    char line[HDR_MAX_HEADER_LINE];

    // TODO: Handle crappy line terminators (\r\n and \r).

    if (!fgets(line, sizeof(line), fp))
    {
        fprintf(stderr, "HDRFileReadHeader: can't read file signature\n");
        return false;
    }

    if (strcmp(line, "#?RADIANCE\n") != 0)
    {
        fprintf(stderr, "HDRFileReadHeader: invalid file signature\n");
        return false;
    }

    for (;;)
    {
        if (!fgets(line, sizeof(line), fp))
        {
            fprintf(stderr, "HDRFileReadHeader: can't read file header\n");
            return false;
        }

        // TODO: At least parse FORMAT and EXPOSURE.

        if (strcmp(line, "\n") == 0)
            break;
    }

    if (!fgets(line, sizeof(line), fp))
    {
        fprintf(stderr, "HDRFileReadHeader: can't read resolution string\n");
        return false;
    }

    // TODO: Support other resolution strings.

    int width = 0, height = 0;
    if (sscanf(line, " -Y %d +X %d", &height, &width) != 2)
    {
        fprintf(stderr, "HDRFileReadHeader: unsupported resolution string\n");
        return false;
    }

    *_width = width;
    *_height = height;

    return true;
}

static bool HDRFileReadNewScanline(FILE* fp, rgbe_t* scanline, int scanline_length, unsigned char lookahead[4])
{
    if (((lookahead[2] << 8) | lookahead[3]) != scanline_length)
        return false;

    for (int component = 0; component < 4; ++component)
    {
        int next_pixel = 0;
        while (next_pixel < scanline_length)
        {
            uint8_t count = 0;
            if (fread(&count, 1, 1, fp) != 1)
                return false;

            if (count > 128)
            {
                count -= 128;

                if (next_pixel + count > scanline_length)
                    return false;

                uint8_t value = 0;
                if (fread(&value, 1, 1, fp) != 1)
                    return false;

                for (int i = 0; i < count; ++i)
                {
                    ((uint8_t*)(scanline + next_pixel + i))[component] = value;
                }

                next_pixel += count;
            }
            else
            {

                if (next_pixel + count > scanline_length)
                    return false;

                for (int i = 0; i < count; ++i)
                {
                    uint8_t value = 0;
                    if (fread(&value, 1, 1, fp) != 1)
                        return false;

                    ((uint8_t*)(scanline + next_pixel + i))[component] = value;
                }

                next_pixel += count;
            }
        }

        assert(next_pixel == scanline_length);
    }

    return true;
}

static bool HDRFileReadOldScanline(FILE* fp, rgbe_t* scanline, int scanline_length, unsigned char lookahead[4])
{
    int shift = 0;

    int next_pixel = 0;

    unsigned char buf[4] = {lookahead[0], lookahead[1], lookahead[2], lookahead[3]};

    goto skip_initial_read;

    while (next_pixel < scanline_length)
    {
        if (fread(buf, 4, 1, fp) != 1)
            return false;

    skip_initial_read:

        if (buf[0] == 1 && buf[1] == 1 && buf[2] == 1)
        {
            int count = ((int) buf[3]) << shift;
            if (next_pixel + count > scanline_length)
                return false;

            rgbe_t pixel = {};
            if (next_pixel > 0)
                pixel = scanline[next_pixel - 1];
            for (int i = 0; i < count; ++i)
                scanline[next_pixel++] = pixel;

            shift += 8;
        }
        else
        {
            scanline[next_pixel].r = buf[0];
            scanline[next_pixel].g = buf[1];
            scanline[next_pixel].b = buf[2];
            scanline[next_pixel].e = buf[3];

            ++next_pixel;

            shift = 0;
        }
    }

    assert(next_pixel == scanline_length);

    return true;
}

static bool HDRFileReadScanline(FILE* fp, rgbe_t* scanline, int scanline_length)
{
    unsigned char lookahead[4];
    if (fread(lookahead, 4, 1, fp) != 1)
        return false;

    if (scanline_length >= HDR_MIN_RLE_SCANLINE && scanline_length <= HDR_MAX_RLE_SCANLINE)
    {
        if (lookahead[0] == 2 && lookahead[1] == 2 && !(lookahead[2] & 0x80))
        {
            if (!HDRFileReadNewScanline(fp, scanline, scanline_length, lookahead))
                return false;
        }
        else
        {
            if (!HDRFileReadOldScanline(fp, scanline, scanline_length, lookahead))
                return false;
        }
    }
    else
    {
        if (!HDRFileReadOldScanline(fp, scanline, scanline_length, lookahead))
            return false;
    }

    return true;
}

bool HDRFileRead(const char* filename, float** _image, int* _width, int* _height)
{
    FILE* fp = fopen(filename, "rb");
    if (!fp)
    {
        fprintf(stderr, "HDRFileRead: can't open file\n");
        return false;
    }

    int width = 0, height = 0;
    if (!HDRFileReadHeader(fp, &width, &height))
    {
        fprintf(stderr, "HDRFileRead: can't read file header\n");
        fclose(fp);
        return false;
    }

    rgbe_t* rgbe = new rgbe_t[width*height];
    for (int i = 0; i < height; ++i)
    {
        if (!HDRFileReadScanline(fp, rgbe+i*width, width))
        {
            fprintf(stderr, "HDRFileRead: can't read scanline %d\n", i);
            delete[] rgbe;
            fclose(fp);
            return false;
        }
    }

    float* image = new float[width * height * 3];
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            float* pixel = image + (i * width + j) * 3;
            rgbe_unpack(rgbe[i * width + j], &pixel[0], &pixel[1], &pixel[2]);
        }
    }

    *_image = image;
    *_width = width;
    *_height = height;

    delete[] rgbe;
    fclose(fp);
    return true;
}

static void HDRFileWriteHeader(FILE* fp, int width, int height)
{
    fprintf(fp, "#?RADIANCE\n");
    fprintf(fp, "# Created by EnvMapTool\n");
    fprintf(fp, "FORMAT=32-bit_rle_rgbe\n");
    fprintf(fp, "EXPOSURE=%g\n", 1.0);
    fprintf(fp, "\n");
    fprintf(fp, "-Y %d +X %d\n", height, width);
}

static void HDRFileWriteScanline(FILE* fp, const rgbe_t* scanline, int scanline_length)
{
    // TODO: RLE compression

    if (scanline_length >= HDR_MIN_RLE_SCANLINE && scanline_length <= HDR_MAX_RLE_SCANLINE)
    {
        // new style scanline

        unsigned char magic[4] = {2, 2, 0, 0};
        magic[2] = (scanline_length >> 8) & 0xFF;
        magic[3] = scanline_length & 0xFF;
        fwrite(magic, sizeof(magic), 1, fp);

        for (int component = 0; component < 4; ++component)
        {
            int next_pixel = 0;
            while (next_pixel < scanline_length)
            {
                unsigned char count = (scanline_length - next_pixel) <= 128 ? (scanline_length - next_pixel) : 128;
                fwrite(&count, 1, 1, fp);

                for (int i = 0; i < count; ++i)
                    fwrite((unsigned char *) (scanline + next_pixel + i) + component, 1, 1, fp);

                next_pixel += count;
            }
        }
    }
    else
    {
        // old style scanline

        fwrite(scanline, sizeof(rgbe_t), scanline_length, fp);
    }
}

bool HDRFileWrite(const char* filename, const float* image, int width, int height)
{
    FILE* fp = fopen(filename, "wb");
    if (!fp)
    {
        fprintf(stderr, "HDRFileWrite: can't open file\n");
        return false;
    }

    rgbe_t* rgbe = new rgbe_t[width*height];
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            const float* pixel = image + (i * width + j) * 3;
            rgbe_pack(rgbe + i * width + j, pixel[0], pixel[1], pixel[2]);
        }
    }

    HDRFileWriteHeader(fp, width, height);

    for (int i = 0; i < height; ++i)
        HDRFileWriteScanline(fp, rgbe + i * width, width);

    delete[] rgbe;
    fclose(fp);
    return true;
}
