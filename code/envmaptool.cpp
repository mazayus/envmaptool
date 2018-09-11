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

#include <SDL.h>

#include "common.h"
#include "hdr_file.h"
#include "imgui.h"
#include "math.h"
#include "opengl.h"

#define WINDOW_TITLE    "EnvMapTool"

#define INITIAL_WINDOW_WIDTH    1366
#define INITIAL_WINDOW_HEIGHT   768

#define MILLISECONDS_PER_FRAME  16

// NOTE: Ideally FILTER_MS_PER_FRAME should be less than MILLISECONDS_PER_FRAME, but this leads to
// a noticable increase in the total running time of the filter. Higher values of FILTER_MS_PER_FRAME
// result in faster filtering, but less responsive UI.
#define FILTER_MS_PER_FRAME     31

static SDL_Window* sdl_window;
static SDL_GLContext sdl_glcontext;
static int window_width = INITIAL_WINDOW_WIDTH;
static int window_height = INITIAL_WINDOW_HEIGHT;

static bool Init();
static void Shutdown();

static GLProgram imgui_program;
static GLuint imgui_texture;
static GLuint imgui_vao;
static GLuint imgui_vertex_buffer;
static GLuint imgui_index_buffer;

static void InitImGui();
static void SendEventToImGui(const SDL_Event* event);
static void RenderImGui();

////////////////////////////////////////////////////////////////////////////////

struct CubeMap
{
    int             face_size;
    float*          face_pixels[6];

    GLuint          gl_id; // NOTE: all cube maps in a mip chain should have the same gl_id

    CubeMap*        next_mip;
};

CubeMap* AllocateCubeMap(int mip_map_count, int face_size)
{
    assert(mip_map_count >= 1);
    assert(face_size >= 1);

    CubeMap* result = new CubeMap;

    glActiveTexture(GL_TEXTURE15);
    glGenTextures(1, &result->gl_id);
    glBindTexture(GL_TEXTURE_CUBE_MAP, result->gl_id);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, mip_map_count == 1 ? GL_LINEAR : GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    CubeMap* cubemap = result;

    for (int mip = 0; mip < mip_map_count; ++mip)
    {
        cubemap->face_size = face_size;

        for (int face = 0; face < 6; ++face)
            cubemap->face_pixels[face] = new float[face_size*face_size*3];

        cubemap->gl_id = result->gl_id;

        if (face_size == 1 || mip == mip_map_count - 1)
        {
            cubemap->next_mip = NULL;
            break;
        }

        cubemap->next_mip = new CubeMap;
        cubemap = cubemap->next_mip;

        face_size /= 2;
    }

    glBindTexture(GL_TEXTURE_CUBE_MAP, 0);

    return result;
}

void FreeCubeMap(CubeMap* cubemap)
{
    glDeleteTextures(1, &cubemap->gl_id);

    while (cubemap)
    {
        for (int face = 0; face < 6; ++face)
            delete[] cubemap->face_pixels[face];

        CubeMap* next_mip = cubemap->next_mip;

        delete cubemap;

        cubemap = next_mip;
    }
}

void UploadCubeMapToGPU(CubeMap* cubemap)
{
    glActiveTexture(GL_TEXTURE15);
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap->gl_id);

    int mip = 0;
    while (cubemap)
    {
        for (int face = 0; face < 6; ++face)
        {
            glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + face, mip, GL_RGBA32F,
                         cubemap->face_size, cubemap->face_size, 0,
                         GL_RGB, GL_FLOAT, cubemap->face_pixels[face]);
        }

        ++mip;
        cubemap = cubemap->next_mip;
    }

    glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
}

void DownloadCubeMapFromGPU(CubeMap* cubemap)
{
    glActiveTexture(GL_TEXTURE15);
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap->gl_id);

    int mip = 0;
    while (cubemap)
    {
        for (int face = 0; face < 6; ++face)
        {
            glGetTexImage(GL_TEXTURE_CUBE_MAP_POSITIVE_X + face, mip,
                          GL_RGB, GL_FLOAT, cubemap->face_pixels[face]);
        }

        ++mip;
        cubemap = cubemap->next_mip;
    }

    glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
}

Vector3 SampleCubeMap(const CubeMap* cubemap, float theta, float phi)
{
    float cos_theta = cos(theta);
    float sin_theta = sin(theta);
    float cos_phi = cos(phi);
    float sin_phi = sin(phi);

    Vector3 n = Vector3(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
    Vector3 abs_n = Vector3(Math::Abs(n.x), Math::Abs(n.y), Math::Abs(n.z));

    int face = 0;
    float face_u = 0;
    float face_v = 0;

    if (abs_n.x >= abs_n.y && abs_n.x >= abs_n.z)
    {
        if (n.x > 0) // +X
        {
            face = 0;
            face_u = -n.z / abs_n.x;
            face_v = -n.y / abs_n.x;
        }
        else // -X
        {
            face = 1;
            face_u =  n.z / abs_n.x;
            face_v = -n.y / abs_n.x;
        }
    }
    else if (abs_n.y >= abs_n.z)
    {
        if (n.y > 0) // +Y
        {
            face = 2;
            face_u =  n.x / abs_n.y;
            face_v =  n.z / abs_n.y;
        }
        else // -Y
        {
            face = 3;
            face_u =  n.x / abs_n.y;
            face_v = -n.z / abs_n.y;
        }
    }
    else
    {
        if (n.z > 0) // +Z
        {
            face = 4;
            face_u =  n.x / abs_n.z;
            face_v = -n.y / abs_n.z;
        }
        else // -Z
        {
            face = 5;
            face_u = -n.x / abs_n.z;
            face_v = -n.y / abs_n.z;
        }
    }

    face_u = face_u * 0.5 + 0.5;
    face_v = face_v * 0.5 + 0.5;

    int face_x = face_u * (cubemap->face_size - 1);
    int face_y = face_v * (cubemap->face_size - 1);

    const float* pixel = &cubemap->face_pixels[face][(face_y * cubemap->face_size + face_x) * 3];
    return Vector3(pixel[0], pixel[1], pixel[2]);
}

void FillCubeMapWithFaceColors(CubeMap* cubemap, const float face_colors[6][3])
{
    while (cubemap)
    {
        for (int face = 0; face < 6; ++face)
        {
            float* pixels = cubemap->face_pixels[face];

            for (int i = 0; i < cubemap->face_size; ++i)
            {
                for (int j = 0; j < cubemap->face_size; ++j)
                {
                    pixels[(i * cubemap->face_size + j) * 3    ] = face_colors[face][0];
                    pixels[(i * cubemap->face_size + j) * 3 + 1] = face_colors[face][1];
                    pixels[(i * cubemap->face_size + j) * 3 + 2] = face_colors[face][2];
                }
            }
        }

        cubemap = cubemap->next_mip;
    }
}

enum CubeMapLayout
{
    CUBEMAP_LAYOUT_INVALID,
    CUBEMAP_LAYOUT_HORIZONTAL_CROSS,
    CUBEMAP_LAYOUT_VERTICAL_CROSS,
    CUBEMAP_LAYOUT_HORIZONTAL_STRIP,
    CUBEMAP_LAYOUT_VERTICAL_STRIP,
};

static void CopyImageBlockToCubeMap(CubeMap* cubemap, int face, const float* image, int stride,
                                    int block_x_index, int block_y_index, bool flip_x, bool flip_y)
{
    int x_offset = flip_x ? ((block_x_index + 1) * cubemap->face_size - 1) : (block_x_index * cubemap->face_size);
    int y_offset = flip_y ? ((block_y_index + 1) * cubemap->face_size - 1) : (block_y_index * cubemap->face_size);
    int x_scale = flip_x ? -1 : 1;
    int y_scale = flip_y ? -1 : 1;

    for (int i = 0; i < cubemap->face_size; ++i)
    {
        for (int j = 0; j < cubemap->face_size; ++j)
        {
            int dest = i * cubemap->face_size + j;
            int src = (y_offset + y_scale * i) * stride + (x_offset + x_scale * j);

            cubemap->face_pixels[face][dest*3    ] = image[src*3    ];
            cubemap->face_pixels[face][dest*3 + 1] = image[src*3 + 1];
            cubemap->face_pixels[face][dest*3 + 2] = image[src*3 + 2];
        }
    }
}

static void CopyCubeMapFaceToImage(const CubeMap* cubemap, int face, float* image, int stride,
                                   int block_x_index, int block_y_index, bool flip_x, bool flip_y)
{
    int x_offset = flip_x ? ((block_x_index + 1) * cubemap->face_size - 1) : (block_x_index * cubemap->face_size);
    int y_offset = flip_y ? ((block_y_index + 1) * cubemap->face_size - 1) : (block_y_index * cubemap->face_size);
    int x_scale = flip_x ? -1 : 1;
    int y_scale = flip_y ? -1 : 1;

    for (int i = 0; i < cubemap->face_size; ++i)
    {
        for (int j = 0; j < cubemap->face_size; ++j)
        {
            int dest = i * cubemap->face_size + j;
            int src = (y_offset + y_scale * i) * stride + (x_offset + x_scale * j);

            image[src*3    ] = cubemap->face_pixels[face][dest*3    ];
            image[src*3 + 1] = cubemap->face_pixels[face][dest*3 + 1];
            image[src*3 + 2]  =cubemap->face_pixels[face][dest*3 + 2];
        }
    }
}

CubeMap* LoadCubeMap(const char* filename, CubeMapLayout* layout)
{
    if (!filename[0])
    {
        fprintf(stderr, "LoadCubeMap: empty filename\n");
        return NULL;
    }

    float* image;
    int width, height;
    if (!HDRFileRead(filename, &image, &width, &height))
    {
        fprintf(stderr, "LoadCubeMap: can't load image from '%s'\n", filename);
        return NULL;
    }

    if (width == 6 * height)
    {
        if (layout) *layout = CUBEMAP_LAYOUT_HORIZONTAL_STRIP;

        int face_size = height;

        CubeMap* cubemap = AllocateCubeMap(1, face_size);
        CopyImageBlockToCubeMap(cubemap, 0, image, width, 0, 0, 1, 0);
        CopyImageBlockToCubeMap(cubemap, 1, image, width, 1, 0, 1, 0);
        CopyImageBlockToCubeMap(cubemap, 2, image, width, 2, 0, 0, 1);
        CopyImageBlockToCubeMap(cubemap, 3, image, width, 3, 0, 0, 1);
        CopyImageBlockToCubeMap(cubemap, 4, image, width, 4, 0, 1, 0);
        CopyImageBlockToCubeMap(cubemap, 5, image, width, 5, 0, 1, 0);

        delete[] image;

        UploadCubeMapToGPU(cubemap);

        return cubemap;

    }
    else if (height == 6 * width)
    {
        if (layout) *layout = CUBEMAP_LAYOUT_VERTICAL_STRIP;

        int face_size = width;

        CubeMap* cubemap = AllocateCubeMap(1, face_size);
        CopyImageBlockToCubeMap(cubemap, 0, image, width, 0, 0, 1, 0);
        CopyImageBlockToCubeMap(cubemap, 1, image, width, 0, 1, 1, 0);
        CopyImageBlockToCubeMap(cubemap, 2, image, width, 0, 2, 0, 1);
        CopyImageBlockToCubeMap(cubemap, 3, image, width, 0, 3, 0, 1);
        CopyImageBlockToCubeMap(cubemap, 4, image, width, 0, 4, 1, 0);
        CopyImageBlockToCubeMap(cubemap, 5, image, width, 0, 5, 1, 0);

        delete[] image;

        UploadCubeMapToGPU(cubemap);

        return cubemap;
    }
    else if (3 * width == 4 * height)
    {
        if (layout) *layout = CUBEMAP_LAYOUT_HORIZONTAL_CROSS;

        int face_size = height / 3;

        CubeMap* cubemap = AllocateCubeMap(1, face_size);
        CopyImageBlockToCubeMap(cubemap, 0, image, width, 2, 1, 1, 0);
        CopyImageBlockToCubeMap(cubemap, 1, image, width, 0, 1, 1, 0);
        CopyImageBlockToCubeMap(cubemap, 2, image, width, 1, 0, 0, 1);
        CopyImageBlockToCubeMap(cubemap, 3, image, width, 1, 2, 0, 1);
        CopyImageBlockToCubeMap(cubemap, 4, image, width, 3, 1, 1, 0);
        CopyImageBlockToCubeMap(cubemap, 5, image, width, 1, 1, 1, 0);

        delete[] image;

        UploadCubeMapToGPU(cubemap);

        return cubemap;
    }
    else if (4 * width == 3 * height)
    {
        if (layout) *layout = CUBEMAP_LAYOUT_VERTICAL_CROSS;

        int face_size = width / 3;

        CubeMap* cubemap = AllocateCubeMap(1, face_size);
        CopyImageBlockToCubeMap(cubemap, 0, image, width, 2, 1, 1, 0);
        CopyImageBlockToCubeMap(cubemap, 1, image, width, 0, 1, 1, 0);
        CopyImageBlockToCubeMap(cubemap, 2, image, width, 1, 0, 0, 1);
        CopyImageBlockToCubeMap(cubemap, 3, image, width, 1, 2, 0, 1);
        CopyImageBlockToCubeMap(cubemap, 4, image, width, 1, 3, 0, 1);
        CopyImageBlockToCubeMap(cubemap, 5, image, width, 1, 1, 1, 0);

        delete[] image;

        UploadCubeMapToGPU(cubemap);

        return cubemap;
    }
    else
    {
        delete[] image;
        fprintf(stderr, "LoadCubeMap: unsupported image layout\n");
        return NULL;
    }
}

void SaveCubeMap(const CubeMap* cubemap, const char* filename, CubeMapLayout layout)
{
    // TODO: handle crappy path separators (\)
    const char* basename = strrchr(filename, '/');
    basename = basename ? basename + 1 : filename;
    if (!basename[0])
    {
        fprintf(stderr, "SaveCubeMap: empty filename\n");
        return;
    }

    const char* extension = strrchr(basename, '.');
    if (!extension) extension = filename + strlen(filename);

    // NOTE: Allocate enough memory to hold any cubemap layout.
    float* image = new float[cubemap->face_size * cubemap->face_size * 12 * 3];
    memset(image, 0, cubemap->face_size * cubemap->face_size * 12 * 3 * sizeof(float));

    int mip = 0;
    while (cubemap)
    {
        int width, height;

        switch (layout)
        {
        case CUBEMAP_LAYOUT_HORIZONTAL_CROSS:
        {
            width = cubemap->face_size * 4;
            height = cubemap->face_size * 3;

            CopyCubeMapFaceToImage(cubemap, 0, image, width, 2, 1, 1, 0);
            CopyCubeMapFaceToImage(cubemap, 1, image, width, 0, 1, 1, 0);
            CopyCubeMapFaceToImage(cubemap, 2, image, width, 1, 0, 0, 1);
            CopyCubeMapFaceToImage(cubemap, 3, image, width, 1, 2, 0, 1);
            CopyCubeMapFaceToImage(cubemap, 4, image, width, 3, 1, 1, 0);
            CopyCubeMapFaceToImage(cubemap, 5, image, width, 1, 1, 1, 0);

            break;
        }
        case CUBEMAP_LAYOUT_VERTICAL_CROSS:
        {
            width = cubemap->face_size * 3;
            height = cubemap->face_size * 4;

            CopyCubeMapFaceToImage(cubemap, 0, image, width, 2, 1, 1, 0);
            CopyCubeMapFaceToImage(cubemap, 1, image, width, 0, 1, 1, 0);
            CopyCubeMapFaceToImage(cubemap, 2, image, width, 1, 0, 0, 1);
            CopyCubeMapFaceToImage(cubemap, 3, image, width, 1, 2, 0, 1);
            CopyCubeMapFaceToImage(cubemap, 4, image, width, 1, 3, 0, 1);
            CopyCubeMapFaceToImage(cubemap, 5, image, width, 1, 1, 1, 0);

            break;
        }
        case CUBEMAP_LAYOUT_HORIZONTAL_STRIP:
        {
            width = cubemap->face_size * 6;
            height = cubemap->face_size;

            CopyCubeMapFaceToImage(cubemap, 0, image, width, 0, 0, 1, 0);
            CopyCubeMapFaceToImage(cubemap, 1, image, width, 1, 0, 1, 0);
            CopyCubeMapFaceToImage(cubemap, 2, image, width, 2, 0, 0, 1);
            CopyCubeMapFaceToImage(cubemap, 3, image, width, 3, 0, 0, 1);
            CopyCubeMapFaceToImage(cubemap, 4, image, width, 4, 0, 1, 0);
            CopyCubeMapFaceToImage(cubemap, 5, image, width, 5, 0, 1, 0);

            break;
        }
        case CUBEMAP_LAYOUT_VERTICAL_STRIP:
        {
            width = cubemap->face_size;
            height = cubemap->face_size * 6;

            CopyCubeMapFaceToImage(cubemap, 0, image, width, 0, 0, 1, 0);
            CopyCubeMapFaceToImage(cubemap, 1, image, width, 0, 1, 1, 0);
            CopyCubeMapFaceToImage(cubemap, 2, image, width, 0, 2, 0, 1);
            CopyCubeMapFaceToImage(cubemap, 3, image, width, 0, 3, 0, 1);
            CopyCubeMapFaceToImage(cubemap, 4, image, width, 0, 4, 1, 0);
            CopyCubeMapFaceToImage(cubemap, 5, image, width, 0, 5, 1, 0);

            break;
        }
        default:
            INVALID_CODE_PATH;
        }

        char fullname[1024];
        snprintf(fullname, sizeof(fullname), "%.*s_mip%d%s", (int)(extension-filename), filename, mip, extension);

        if (!HDRFileWrite(fullname, image, width, height))
        {
            fprintf(stderr, "SaveCubeMap: can't save image to '%s'\n", fullname);
            // NOTE: Keep going even if we failed to write a mip map.
        }

        cubemap = cubemap->next_mip;
        ++mip;
    }

    delete[] image;
}

struct CubeMapFace
{
    Vector3 x;
    Vector3 y;
    Vector3 z;
};

static const CubeMapFace g_CubeMapFaces[6] = {
    {Vector3( 0,  0, -1), Vector3( 0, -1,  0), Vector3( 1,  0,  0)},
    {Vector3( 0,  0,  1), Vector3( 0, -1,  0), Vector3(-1,  0,  0)},
    {Vector3( 1,  0,  0), Vector3( 0,  0,  1), Vector3( 0,  1,  0)},
    {Vector3( 1,  0,  0), Vector3( 0,  0, -1), Vector3( 0, -1,  0)},
    {Vector3( 1,  0,  0), Vector3( 0, -1,  0), Vector3( 0,  0,  1)},
    {Vector3(-1,  0,  0), Vector3( 0, -1,  0), Vector3( 0,  0, -1)}
};

static float RadicalInverse(uint32_t bits)
{
    bits = (bits << 16) | (bits >> 16);
    bits = ((bits & 0x00FF00FFu) << 8) | ((bits & 0xFF00FF00u) >> 8);
    bits = ((bits & 0x0F0F0F0Fu) << 4) | ((bits & 0xF0F0F0F0u) >> 4);
    bits = ((bits & 0x33333333u) << 2) | ((bits & 0xCCCCCCCCu) >> 2);
    bits = ((bits & 0x55555555u) << 1) | ((bits & 0xAAAAAAAAu) >> 1);

    const float InverseMaxU32 = 2.3283064365386963e-10;
    return float(bits) * InverseMaxU32;
}

static void Hammersley(uint32_t i, uint32_t N, float* xi1, float* xi2)
{
    *xi1 = float(i) / float(N);
    *xi2 = RadicalInverse(i);
}

static Vector3 CartesianToSpherical(Vector3 p)
{
    float r = Math::Length(p);
    if (Math::Abs(r) < 1e-5) return Vector3(0, 0, 0);

    float theta = acos(p.z / r);
    float phi = atan2(p.y / r, p.x / r);
    if (phi < 0) phi += 2 * Math::PI;

    return Vector3(theta, phi, r);
}

static float Factorial(int n)
{
    float result = 1.0f;
    for (int i = 2; i <= n; ++i)
        result *= i;
    return result;
}

static float sh_P(int l, int m, float x)
{
    assert(l >= 0 && m >= 0 && m <= l);

    float pmm = 1.0;

    if (m != 0)
    {
        float s = sqrt(1 - x*x);
        for (int i = 1; i <= m; ++i)
            pmm *= -(2*i-1) * s;
    }

    if (l == m)
        return pmm;

    float pm1m = x * (2*m+1) * pmm;

    if (l == m + 1)
        return pm1m;

    float plm = 0.0;
    for (int i = m + 2; i <= l; ++i)
    {
        plm = ((2*i-1)*x*pm1m - (i-1+m)*pmm) / (i-m);
        pmm = pm1m;
        pm1m = plm;
    }

    return plm;
}

static float sh_K(int l, int m)
{
    assert(l >= 0 && m >= 0 && m <= l);

    return sqrt((2*l+1) / (4*Math::PI) * Factorial(l-m) / Factorial(l+m));
}

static float sh_Y(int l, int m, float theta, float phi)
{
    assert(l >= 0 && -l <= m && m <= l);

    if (m > 0)
        return Math::SQRT_2 * sh_K(l, m) * cos(m*phi) * sh_P(l, m, cos(theta));
    else if (m < 0)
        return Math::SQRT_2 * sh_K(l, -m) * sin(-m*phi) * sh_P(l, -m, cos(theta));
    else
        return sh_K(l, 0) * sh_P(l, 0, cos(theta));
}

static int sh_Index(int l, int m)
{
    return l * l + l + m;
}

static int sh_Count(int lmax)
{
    return (lmax + 1) * (lmax + 1);
}

void ProjectCubeMapOntoSH(CubeMap* cubemap, int degree, float* shcoeff_r, float* shcoeff_g, float* shcoeff_b)
{
    int num_coeffs = sh_Count(degree);
    for (int i = 0; i < num_coeffs; ++i)
    {
        shcoeff_r[i] = 0.0f;
        shcoeff_g[i] = 0.0f;
        shcoeff_b[i] = 0.0f;
    }

    const int num_samples = 10000;
    for (int i = 0; i < num_samples; ++i)
    {
        float xi1, xi2;
        Hammersley(i, num_samples, &xi1, &xi2);

        float theta = acos(1.0 - 2.0 * xi1);
        float phi = 2.0 * Math::PI * xi2;

        Vector3 pixel = SampleCubeMap(cubemap, theta, phi);

        for (int l = 0; l <= degree; ++l)
        {
            for (int m = -l; m <= l; ++m)
            {
                float y_lm = sh_Y(l, m, theta, phi);
                shcoeff_r[sh_Index(l, m)] += pixel.r * y_lm;
                shcoeff_g[sh_Index(l, m)] += pixel.g * y_lm;
                shcoeff_b[sh_Index(l, m)] += pixel.b * y_lm;
            }
        }
    }

    for (int l = 0; l <= degree; ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            float weight = 4 * Math::PI / num_samples;
            shcoeff_r[sh_Index(l, m)] *= weight;
            shcoeff_g[sh_Index(l, m)] *= weight;
            shcoeff_b[sh_Index(l, m)] *= weight;
        }
    }
}

void ReconstructCubeMapFromSH(CubeMap* cubemap, int degree, float* shcoeff_r, float* shcoeff_g, float* shcoeff_b)
{
    while (cubemap)
    {
        for (int face = 0; face < 6; ++face)
        {
            for (int i = 0; i < cubemap->face_size; ++i)
            {
                for (int j = 0; j < cubemap->face_size; ++j)
                {
                    Vector3 v = g_CubeMapFaces[face].z
                        + g_CubeMapFaces[face].y * ((i+0.5) / (float) cubemap->face_size * 2 - 1)
                        + g_CubeMapFaces[face].x * ((j+0.5) / (float) cubemap->face_size * 2 - 1);

                    Vector3 spherical = CartesianToSpherical(Math::Normalize(v));

                    float theta = spherical.x;
                    float phi = spherical.y;

                    float r = 0.0f;
                    float g = 0.0f;
                    float b = 0.0f;

                    for (int l = 0; l <= degree; ++l)
                    {
                        for (int m = -l; m <= l; ++m)
                        {
                            float y_lm = sh_Y(l, m, theta, phi);
                            r += shcoeff_r[sh_Index(l, m)] * y_lm;
                            g += shcoeff_g[sh_Index(l, m)] * y_lm;
                            b += shcoeff_b[sh_Index(l, m)] * y_lm;
                        }
                    }

                    float* pixel = &cubemap->face_pixels[face][(i * cubemap->face_size + j) * 3];

                    pixel[0] = r > 0 ? r : 0;
                    pixel[1] = g > 0 ? g : 0;
                    pixel[2] = b > 0 ? b : 0;
                }
            }
        }

        cubemap = cubemap->next_mip;
    }
}

////////////////////////////////////////////////////////////////////////////////

struct Camera
{
    Transform       transform;
    float           fovy, aspect, znear, zfar;
};

struct Mesh
{
    GLuint vao;

    int num_vertices;
    GLuint vertex_buffer;

    int num_indices;
    GLuint index_buffer;
};

enum FilterType
{
    FILTER_TYPE_COSINE,
    FILTER_TYPE_PHONG_COSINE,
    FILTER_TYPE_BECKMANN_COSINE,
    FILTER_TYPE_GGX_COSINE,
};

struct FilterSettings
{
    FilterType type;

    int face_size;
    int max_mip_level;

    bool sh_approximation;
    int sh_max_degree;
    float sh_window_strength;

    int num_samples;

    float specular_exponent;
    float specular_exponent_divisor;

    float roughness;
    float roughness_increment;
};

struct DisplaySettings
{
    int env_map;        // 0 - source, 1 - filtered
    int tex_coord;      // 0 - normal, 1 - reflection

    float lod;

    float exposure;
};

struct FilterTaskProgress
{
    FilterSettings filter_params;

    bool is_running;
    int num_samples_completed;
    int num_samples_per_frame;

    GLuint frame_time_query;
};

struct EnvMapTool
{
    Camera camera;
    FilterSettings filter;
    DisplaySettings display;

    GLProgram draw_mesh_program;
    GLProgram filter_cosine_program;
    GLProgram filter_phong_cosine_program;
    GLProgram filter_beckmann_cosine_program;
    GLProgram filter_ggx_cosine_program;

    GLuint filter_vao;
    GLuint filter_fbo;

    CubeMap* source_env_map;
    CubeMap* filtered_env_map;

    CubeMapLayout layout;

    Mesh sphere_mesh;

    FilterTaskProgress filter_task;
};

static void InitEnvMapTool(EnvMapTool* tool);
static void UpdateEnvMapTool(EnvMapTool* tool, int buttons, int x, int y, int dwheel, int dx, int dy);

int main(int argc, char* argv[])
{
    (void) argc; (void) argv;

    if (!Init())
        return 1;

    InitImGui();

    EnvMapTool env_map_tool = {};
    InitEnvMapTool(&env_map_tool);

    bool should_quit = false;

    uint32_t last_frame_time = SDL_GetTicks();

    int mouse_buttons = 0, mouse_x = 0, mouse_y = 0, mouse_dwheel = 0, mouse_dx = 0, mouse_dy = 0;

    while (!should_quit)
    {
        uint32_t time = SDL_GetTicks();
        uint32_t dt = time - last_frame_time;

        if (dt < MILLISECONDS_PER_FRAME)
        {
            SDL_Delay(MILLISECONDS_PER_FRAME - dt);

            time = SDL_GetTicks();
            dt = time - last_frame_time;
        }

        last_frame_time = time;

        mouse_dwheel = 0; mouse_dx = 0; mouse_dy = 0;

        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            switch (event.type)
            {
            case SDL_QUIT:
                should_quit = true;
                break;
            case SDL_KEYDOWN:
                break;
            case SDL_KEYUP:
                break;
            case SDL_MOUSEBUTTONDOWN:
                switch (event.button.button)
                {
                case SDL_BUTTON_LEFT:
                    mouse_buttons |= BIT(0);
                    break;
                case SDL_BUTTON_RIGHT:
                    mouse_buttons |= BIT(1);
                    break;
                case SDL_BUTTON_MIDDLE:
                    mouse_buttons |= BIT(2);
                    break;
                }
                break;
            case SDL_MOUSEBUTTONUP:
                switch (event.button.button)
                {
                case SDL_BUTTON_LEFT:
                    mouse_buttons &= ~BIT(0);
                    break;
                case SDL_BUTTON_RIGHT:
                    mouse_buttons &= ~BIT(1);
                    break;
                case SDL_BUTTON_MIDDLE:
                    mouse_buttons &= ~BIT(2);
                    break;
                }
                break;
            case SDL_MOUSEMOTION:
                mouse_x = event.motion.x;
                mouse_y = event.motion.y;
                mouse_dx += event.motion.xrel;
                mouse_dy += event.motion.yrel;
                break;
            case SDL_MOUSEWHEEL:
                mouse_dwheel += event.wheel.y * ((event.wheel.direction == SDL_MOUSEWHEEL_NORMAL) ? 1 : -1);
                break;
            case SDL_WINDOWEVENT:
                if (event.window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
                {
                    window_width = event.window.data1;
                    window_height = event.window.data2;
                }
                break;
            }

            SendEventToImGui(&event);
        }

        ImGui::NewFrame();

        ImGui::GetIO().DeltaTime = dt / 1000.0f;

        UpdateEnvMapTool(&env_map_tool, mouse_buttons, mouse_x, mouse_y, mouse_dwheel, mouse_dx, mouse_dy);

        RenderImGui();

        SDL_GL_SwapWindow(sdl_window);
    }

    Shutdown();
    return 0;
}

#if DEBUG_OPENGL

static void APIENTRY GLDebugCallback(GLenum source, GLenum type, GLuint /*id*/, GLenum severity,
                                     GLsizei length, const GLchar* message, const void* /*userdata*/)
{
    const char* severity_str = NULL;
    switch (severity)
    {
    case GL_DEBUG_SEVERITY_HIGH_ARB:
        severity_str = "HIGH";
        break;
    case GL_DEBUG_SEVERITY_MEDIUM_ARB:
        severity_str = "MEDIUM";
        break;
    case GL_DEBUG_SEVERITY_LOW_ARB:
        severity_str = "LOW";
        break;
    default:
        severity_str = "UNKNOWN";
        break;
    }

    const char* source_str = NULL;
    switch (source)
    {
    case GL_DEBUG_SOURCE_API_ARB:
        source_str = "API";
        break;
    case GL_DEBUG_SOURCE_WINDOW_SYSTEM_ARB:
        source_str = "WINDOW_SYSTEM";
        break;
    case GL_DEBUG_SOURCE_SHADER_COMPILER_ARB:
        source_str = "SHADER_COMPILER";
        break;
    case GL_DEBUG_SOURCE_THIRD_PARTY_ARB:
        source_str = "THIRD_PARTY";
        break;
    case GL_DEBUG_SOURCE_APPLICATION_ARB:
        source_str = "APPLICATION";
        break;
    case GL_DEBUG_SOURCE_OTHER_ARB:
        source_str = "OTHER";
        break;
    default:
        source_str = "UNKNOWN";
        break;
    }

    const char* type_str = NULL;
    switch (type)
    {
    case GL_DEBUG_TYPE_ERROR_ARB:
        type_str = "ERROR";
        break;
    case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR_ARB:
        type_str = "DEPRECATED_BEHAVIOR";
        break;
    case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR_ARB:
        type_str = "UNDEFINED_BEHAVIOR";
        break;
    case GL_DEBUG_TYPE_PORTABILITY_ARB:
        type_str = "PORTABILITY";
        break;
    case GL_DEBUG_TYPE_PERFORMANCE_ARB:
        type_str = "PERFORMANCE";
        break;
    case GL_DEBUG_TYPE_OTHER_ARB:
        type_str = "OTHER";
        break;
    default:
        type_str = "UNKNOWN";
        break;
    }

    printf("[GL][%s][%s][%s]: %s", severity_str, source_str, type_str, message);
    if (length == 0 || message[length-1] != '\n')
        putchar('\n');
}

#endif

static bool Init()
{
    if (SDL_Init(SDL_INIT_EVENTS | SDL_INIT_VIDEO) != 0)
    {
        fprintf(stderr, "SDL_Init: %s\n", SDL_GetError());
        return false;
    }

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);

#if DEBUG_OPENGL
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_DEBUG_FLAG);
#endif

    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 0);
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 0);

    sdl_window = SDL_CreateWindow(WINDOW_TITLE, 0, 0, window_width, window_height,
                                  SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
    if (!sdl_window)
    {
        fprintf(stderr, "SDL_CreateWindow: %s\n", SDL_GetError());
        Shutdown();
        return false;
    }

    sdl_glcontext = SDL_GL_CreateContext(sdl_window);
    if (!sdl_glcontext)
    {
        fprintf(stderr, "SDL_GL_CreateContext: %s\n", SDL_GetError());
        Shutdown();
        return false;
    }

    GL_Init(&SDL_GL_GetProcAddress);

#if DEBUG_OPENGL

    printf("GL_VENDOR = %s\n", glGetString(GL_VENDOR));
    printf("GL_RENDERER = %s\n", glGetString(GL_RENDERER));
    printf("GL_VERSION = %s\n", glGetString(GL_VERSION));
    printf("GL_SHADING_LANGUAGE_VERSION = %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));

    if (g_GLInfo.ARB_debug_output)
    {
        glDebugMessageCallbackARB(&GLDebugCallback, NULL);

        glDebugMessageControlARB(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, NULL, GL_FALSE);
        glDebugMessageControlARB(GL_DONT_CARE, GL_DONT_CARE, GL_DEBUG_SEVERITY_LOW_ARB, 0, NULL, GL_TRUE);
        glDebugMessageControlARB(GL_DONT_CARE, GL_DONT_CARE, GL_DEBUG_SEVERITY_MEDIUM_ARB, 0, NULL, GL_TRUE);
        glDebugMessageControlARB(GL_DONT_CARE, GL_DONT_CARE, GL_DEBUG_SEVERITY_HIGH_ARB, 0, NULL, GL_TRUE);

        glEnable(GL_DEBUG_OUTPUT);
        glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
    }
    else
    {
        fprintf(stderr, "ARB_debug_output not supported\n");
    }

#endif

    // NOTE: This fixes flickering when resizing the window on Linux.
    SDL_GL_SetSwapInterval(0);

    return true;
}

static void Shutdown()
{
    if (sdl_glcontext)
    {
        SDL_GL_DeleteContext(sdl_glcontext);
        sdl_glcontext = NULL;
    }

    if (sdl_window)
    {
        SDL_DestroyWindow(sdl_window);
        sdl_window = NULL;
    }

    SDL_Quit();
}

static void InitImGui()
{
    ImGuiIO& io = ImGui::GetIO();

    io.DisplaySize.x = window_width;
    io.DisplaySize.y = window_height;

    io.KeyMap[ImGuiKey_Tab] = SDL_SCANCODE_TAB;
    io.KeyMap[ImGuiKey_LeftArrow] = SDL_SCANCODE_LEFT;
    io.KeyMap[ImGuiKey_RightArrow] = SDL_SCANCODE_RIGHT;
    io.KeyMap[ImGuiKey_UpArrow] = SDL_SCANCODE_UP;
    io.KeyMap[ImGuiKey_DownArrow] = SDL_SCANCODE_DOWN;
    io.KeyMap[ImGuiKey_PageUp] = SDL_SCANCODE_PAGEUP;
    io.KeyMap[ImGuiKey_PageDown] = SDL_SCANCODE_PAGEDOWN;
    io.KeyMap[ImGuiKey_Home] = SDL_SCANCODE_HOME;
    io.KeyMap[ImGuiKey_End] = SDL_SCANCODE_END;
    io.KeyMap[ImGuiKey_Delete] = SDL_SCANCODE_DELETE;
    io.KeyMap[ImGuiKey_Backspace] = SDL_SCANCODE_BACKSPACE;
    io.KeyMap[ImGuiKey_Enter] = SDL_SCANCODE_RETURN;
    io.KeyMap[ImGuiKey_Escape] = SDL_SCANCODE_ESCAPE;
    io.KeyMap[ImGuiKey_A] = SDL_SCANCODE_A;
    io.KeyMap[ImGuiKey_C] = SDL_SCANCODE_C;
    io.KeyMap[ImGuiKey_V] = SDL_SCANCODE_V;
    io.KeyMap[ImGuiKey_X] = SDL_SCANCODE_X;
    io.KeyMap[ImGuiKey_Y] = SDL_SCANCODE_Y;
    io.KeyMap[ImGuiKey_Z] = SDL_SCANCODE_Z;

    imgui_program = {
        "imgui",
        {
            {"shaders/imgui.vert", GL_VERTEX_SHADER},
            {"shaders/imgui.frag", GL_FRAGMENT_SHADER}
        }
    };
    GL_InitShaderProgram(&imgui_program);

    glUseProgram(imgui_program.id);

    glUniform1i(glGetUniformLocation(imgui_program.id, "u_Font"), 0);

    glUseProgram(0);

    glGenVertexArrays(1, &imgui_vao);
    glBindVertexArray(imgui_vao);

    glGenBuffers(1, &imgui_index_buffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, imgui_index_buffer);

    glGenBuffers(1, &imgui_vertex_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, imgui_vertex_buffer);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(ImDrawVert), (void*) offsetof(ImDrawVert, pos));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(ImDrawVert), (void*) offsetof(ImDrawVert, uv));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(ImDrawVert), (void*) offsetof(ImDrawVert, col));
    glEnableVertexAttribArray(2);

    unsigned char* font_data = NULL;
    int font_width = 0, font_height = 0;
    io.Fonts->GetTexDataAsAlpha8(&font_data, &font_width, &font_height);

    glGenTextures(1, &imgui_texture);
    glBindTexture(GL_TEXTURE_2D, imgui_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, font_width, font_height, 0, GL_RED, GL_UNSIGNED_BYTE, font_data);

    io.Fonts->TexID = (void*)(uintptr_t)imgui_texture;

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    const GLint font_swizzle[4] = {GL_RED, GL_RED, GL_RED, GL_RED};
    glTexParameteriv(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_RGBA, font_swizzle);
}

static void SendEventToImGui(const SDL_Event* event)
{
    ImGuiIO& io = ImGui::GetIO();

    switch (event->type)
    {
    case SDL_KEYDOWN:
    case SDL_KEYUP:
        switch (event->key.keysym.scancode)
        {
        case SDL_SCANCODE_LCTRL:
        case SDL_SCANCODE_RCTRL:
            io.KeyCtrl = (event->key.state == SDL_PRESSED);
            break;
        case SDL_SCANCODE_LSHIFT:
        case SDL_SCANCODE_RSHIFT:
            io.KeyShift = (event->key.state == SDL_PRESSED);
            break;
        case SDL_SCANCODE_LALT:
        case SDL_SCANCODE_RALT:
            io.KeyAlt = (event->key.state == SDL_PRESSED);
            break;
        case SDL_SCANCODE_LGUI:
        case SDL_SCANCODE_RGUI:
            io.KeySuper = (event->key.state == SDL_PRESSED);
            break;
        default:
            break;
        }
        io.KeysDown[event->key.keysym.scancode] = (event->key.state == SDL_PRESSED);
        break;
    case SDL_MOUSEBUTTONDOWN:
    case SDL_MOUSEBUTTONUP:
        switch (event->button.button)
        {
        case SDL_BUTTON_LEFT:
            io.MouseDown[0] = (event->button.state == SDL_PRESSED);
            break;
        case SDL_BUTTON_RIGHT:
            io.MouseDown[1] = (event->button.state == SDL_PRESSED);
            break;
        case SDL_BUTTON_MIDDLE:
            io.MouseDown[2] = (event->button.state == SDL_PRESSED);
            break;
        default:
            break;
        }
        break;
    case SDL_MOUSEWHEEL:
        if (event->wheel.y > 0)
            io.MouseWheel = 1;
        if (event->wheel.y < 0)
            io.MouseWheel = -1;
        break;
    case SDL_MOUSEMOTION:
        io.MousePos.x = event->motion.x;
        io.MousePos.y = event->motion.y;
        break;
    case SDL_TEXTINPUT:
        io.AddInputCharactersUTF8(event->text.text);
        break;
    case SDL_WINDOWEVENT:
        if (event->window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
        {
            io.DisplaySize.x = event->window.data1;
            io.DisplaySize.y = event->window.data2;
        }
        break;
    }
}

static void RenderImGui()
{
    ImGui::Render();

    glViewport(0, 0, window_width, window_height);

    glUseProgram(imgui_program.id);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, imgui_texture);

    glBindVertexArray(imgui_vao);
    glBindBuffer(GL_ARRAY_BUFFER, imgui_vertex_buffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, imgui_index_buffer);

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_SCISSOR_TEST);
    glDisable(GL_CULL_FACE);

    glUniformMatrix4fv(glGetUniformLocation(imgui_program.id, "u_ProjectionMatrix"), 1, GL_FALSE,
                       Matrix4::MakeOrtho(0, window_width, window_height, 0, -1, 1).data);

    ImDrawData* draw_data = ImGui::GetDrawData();
    for (int draw_list_index = 0; draw_list_index < draw_data->CmdListsCount; ++draw_list_index)
    {
        ImDrawList* draw_list = draw_data->CmdLists[draw_list_index];

        glBufferData(GL_ARRAY_BUFFER,
                     draw_list->VtxBuffer.Size * sizeof(ImDrawVert), draw_list->VtxBuffer.Data,
                     GL_STREAM_DRAW);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                     draw_list->IdxBuffer.Size * sizeof(ImDrawIdx), draw_list->IdxBuffer.Data,
                     GL_STREAM_DRAW);

        unsigned int first_index = 0;
        for (int draw_cmd_index = 0; draw_cmd_index < draw_list->CmdBuffer.Size; ++draw_cmd_index)
        {
            ImDrawCmd* draw_cmd = &draw_list->CmdBuffer.Data[draw_cmd_index];
            glBindTexture(GL_TEXTURE_2D, (GLuint)(uintptr_t)draw_cmd->TextureId);

            glScissor(draw_cmd->ClipRect.x,
                      window_height - draw_cmd->ClipRect.w,
                      draw_cmd->ClipRect.z - draw_cmd->ClipRect.x,
                      draw_cmd->ClipRect.w - draw_cmd->ClipRect.y);
            glDrawElements(GL_TRIANGLES, draw_cmd->ElemCount, GL_UNSIGNED_SHORT,
                           (void*) (first_index * sizeof(ImDrawIdx)));

            first_index += draw_cmd->ElemCount;
        }
    }

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_SCISSOR_TEST);
}

static void SubdivideSphereFace(int subdiv, Vector3* vertices, int* num_vertices, GLuint* indices, int* num_indices)
{
    if (subdiv == 0)
        return;

    GLuint i0 = indices[*num_indices - 3];
    GLuint i1 = indices[*num_indices - 2];
    GLuint i2 = indices[*num_indices - 1];
    Vector3 v0 = vertices[i0];
    Vector3 v1 = vertices[i1];
    Vector3 v2 = vertices[i2];

    GLuint i01 = (*num_vertices)++;
    GLuint i12 = (*num_vertices)++;
    GLuint i20 = (*num_vertices)++;
    Vector3 v01 = Math::Normalize(v0 + v1);
    Vector3 v12 = Math::Normalize(v1 + v2);
    Vector3 v20 = Math::Normalize(v2 + v0);

    vertices[i01] = v01;
    vertices[i12] = v12;
    vertices[i20] = v20;

    indices[*num_indices - 3] = i01;
    indices[*num_indices - 2] = i12;
    indices[*num_indices - 1] = i20;
    SubdivideSphereFace(subdiv-1, vertices, num_vertices, indices, num_indices);

    *num_indices += 3;
    indices[*num_indices - 3] = i0;
    indices[*num_indices - 2] = i01;
    indices[*num_indices - 1] = i20;
    SubdivideSphereFace(subdiv-1, vertices, num_vertices, indices, num_indices);

    *num_indices += 3;
    indices[*num_indices - 3] = i1;
    indices[*num_indices - 2] = i12;
    indices[*num_indices - 1] = i01;
    SubdivideSphereFace(subdiv-1, vertices, num_vertices, indices, num_indices);

    *num_indices += 3;
    indices[*num_indices - 3] = i2;
    indices[*num_indices - 2] = i20;
    indices[*num_indices - 1] = i12;
    SubdivideSphereFace(subdiv-1, vertices, num_vertices, indices, num_indices);
}

static void GenerateIcoSphereMesh(Mesh* mesh, int subdiv)
{
    int total_num_vertices = 12;
    int total_num_faces = 20;
    for (int i = 0; i < subdiv; ++i)
    {
        total_num_vertices += total_num_faces * 3;
        total_num_faces *= 4;
    }

    glGenVertexArrays(1, &mesh->vao);
    glBindVertexArray(mesh->vao);

    glGenBuffers(1, &mesh->vertex_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, mesh->vertex_buffer);
    glBufferData(GL_ARRAY_BUFFER, total_num_vertices * sizeof(Vector3), NULL, GL_STATIC_DRAW);
    Vector3* vertices = (Vector3*) glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(0);

    glGenBuffers(1, &mesh->index_buffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->index_buffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, total_num_faces * 3 * sizeof(GLuint), NULL, GL_STATIC_DRAW);
    GLuint* indices = (GLuint*) glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);

    static const Vector3 base_vertices[12] = {
        { 0.000000,  0.000000, -1.000000},
        { 0.723607, -0.525725, -0.447220},
        {-0.276388, -0.850649, -0.447220},
        {-0.894426,  0.000000, -0.447216},
        {-0.276388,  0.850649, -0.447220},
        { 0.723607,  0.525725, -0.447220},
        { 0.276388, -0.850649,  0.447220},
        {-0.723607, -0.525725,  0.447220},
        {-0.723607,  0.525725,  0.447220},
        { 0.276388,  0.850649,  0.447220},
        { 0.894426,  0.000000,  0.447216},
        { 0.000000,  0.000000,  1.000000},
    };

    int num_vertices = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(base_vertices); ++i)
    {
        vertices[i] = base_vertices[i];
        ++num_vertices;
    }

    static const GLuint base_indices[20][3] = {
        { 0,  1,  2},
        { 1,  0,  5},
        { 0,  2,  3},
        { 0,  3,  4},
        { 0,  4,  5},
        { 1,  5, 10},
        { 2,  1,  6},
        { 3,  2,  7},
        { 4,  3,  8},
        { 5,  4,  9},
        { 1, 10,  6},
        { 2,  6,  7},
        { 3,  7,  8},
        { 4,  8,  9},
        { 5,  9, 10},
        { 6, 10, 11},
        { 7,  6, 11},
        { 8,  7, 11},
        { 9,  8, 11},
        {10,  9, 11},
    };

    int num_indices = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(base_indices); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            indices[num_indices] = base_indices[i][j];
            ++num_indices;
        }

        SubdivideSphereFace(subdiv, vertices, &num_vertices, indices, &num_indices);
    }

    assert(num_vertices == total_num_vertices);
    assert(num_indices == total_num_faces * 3);

    mesh->num_vertices = num_vertices;
    mesh->num_indices = num_indices;

    glUnmapBuffer(GL_ARRAY_BUFFER);
    glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
}

static void InitEnvMapTool(EnvMapTool* tool)
{
    glClearColor(0, 0, 0, 0);

    glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);

    tool->camera.fovy = Math::PI / 3;
    tool->camera.aspect = (float) (window_width - window_width / 4) / (float) window_height;
    tool->camera.znear = 0.1f;
    tool->camera.zfar = 100.0f;
    tool->camera.transform.translation = Vector3(0, 0, 3);
    tool->camera.transform.rotation = Quaternion::MakeIdentity();
    tool->camera.transform.scale = 1;

    tool->filter.type = FILTER_TYPE_PHONG_COSINE;
    tool->filter.face_size = 256;
    tool->filter.max_mip_level = 4;
    tool->filter.sh_approximation = true;
    tool->filter.sh_max_degree = 4;
    tool->filter.sh_window_strength = 0;
    tool->filter.num_samples = 10000;
    tool->filter.specular_exponent = 8192;
    tool->filter.specular_exponent_divisor = 8;
    tool->filter.roughness = 0;
    tool->filter.roughness_increment = 0.2;

    tool->display.env_map = 0;
    tool->display.tex_coord = 0;
    tool->display.lod = 0;
    tool->display.exposure = 0;

    tool->draw_mesh_program = {
        "draw_mesh_program",
        {
            {"shaders/draw_mesh.vert", GL_VERTEX_SHADER},
            {"shaders/draw_mesh.frag", GL_FRAGMENT_SHADER}
        }
    };
    GL_InitShaderProgram(&tool->draw_mesh_program);

    tool->filter_cosine_program = {
        "filter_cosine",
        {
            {"shaders/filter.vert", GL_VERTEX_SHADER},
            {"shaders/filter_cosine.frag", GL_FRAGMENT_SHADER}
        }
    };
    GL_InitShaderProgram(&tool->filter_cosine_program);

    tool->filter_phong_cosine_program = {
        "filter_phong_cosine",
        {
            {"shaders/filter.vert", GL_VERTEX_SHADER},
            {"shaders/filter_phong_cosine.frag", GL_FRAGMENT_SHADER}
        }
    };
    GL_InitShaderProgram(&tool->filter_phong_cosine_program);

    tool->filter_beckmann_cosine_program = {
        "filter_beckmann_cosine",
        {
            {"shaders/filter.vert", GL_VERTEX_SHADER},
            {"shaders/filter_beckmann_cosine.frag", GL_FRAGMENT_SHADER}
        }
    };
    GL_InitShaderProgram(&tool->filter_beckmann_cosine_program);

    tool->filter_ggx_cosine_program = {
        "filter_ggx_cosine",
        {
            {"shaders/filter.vert", GL_VERTEX_SHADER},
            {"shaders/filter_ggx_cosine.frag", GL_FRAGMENT_SHADER}
        }
    };
    GL_InitShaderProgram(&tool->filter_ggx_cosine_program);

    glGenVertexArrays(1, &tool->filter_vao);
    glBindVertexArray(tool->filter_vao);
    glBindVertexArray(0);

    glGenFramebuffers(1, &tool->filter_fbo);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, tool->filter_fbo);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    const float default_face_colors[6][3] = {
        {1, 0, 0}, {0, 1, 1},
        {0, 1, 0}, {1, 0, 1},
        {0, 0, 1}, {1, 1, 0},
    };

    tool->source_env_map = AllocateCubeMap(1, 256);
    tool->filtered_env_map = AllocateCubeMap(1, 256);

    FillCubeMapWithFaceColors(tool->source_env_map, default_face_colors);
    FillCubeMapWithFaceColors(tool->filtered_env_map, default_face_colors);

    UploadCubeMapToGPU(tool->source_env_map);
    UploadCubeMapToGPU(tool->filtered_env_map);

    tool->layout = CUBEMAP_LAYOUT_HORIZONTAL_CROSS;

    GenerateIcoSphereMesh(&tool->sphere_mesh, 3);

    tool->filter_task.filter_params = tool->filter;
    tool->filter_task.is_running = false;
    tool->filter_task.num_samples_completed = 0;
    tool->filter_task.num_samples_per_frame = 1;

    glGenQueries(1, &tool->filter_task.frame_time_query);
}

#define SH_MAX_DEGREE 9

static void FilterCubeMap_Cosine_SH(EnvMapTool* tool)
{
    int face_size = tool->filter.face_size;
    int sh_max_degree = tool->filter.sh_max_degree;
    float sh_window_strength = tool->filter.sh_window_strength;

    if (tool->filtered_env_map) FreeCubeMap(tool->filtered_env_map);
    tool->filtered_env_map = AllocateCubeMap(1, face_size);

    float shcoeff_r[sh_Count(SH_MAX_DEGREE)], shcoeff_g[sh_Count(SH_MAX_DEGREE)], shcoeff_b[sh_Count(SH_MAX_DEGREE)];

    DownloadCubeMapFromGPU(tool->source_env_map);
    ProjectCubeMapOntoSH(tool->source_env_map, sh_max_degree, shcoeff_r, shcoeff_g, shcoeff_b);

    static const float CosineLobeZH[] = {
        0.886226925, 1.023326708, 0.495415912, 0.0, -0.110778366, 0.0, 0.049927135, 0.0, -0.028546931, 0.0
    };

    // NOTE: See "An Efficient Representation for Irradiance Environment Maps" by Ramamoorthi and Hanrahan.

    for (int l = 0; l <= sh_max_degree; ++l)
    {
        float f = sqrt(4*Math::PI / (2*l+1)) * CosineLobeZH[l];

        // NOTE: Blend between original and windowed SH coefficients based on sh_window_strength.
        float hann = (1 + cos(Math::PI*l / (sh_max_degree+1))) / 2;
        f *= Math::Lerp(1, hann, sh_window_strength);

        for (int m = -l; m <= l; ++m)
        {
            shcoeff_r[sh_Index(l, m)] *= f;
            shcoeff_g[sh_Index(l, m)] *= f;
            shcoeff_b[sh_Index(l, m)] *= f;
        }
    }

    ReconstructCubeMapFromSH(tool->filtered_env_map, sh_max_degree, shcoeff_r, shcoeff_g, shcoeff_b);
    UploadCubeMapToGPU(tool->filtered_env_map);
}

static void BeginFilterCubeMap(EnvMapTool* tool, int face_size, int max_mip_level)
{
    if (tool->filtered_env_map) FreeCubeMap(tool->filtered_env_map);
    tool->filtered_env_map = AllocateCubeMap(max_mip_level + 1, face_size);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_CUBE_MAP, tool->filtered_env_map->gl_id);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAX_LEVEL, max_mip_level);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, tool->filter_fbo);

    for (int mip = 0; mip <= max_mip_level; mip++)
    {
        glViewport(0, 0, face_size, face_size);

        for (int face = 0; face < 6; ++face)
        {
            glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + face, mip, GL_RGBA32F, face_size, face_size, 0,
                         GL_RGB, GL_FLOAT, NULL);

            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                   GL_TEXTURE_CUBE_MAP_POSITIVE_X + face,
                                   tool->filtered_env_map->gl_id, mip);

            glClear(GL_COLOR_BUFFER_BIT);
        }

        face_size /= 2;
        if (face_size == 0)
            break;
    }

    tool->filter_task.filter_params = tool->filter;
    tool->filter_task.is_running = true;
    tool->filter_task.num_samples_completed = 0;
    tool->filter_task.num_samples_per_frame = 1;

    // NOTE: Force OpenGL to create a query, so that subsequent calls to glGetQueryObject won't fail.
    glBeginQuery(GL_TIME_ELAPSED, tool->filter_task.frame_time_query);
    glEndQuery(GL_TIME_ELAPSED);
}

static void UpdateFilterCubeMap_Cosine(EnvMapTool* tool)
{
    const FilterSettings params = tool->filter_task.filter_params;

    int face_size = params.face_size;
    int num_samples = params.num_samples;

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);

    glUseProgram(tool->filter_cosine_program.id);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_CUBE_MAP, tool->source_env_map->gl_id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_CUBE_MAP, tool->filtered_env_map->gl_id);

    glUniform1i(glGetUniformLocation(tool->filter_cosine_program.id, "u_EnvMap"), 0);

    int start_sample = tool->filter_task.num_samples_completed;
    int end_sample = start_sample + tool->filter_task.num_samples_per_frame;
    if (end_sample > num_samples) end_sample = num_samples;

    glUniform1i(glGetUniformLocation(tool->filter_cosine_program.id, "u_NumSamples"), num_samples);
    glUniform1i(glGetUniformLocation(tool->filter_cosine_program.id, "u_StartSample"), start_sample);
    glUniform1i(glGetUniformLocation(tool->filter_cosine_program.id, "u_EndSample"), end_sample);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, tool->filter_fbo);

    glViewport(0, 0, face_size, face_size);

    for (int face = 0; face < 6; ++face)
    {
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                               GL_TEXTURE_CUBE_MAP_POSITIVE_X + face,
                               tool->filtered_env_map->gl_id, 0);

        glUniformMatrix3fv(glGetUniformLocation(tool->filter_cosine_program.id, "u_FaceBasis"), 1, GL_FALSE,
                           (const float*) (g_CubeMapFaces + face));

        glBindVertexArray(tool->filter_vao);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    }

    tool->filter_task.num_samples_completed = end_sample;
    if (tool->filter_task.num_samples_completed == num_samples)
        tool->filter_task.is_running = false;
}

static void UpdateFilterCubeMap_PhongCosine(EnvMapTool* tool)
{
    const FilterSettings params = tool->filter_task.filter_params;

    int face_size = params.face_size;
    int max_mip_level = params.max_mip_level;
    int num_samples = params.num_samples;
    float specular_exponent = params.specular_exponent;
    float specular_exponent_divisor = params.specular_exponent_divisor;

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);

    glUseProgram(tool->filter_phong_cosine_program.id);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_CUBE_MAP, tool->source_env_map->gl_id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_CUBE_MAP, tool->filtered_env_map->gl_id);

    glUniform1i(glGetUniformLocation(tool->filter_phong_cosine_program.id, "u_EnvMap"), 0);

    int start_sample = tool->filter_task.num_samples_completed;
    int end_sample = start_sample + tool->filter_task.num_samples_per_frame;
    if (end_sample > num_samples) end_sample = num_samples;

    glUniform1i(glGetUniformLocation(tool->filter_phong_cosine_program.id, "u_NumSamples"), num_samples);
    glUniform1i(glGetUniformLocation(tool->filter_phong_cosine_program.id, "u_StartSample"), start_sample);
    glUniform1i(glGetUniformLocation(tool->filter_phong_cosine_program.id, "u_EndSample"), end_sample);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, tool->filter_fbo);

    for (int mip = 0; mip <= max_mip_level; mip++)
    {
        glViewport(0, 0, face_size, face_size);

        glUniform1f(glGetUniformLocation(tool->filter_phong_cosine_program.id, "u_SpecularExponent"), specular_exponent);

        for (int face = 0; face < 6; ++face)
        {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                   GL_TEXTURE_CUBE_MAP_POSITIVE_X + face,
                                   tool->filtered_env_map->gl_id, mip);

            glUniformMatrix3fv(glGetUniformLocation(tool->filter_phong_cosine_program.id, "u_FaceBasis"), 1, GL_FALSE,
                               (const float*) (g_CubeMapFaces + face));

            glBindVertexArray(tool->filter_vao);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        }

        face_size /= 2;
        if (face_size == 0)
            break;

        specular_exponent /= specular_exponent_divisor;
    }

    tool->filter_task.num_samples_completed = end_sample;
    if (tool->filter_task.num_samples_completed == num_samples)
        tool->filter_task.is_running = false;
}

static void UpdateFilterCubeMap_BeckmannCosine(EnvMapTool* tool)
{
    const FilterSettings params = tool->filter_task.filter_params;

    int face_size = params.face_size;
    int max_mip_level = params.max_mip_level;
    int num_samples = params.num_samples;
    float roughness = params.roughness;
    float roughness_increment = params.roughness_increment;

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);

    glUseProgram(tool->filter_beckmann_cosine_program.id);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_CUBE_MAP, tool->source_env_map->gl_id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_CUBE_MAP, tool->filtered_env_map->gl_id);

    glUniform1i(glGetUniformLocation(tool->filter_beckmann_cosine_program.id, "u_EnvMap"), 0);

    int start_sample = tool->filter_task.num_samples_completed;
    int end_sample = start_sample + tool->filter_task.num_samples_per_frame;
    if (end_sample > num_samples) end_sample = num_samples;

    glUniform1i(glGetUniformLocation(tool->filter_beckmann_cosine_program.id, "u_NumSamples"), num_samples);
    glUniform1i(glGetUniformLocation(tool->filter_beckmann_cosine_program.id, "u_StartSample"), start_sample);
    glUniform1i(glGetUniformLocation(tool->filter_beckmann_cosine_program.id, "u_EndSample"), end_sample);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, tool->filter_fbo);

    for (int mip = 0; mip <= max_mip_level; mip++)
    {
        glViewport(0, 0, face_size, face_size);

        roughness = Math::Clamp(roughness, 0, 1);
        glUniform1f(glGetUniformLocation(tool->filter_beckmann_cosine_program.id, "u_Roughness"), roughness*roughness);

        for (int face = 0; face < 6; ++face)
        {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                   GL_TEXTURE_CUBE_MAP_POSITIVE_X + face,
                                   tool->filtered_env_map->gl_id, mip);

            glUniformMatrix3fv(glGetUniformLocation(tool->filter_beckmann_cosine_program.id, "u_FaceBasis"), 1, GL_FALSE,
                               (const float*) (g_CubeMapFaces + face));

            glBindVertexArray(tool->filter_vao);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        }

        face_size /= 2;
        if (face_size == 0)
            break;

        roughness += roughness_increment;
    }

    tool->filter_task.num_samples_completed = end_sample;
    if (tool->filter_task.num_samples_completed == num_samples)
        tool->filter_task.is_running = false;
}

static void UpdateFilterCubeMap_GGXCosine(EnvMapTool* tool)
{
    const FilterSettings params = tool->filter_task.filter_params;

    int face_size = params.face_size;
    int max_mip_level = params.max_mip_level;
    int num_samples = params.num_samples;
    float roughness = params.roughness;
    float roughness_increment = params.roughness_increment;

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);

    glUseProgram(tool->filter_ggx_cosine_program.id);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_CUBE_MAP, tool->source_env_map->gl_id);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_CUBE_MAP, tool->filtered_env_map->gl_id);

    glUniform1i(glGetUniformLocation(tool->filter_ggx_cosine_program.id, "u_EnvMap"), 0);

    int start_sample = tool->filter_task.num_samples_completed;
    int end_sample = start_sample + tool->filter_task.num_samples_per_frame;
    if (end_sample > num_samples) end_sample = num_samples;

    glUniform1i(glGetUniformLocation(tool->filter_ggx_cosine_program.id, "u_NumSamples"), num_samples);
    glUniform1i(glGetUniformLocation(tool->filter_ggx_cosine_program.id, "u_StartSample"), start_sample);
    glUniform1i(glGetUniformLocation(tool->filter_ggx_cosine_program.id, "u_EndSample"), end_sample);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, tool->filter_fbo);

    for (int mip = 0; mip <= max_mip_level; mip++)
    {
        glViewport(0, 0, face_size, face_size);

        roughness = Math::Clamp(roughness, 0, 1);
        glUniform1f(glGetUniformLocation(tool->filter_ggx_cosine_program.id, "u_Roughness"), roughness*roughness);

        for (int face = 0; face < 6; ++face)
        {
            glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                   GL_TEXTURE_CUBE_MAP_POSITIVE_X + face,
                                   tool->filtered_env_map->gl_id, mip);

            glUniformMatrix3fv(glGetUniformLocation(tool->filter_ggx_cosine_program.id, "u_FaceBasis"), 1, GL_FALSE,
                               (const float*) (g_CubeMapFaces + face));

            glBindVertexArray(tool->filter_vao);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        }

        face_size /= 2;
        if (face_size == 0)
            break;

        roughness += roughness_increment;
    }

    tool->filter_task.num_samples_completed = end_sample;
    if (tool->filter_task.num_samples_completed == num_samples)
        tool->filter_task.is_running = false;
}

static void UpdateEnvMapTool(EnvMapTool* tool, int buttons, int x, int y, int dwheel, int dx, int dy)
{
    if (tool->filter_task.is_running)
    {
        FilterTaskProgress* task = &tool->filter_task;

        GLint query_result_available = 0;
        glGetQueryObjectiv(task->frame_time_query, GL_QUERY_RESULT_AVAILABLE, &query_result_available);

        if (query_result_available)
        {
            GLuint64 elapsed_time = 0;
            glGetQueryObjectui64v(task->frame_time_query, GL_QUERY_RESULT, &elapsed_time);

            float f = (elapsed_time > 1.0e6) ? Math::Clamp(FILTER_MS_PER_FRAME / (elapsed_time * 1.0e-6), 0.5, 2) : 2;
            task->num_samples_per_frame *= f;
            if (task->num_samples_per_frame < 1)
                task->num_samples_per_frame = 1;
        }



        if (query_result_available) glBeginQuery(GL_TIME_ELAPSED, tool->filter_task.frame_time_query);

        switch (tool->filter_task.filter_params.type)
        {
        case FILTER_TYPE_COSINE:
            if (tool->filter_task.filter_params.sh_approximation)
                INVALID_CODE_PATH; // NOTE: SH approximation is done synchronously on the CPU side.
            else
                UpdateFilterCubeMap_Cosine(tool);
            break;
        case FILTER_TYPE_PHONG_COSINE:
            UpdateFilterCubeMap_PhongCosine(tool);
            break;
        case FILTER_TYPE_BECKMANN_COSINE:
            UpdateFilterCubeMap_BeckmannCosine(tool);
            break;
        case FILTER_TYPE_GGX_COSINE:
            UpdateFilterCubeMap_GGXCosine(tool);
            break;
        default:
            INVALID_CODE_PATH;
        }

        if (query_result_available) glEndQuery(GL_TIME_ELAPSED);
    }

    // draw main panel

    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0);

    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(ImVec2(window_width / 4, window_height));
    unsigned long window_flags = ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse;
    if (ImGui::Begin("EnvMapTool", NULL, window_flags))
    {
        if (ImGui::CollapsingHeader("Import", ImGuiTreeNodeFlags_DefaultOpen))
        {
            static char filename[1024] = {};
            ImGui::InputText("filename##import", filename, sizeof(filename));

            if (ImGui::Button("Load Cube Map"))
            {
                CubeMap* cubemap = LoadCubeMap(filename, &tool->layout);
                if (cubemap)
                {
                    if (tool->source_env_map) FreeCubeMap(tool->source_env_map);
                    tool->source_env_map = cubemap;
                    tool->display.env_map = 0;
                }

                tool->filter_task.is_running = false;
            }
        }

        if (ImGui::CollapsingHeader("Export", ImGuiTreeNodeFlags_DefaultOpen))
        {
            static char filename[1024] = {};
            ImGui::InputText("filename##export", filename, sizeof(filename));

            static const char* layout_items[] = {
                "Auto",
                "Horizontal Cross",
                "Vertical Cross",
                "Horizontal Strip",
                "Vertical Strip",
            };

            static int layout = 0;
            ImGui::Combo("layout##export", &layout, layout_items, ARRAY_SIZE(layout_items));

            if (ImGui::Button("Save Cube Map"))
            {
                DownloadCubeMapFromGPU(tool->filtered_env_map);

                if (layout == 0)
                    SaveCubeMap(tool->filtered_env_map, filename, tool->layout);
                else
                    SaveCubeMap(tool->filtered_env_map, filename, (CubeMapLayout) layout);

                tool->filter_task.is_running = false;
            }
        }

        if (ImGui::CollapsingHeader("Filter", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::PushItemWidth(ImGui::GetWindowWidth() / 2);

            static const char* filter_type_items[] = {
                "Cosine",
                "Phong NDF with Cosine",
                "Beckmann NDF with Cosine",
                "GGX NDF with Cosine"
            };

            int filter_type = tool->filter.type;
            ImGui::Combo("type", &filter_type, filter_type_items, ARRAY_SIZE(filter_type_items));
            tool->filter.type = (FilterType) filter_type;

            ImGui::InputInt("face_size", &tool->filter.face_size);
            if (tool->filter.face_size < 1) tool->filter.face_size = 1;

            if (tool->filter.type != FILTER_TYPE_COSINE)
            {
                ImGui::InputInt("max_mip_level", &tool->filter.max_mip_level);
                if (tool->filter.max_mip_level < 0) tool->filter.max_mip_level = 0;
            }

            if (tool->filter.type == FILTER_TYPE_COSINE)
            {
                ImGui::Checkbox("sh_approximation", &tool->filter.sh_approximation);

                if (tool->filter.sh_approximation)
                {
                    ImGui::InputInt("sh_max_degree", &tool->filter.sh_max_degree);
                    if (tool->filter.sh_max_degree < 0) tool->filter.sh_max_degree = 0;
                    if (tool->filter.sh_max_degree > SH_MAX_DEGREE) tool->filter.sh_max_degree = SH_MAX_DEGREE;
                    ImGui::DragFloat("sh_window_strength", &tool->filter.sh_window_strength, 0.01f, 0, 1);
                    if (tool->filter.sh_window_strength < 0) tool->filter.sh_window_strength = 0;
                    if (tool->filter.sh_window_strength > 1) tool->filter.sh_window_strength = 1;
                }
            }

            if (tool->filter.type != FILTER_TYPE_COSINE || !tool->filter.sh_approximation)
            {
                ImGui::InputInt("num_samples", &tool->filter.num_samples);
                if (tool->filter.num_samples < 1) tool->filter.num_samples = 1;
            }

            if (tool->filter.type == FILTER_TYPE_PHONG_COSINE)
            {
                ImGui::InputFloat("specular_exponent", &tool->filter.specular_exponent);
                if (tool->filter.specular_exponent < 0) tool->filter.specular_exponent = 0;
                ImGui::InputFloat("specular_exponent_divisor", &tool->filter.specular_exponent_divisor);
                if (tool->filter.specular_exponent_divisor < 1) tool->filter.specular_exponent_divisor = 1;
            }

            if (tool->filter.type == FILTER_TYPE_BECKMANN_COSINE || tool->filter.type == FILTER_TYPE_GGX_COSINE)
            {
                ImGui::DragFloat("roughness", &tool->filter.roughness, 0.01f, 0, 1);
                if (tool->filter.roughness < 0) tool->filter.roughness = 0;
                if (tool->filter.roughness > 1) tool->filter.roughness = 1;
                ImGui::DragFloat("roughness_increment", &tool->filter.roughness_increment, 0.01f, 0, 1);
                if (tool->filter.roughness_increment < 0) tool->filter.roughness_increment = 0;
                if (tool->filter.roughness_increment > 1) tool->filter.roughness_increment = 1;
            }

            ImGui::PopItemWidth();

            if (ImGui::Button("Filter Cube Map"))
            {
                tool->filter_task.is_running = false;

                switch (tool->filter.type)
                {
                case FILTER_TYPE_COSINE:
                    if (tool->filter.sh_approximation)
                        FilterCubeMap_Cosine_SH(tool);
                    else
                        BeginFilterCubeMap(tool, tool->filter.face_size, 0);
                    break;
                case FILTER_TYPE_PHONG_COSINE:
                    BeginFilterCubeMap(tool, tool->filter.face_size, tool->filter.max_mip_level);
                    break;
                case FILTER_TYPE_BECKMANN_COSINE:
                    BeginFilterCubeMap(tool, tool->filter.face_size, tool->filter.max_mip_level);
                    break;
                case FILTER_TYPE_GGX_COSINE:
                    BeginFilterCubeMap(tool, tool->filter.face_size, tool->filter.max_mip_level);
                    break;
                default:
                    INVALID_CODE_PATH;
                }

                tool->display.env_map = 1;
            }
        }

        if (ImGui::CollapsingHeader("Display", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::PushItemWidth(ImGui::GetWindowWidth() / 2);

            static const char* display_env_map_items[] = {"Source", "Filtered"};
            ImGui::Combo("env_map", &tool->display.env_map,
                         display_env_map_items, ARRAY_SIZE(display_env_map_items));

            static const char* display_tex_coord_items[] = {"Normal Vector", "Reflection Vector"};
            ImGui::Combo("tex_coord", &tool->display.tex_coord,
                         display_tex_coord_items, ARRAY_SIZE(display_tex_coord_items));

            ImGui::DragFloat("lod", &tool->display.lod, 0.1f, 0.0f, 100.0f);

            ImGui::DragFloat("exposure", &tool->display.exposure, 0.1f, -100.0f, 100.0f);

            ImGui::PopItemWidth();
        }

        // draw filtering progress overlay

        {
            char text[32] = {};
            uint32_t color = 0;
            if (tool->filter_task.is_running)
            {
                int progress = 100.0 * tool->filter_task.num_samples_completed / tool->filter_task.filter_params.num_samples;
                snprintf(text, sizeof(text), "%d %%", progress);
                color = 0xFF00FFFF;
            }
            else
            {
                snprintf(text, sizeof(text), "DONE");
                color = 0xFF00FF00;
            }

            float text_width = ImGui::CalcTextSize(text).x;

            const float padding_x = 10;
            const float padding_y = 5;

            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->PushClipRectFullScreen();
            draw_list->AddText(ImVec2(window_width - (text_width + padding_x), padding_y), color, text);
            draw_list->PopClipRect();
        }
    }
    ImGui::End();

    ImGui::PopStyleVar();

    // update camera

    tool->camera.aspect = (float) (window_width - window_width / 4) / (float) window_height;

    if (!ImGui::GetIO().WantCaptureMouse)
    {
        Matrix4 matrix = tool->camera.transform.GetMatrix();
        Vector3 local_x_axis = Vector3(matrix._11, matrix._21, matrix._31);
        Vector3 local_y_axis = Vector3(matrix._12, matrix._22, matrix._32);
        Vector3 local_z_axis = Vector3(matrix._13, matrix._23, matrix._33);

        if (buttons & BIT(0))
        {
            Transform tmp = Transform::MakeIdentity();
            tmp.rotation = AxisAngleToQuaternion(local_x_axis, -dy/256.0f);
            tool->camera.transform = tmp * tool->camera.transform;
            tmp.rotation = AxisAngleToQuaternion(local_y_axis, -dx/256.0f);
            tool->camera.transform = tmp * tool->camera.transform;
        }

        if (buttons & BIT(1))
        {
            int cx = window_width * 5 / 8; // window_width / 4 + (window_width - window_width / 4) / 2;
            int cy = window_height / 2;

            float rz = (-dx/256.0f) * (cy<y ? 1 : -1) + (-dy/256.0f) * (cx<x ? -1 : 1); // black magic

            Transform tmp = Transform::MakeIdentity();
            tmp.rotation = AxisAngleToQuaternion(local_z_axis, rz);
            tool->camera.transform = tmp * tool->camera.transform;
        }

        tool->camera.transform.translation = tool->camera.transform.translation * exp(-dwheel/10.0f);

        const float min_dist = 2;
        const float max_dist = 20;
        float dist_to_origin = Math::Length(tool->camera.transform.translation);
        float new_dist_to_origin = Math::Clamp(dist_to_origin, min_dist, max_dist);
        tool->camera.transform.translation = Math::Normalize(tool->camera.transform.translation) * new_dist_to_origin;

    }

    // draw sphere

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    glViewport(window_width / 4, 0, window_width - window_width / 4, window_height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glDisable(GL_BLEND);

    glUseProgram(tool->draw_mesh_program.id);

    glUniform3fv(glGetUniformLocation(tool->draw_mesh_program.id, "u_CameraPosition"), 1,
                 tool->camera.transform.translation.data);

    glUniform1i(glGetUniformLocation(tool->draw_mesh_program.id, "u_DisplayMode"), tool->display.tex_coord);

    glActiveTexture(GL_TEXTURE0);
    glUniform1i(glGetUniformLocation(tool->draw_mesh_program.id, "u_EnvMap"), 0);
    if (tool->display.env_map == 0)
        glBindTexture(GL_TEXTURE_CUBE_MAP, tool->source_env_map->gl_id);
    else
        glBindTexture(GL_TEXTURE_CUBE_MAP, tool->filtered_env_map->gl_id);

    glUniform1f(glGetUniformLocation(tool->draw_mesh_program.id, "u_MipLevel"), tool->display.lod);
    glUniform1f(glGetUniformLocation(tool->draw_mesh_program.id, "u_Exposure"), tool->display.exposure);

    Matrix4 projection_matrix =
        Matrix4::MakePerspective(tool->camera.fovy, tool->camera.aspect, tool->camera.znear, tool->camera.zfar);
    Matrix4 view_matrix =
        MakeWorldToLocalMatrix(tool->camera.transform.translation, tool->camera.transform.rotation, 1);

    glUniformMatrix4fv(glGetUniformLocation(tool->draw_mesh_program.id, "u_WorldToClipMatrix"), 1, GL_FALSE,
                       (projection_matrix * view_matrix).data);
    glUniformMatrix4fv(glGetUniformLocation(tool->draw_mesh_program.id, "u_ObjectToWorldMatrix"), 1, GL_FALSE,
                       Matrix4::MakeIdentity().data);

    glBindVertexArray(tool->sphere_mesh.vao);
    glDrawElements(GL_TRIANGLES, tool->sphere_mesh.num_indices, GL_UNSIGNED_INT, 0);
}
