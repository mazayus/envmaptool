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

#include "common.h"

void DebugPrint(const char* format, ...)
{
    va_list va;
    va_start(va, format);
    DebugPrintVA(format, va);
    va_end(va);
}

void DebugPrintVA(const char* format, va_list va)
{
    vprintf(format, va);
}

long GetFileSize(const char* filename)
{
    FILE* fp = fopen(filename, "rb");
    if (!fp)
    {
        fprintf(stderr, "GetFileSize(\"%s\"): %s\n", filename, strerror(errno));
        return -1;
    }

    if (fseek(fp, 0, SEEK_END) == -1)
    {
        fprintf(stderr, "GetFileSize(\"%s\"): %s\n", filename, strerror(errno));
        return -1;
    }

    long size = ftell(fp);
    if (size == -1)
    {
        fprintf(stderr, "GetFileSize(\"%s\"): %s\n", filename, strerror(errno));
        return -1;
    }

    fclose(fp);

    return size;
}

long GetFileContents(const char* filename, void* buf, long bufsize)
{
    FILE* fp = fopen(filename, "rb");
    if (!fp)
    {
        fprintf(stderr, "GetFileContents(\"%s\"): %s\n", filename, strerror(errno));
        return -1;
    }

    size_t numbytes = fread(buf, 1, bufsize, fp);
    if ((numbytes != (size_t) bufsize) && ferror(fp))
    {
        fprintf(stderr, "GetFileContents(\"%s\"): %s\n", filename, strerror(errno));
        fclose(fp);
        return numbytes;
    }

    fclose(fp);

    return numbytes;
}

#define SCRATCH_SIZE SIZE_MB(16)

static char     scratch_memory[SCRATCH_SIZE];
static size_t   scratch_allocated_size = 0;
static size_t   scratch_total_size = SCRATCH_SIZE;

void* ScratchAlloc(size_t size)
{
    if (scratch_allocated_size + size > scratch_total_size)
    {
        fprintf(stderr, "ScratchAlloc: out of memory\n");
        return NULL;
    }

    if (scratch_allocated_size + size < scratch_allocated_size)
    {
        fprintf(stderr, "ScratchAlloc: memory overflow\n");
        return NULL;
    }

    void* ptr = (char*) scratch_memory + scratch_allocated_size;

    scratch_allocated_size += size;

    return ptr;
}

void ScratchFreeTo(void* ptr)
{
    size_t new_allocated_size = (char*) ptr - (char*) scratch_memory;

    if (new_allocated_size > scratch_allocated_size)
    {
        fprintf(stderr, "ScratchFreeTo: invalid pointer\n");
        return;
    }

    scratch_allocated_size = new_allocated_size;
}

void ScratchClear()
{
    scratch_allocated_size = 0;
}
