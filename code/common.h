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

#ifndef COMMON_H
#define COMMON_H

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ARRAY_SIZE(array) (sizeof(array) / sizeof((array)[0]))

#define FOURCC(a, b, c, d) ((a) | ((b) << 8) | ((c) << 16) | ((d) << 24))

#define BIT(n) (1 << (n))

#define SIZE_KB(n) ((size_t)(n) << 10)
#define SIZE_MB(n) ((size_t)(n) << 20)
#define SIZE_GB(n) ((size_t)(n) << 30)

#define NO_DEFAULT_ASSIGN_COPY_MOVE(Type)                                       \
    Type(const Type&) = delete;                                                 \
    Type& operator=(const Type&) = delete;                                      \
    Type(Type&&) = delete;                                                      \
    Type& operator=(Type&&) = delete;

#define PASTE_LINE_NUMBER(x) PASTE_LINE_NUMBER_(x, __LINE__)
#define PASTE_LINE_NUMBER_(x, line) PASTE_LINE_NUMBER__(x, line)
#define PASTE_LINE_NUMBER__(x, line) x##line

#ifdef __GNUC__

#define ATTRIBUTE_FORMAT_PRINTF(fmt_str_index, arg_index) __attribute__((format(printf, fmt_str_index, arg_index)))
#define DEBUG_BREAK asm("int $3")
#define FORCE_INLINE inline __attribute__((always_inline))

#else

// TODO: not implemented
#define ATTRIBUTE_FORMAT_PRINTF
#define DEBUG_BREAK
#define FORCE_INLINE inline

#endif

void    DebugPrint(const char* format, ...) ATTRIBUTE_FORMAT_PRINTF(1, 2);
void    DebugPrintVA(const char* format, va_list va);

long    GetFileSize(const char* filename);
long    GetFileContents(const char* filename, void* buf, long bufsize);

#define INVALID_CODE_PATH                                                                   \
    do {                                                                                    \
        fprintf(stderr, "INVALID_CODE_PATH: file '%s' line '%d'\n", __FILE__, __LINE__);    \
        DEBUG_BREAK;                                                                        \
    } while (0)


void*   ScratchAlloc(size_t size);
void    ScratchFreeTo(void* ptr);
void    ScratchClear();

#endif
