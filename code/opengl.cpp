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

#include "opengl.h"

#define GLFUNC(type, name) type name;
GLFUNCLIST_CORE
GLFUNCLIST_ARB_DEBUG_OUTPUT
#undef GLFUNC

GLInfo g_GLInfo;

void GL_Init(GL_GetProcAddress_t* glGetProcAddress)
{
#define GLFUNC(type, name) name = (type) glGetProcAddress(#name);

    GLFUNCLIST_CORE

    GLint num_extensions = 0;
    glGetIntegerv(GL_NUM_EXTENSIONS, &num_extensions);

    for (int i = 0; i < num_extensions; ++i)
    {
        const char* extension_name = (const char*) glGetStringi(GL_EXTENSIONS, i);

        if (!g_GLInfo.ARB_debug_output &&
            strcmp(extension_name, "GL_ARB_debug_output") == 0)
        {
            g_GLInfo.ARB_debug_output = true;
            GLFUNCLIST_ARB_DEBUG_OUTPUT
        }
    }

    glGetIntegerv(GL_MAJOR_VERSION, &g_GLInfo.major_version);
    glGetIntegerv(GL_MINOR_VERSION, &g_GLInfo.minor_version);

#undef GLFUNC
}

void GL_InitShaderProgram(GLProgram* program)
{
    program->id = glCreateProgram();
    program->valid = true;

    // compile shaders
    for (int shader_idx = 0; shader_idx < MAX_NUM_SHADERS_PER_PROGRAM; ++shader_idx)
    {
        GLShader* shader = program->shaders + shader_idx;

        if (!shader->name)
            continue;

        shader->id = glCreateShader(shader->stage);
        shader->valid = true;

        long shader_file_size = GetFileSize(shader->name);
        if (shader_file_size < 0)
        {
            fprintf(stderr, "GL_InitShaderProgram: can't read shader '%s'\n", (const char*) shader->name);
            glDeleteShader(shader->id);
            shader->id = 0;
            shader->valid = false;
            continue;
        }

        char* shader_source = (char*) ScratchAlloc(shader_file_size);
        if (GetFileContents(shader->name, shader_source, shader_file_size) != shader_file_size)
        {
            fprintf(stderr, "GL_InitShaderProgram: can't read shader '%s'\n", (const char*) shader->name);
            glDeleteShader(shader->id);
            shader->id = 0;
            shader->valid = false;
            ScratchFreeTo(shader_source);
            continue;
        }

        const GLchar* source_strings[] = {(const GLchar*) shader_source};
        GLint source_string_lengths[] = {(GLint) shader_file_size};
        glShaderSource(shader->id, 1, source_strings, source_string_lengths);

        glCompileShader(shader->id);
        GLint shader_compile_status = 0;
        glGetShaderiv(shader->id, GL_COMPILE_STATUS, &shader_compile_status);
        if (shader_compile_status != GL_TRUE)
        {
            GLint info_log_length = 0;
            glGetShaderiv(shader->id, GL_INFO_LOG_LENGTH, &info_log_length);
            char* info_log = (char*) ScratchAlloc(info_log_length + 1);
            glGetShaderInfoLog(shader->id, info_log_length, NULL, info_log);
            info_log[info_log_length] = '\0';
            fprintf(stderr, "GL_InitShaderProgram: can't compile shader:\n%s\n", info_log);
            ScratchFreeTo(info_log);
            shader->valid = false;
        }

        glAttachShader(program->id, shader->id);
        ScratchFreeTo(shader_source);
    }

    // link program
    glLinkProgram(program->id);
    GLint program_link_status = 0;
    glGetProgramiv(program->id, GL_LINK_STATUS, &program_link_status);
    if (program_link_status != GL_TRUE)
    {
        GLint info_log_length = 0;
        glGetProgramiv(program->id, GL_INFO_LOG_LENGTH, &info_log_length);
        char* info_log = (char*) ScratchAlloc(info_log_length + 1);
        glGetProgramInfoLog(program->id, info_log_length, NULL, info_log);
        info_log[info_log_length] = '\0';
        fprintf(stderr, "GL_InitShaderProgram: can't link program:\n%s\n", info_log);
        ScratchFreeTo(info_log);
        program->valid = false;
    }
}

void GL_FreeShaderProgram(GLProgram* program)
{
    if (program->id)
        glDeleteProgram(program->id);

    program->id = 0;
    program->valid = false;

    for (int shader_idx = 0; shader_idx < MAX_NUM_SHADERS_PER_PROGRAM; ++shader_idx)
    {
        GLShader* shader = program->shaders + shader_idx;

        if (shader->id)
            glDeleteShader(shader->id);

        shader->id = 0;
        shader->valid = false;
    }
}
