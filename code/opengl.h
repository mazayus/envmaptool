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

#ifndef OPENGL_H
#define OPENGL_H

#include "common.h"

#include <GL/glcorearb.h>

#define GLFUNCLIST_CORE                                                                 \
    GLFUNC(PFNGLFLUSHPROC, glFlush)                                                     \
    GLFUNC(PFNGLFINISHPROC, glFinish)                                                   \
    GLFUNC(PFNGLCLEARPROC, glClear)                                                     \
    GLFUNC(PFNGLCLEARCOLORPROC, glClearColor)                                           \
    GLFUNC(PFNGLCLEARDEPTHPROC, glClearDepth)                                           \
    GLFUNC(PFNGLENABLEPROC, glEnable)                                                   \
    GLFUNC(PFNGLDISABLEPROC, glDisable)                                                 \
    GLFUNC(PFNGLGETERRORPROC, glGetError)                                               \
    GLFUNC(PFNGLGETSTRINGPROC, glGetString)                                             \
    GLFUNC(PFNGLGETSTRINGIPROC, glGetStringi)                                           \
    GLFUNC(PFNGLGETINTEGERVPROC, glGetIntegerv)                                         \
    GLFUNC(PFNGLGETINTEGER64VPROC, glGetInteger64v)                                     \
    GLFUNC(PFNGLBLENDCOLORPROC, glBlendColor)                                           \
    GLFUNC(PFNGLBLENDEQUATIONPROC, glBlendEquation)                                     \
    GLFUNC(PFNGLBLENDFUNCPROC, glBlendFunc)                                             \
    GLFUNC(PFNGLDEPTHFUNCPROC, glDepthFunc)                                             \
    GLFUNC(PFNGLVIEWPORTPROC, glViewport)                                               \
    GLFUNC(PFNGLSCISSORPROC, glScissor)                                                 \
    GLFUNC(PFNGLFRONTFACEPROC, glFrontFace)                                             \
    GLFUNC(PFNGLCULLFACEPROC, glCullFace)                                               \
    GLFUNC(PFNGLPOLYGONMODEPROC, glPolygonMode)                                         \
    GLFUNC(PFNGLPOLYGONOFFSETPROC, glPolygonOffset)                                     \
    GLFUNC(PFNGLPIXELSTOREFPROC, glPixelStoref)                                         \
    GLFUNC(PFNGLPIXELSTOREIPROC, glPixelStorei)                                         \
    GLFUNC(PFNGLGENBUFFERSPROC, glGenBuffers)                                           \
    GLFUNC(PFNGLDELETEBUFFERSPROC, glDeleteBuffers)                                     \
    GLFUNC(PFNGLBINDBUFFERPROC, glBindBuffer)                                           \
    GLFUNC(PFNGLBINDBUFFERBASEPROC, glBindBufferBase)                                   \
    GLFUNC(PFNGLBINDBUFFERRANGEPROC, glBindBufferRange)                                 \
    GLFUNC(PFNGLBUFFERDATAPROC, glBufferData)                                           \
    GLFUNC(PFNGLBUFFERSUBDATAPROC, glBufferSubData)                                     \
    GLFUNC(PFNGLMAPBUFFERPROC, glMapBuffer)                                             \
    GLFUNC(PFNGLMAPBUFFERRANGEPROC, glMapBufferRange)                                   \
    GLFUNC(PFNGLUNMAPBUFFERPROC, glUnmapBuffer)                                         \
    GLFUNC(PFNGLACTIVETEXTUREPROC, glActiveTexture)                                     \
    GLFUNC(PFNGLGENTEXTURESPROC, glGenTextures)                                         \
    GLFUNC(PFNGLDELETETEXTURESPROC, glDeleteTextures)                                   \
    GLFUNC(PFNGLBINDTEXTUREPROC, glBindTexture)                                         \
    GLFUNC(PFNGLTEXIMAGE2DPROC, glTexImage2D)                                           \
    GLFUNC(PFNGLTEXSUBIMAGE2DPROC, glTexSubImage2D)                                     \
    GLFUNC(PFNGLGENERATEMIPMAPPROC, glGenerateMipmap)                                   \
    GLFUNC(PFNGLTEXPARAMETERFPROC, glTexParameterf)                                     \
    GLFUNC(PFNGLTEXPARAMETERFVPROC, glTexParameterfv)                                   \
    GLFUNC(PFNGLTEXPARAMETERIPROC, glTexParameteri)                                     \
    GLFUNC(PFNGLTEXPARAMETERIVPROC, glTexParameteriv)                                   \
    GLFUNC(PFNGLGETTEXIMAGEPROC, glGetTexImage)                                         \
    GLFUNC(PFNGLGETTEXPARAMETERIVPROC, glGetTexParameteriv)                             \
    GLFUNC(PFNGLGETTEXPARAMETERFVPROC, glGetTexParameterfv)                             \
    GLFUNC(PFNGLGETTEXLEVELPARAMETERIVPROC, glGetTexLevelParameteriv)                   \
    GLFUNC(PFNGLGETTEXLEVELPARAMETERFVPROC, glGetTexLevelParameterfv)                   \
    GLFUNC(PFNGLCREATESHADERPROC, glCreateShader)                                       \
    GLFUNC(PFNGLDELETESHADERPROC, glDeleteShader)                                       \
    GLFUNC(PFNGLSHADERSOURCEPROC, glShaderSource)                                       \
    GLFUNC(PFNGLCOMPILESHADERPROC, glCompileShader)                                     \
    GLFUNC(PFNGLGETSHADERINFOLOGPROC, glGetShaderInfoLog)                               \
    GLFUNC(PFNGLGETSHADERIVPROC, glGetShaderiv)                                         \
    GLFUNC(PFNGLCREATEPROGRAMPROC, glCreateProgram)                                     \
    GLFUNC(PFNGLDELETEPROGRAMPROC, glDeleteProgram)                                     \
    GLFUNC(PFNGLATTACHSHADERPROC, glAttachShader)                                       \
    GLFUNC(PFNGLDETACHSHADERPROC, glDetachShader)                                       \
    GLFUNC(PFNGLLINKPROGRAMPROC, glLinkProgram)                                         \
    GLFUNC(PFNGLVALIDATEPROGRAMPROC, glValidateProgram)                                 \
    GLFUNC(PFNGLUSEPROGRAMPROC, glUseProgram)                                           \
    GLFUNC(PFNGLGETPROGRAMINFOLOGPROC, glGetProgramInfoLog)                             \
    GLFUNC(PFNGLGETPROGRAMIVPROC, glGetProgramiv)                                       \
    GLFUNC(PFNGLGENVERTEXARRAYSPROC, glGenVertexArrays)                                 \
    GLFUNC(PFNGLDELETEVERTEXARRAYSPROC, glDeleteVertexArrays)                           \
    GLFUNC(PFNGLBINDVERTEXARRAYPROC, glBindVertexArray)                                 \
    GLFUNC(PFNGLENABLEVERTEXATTRIBARRAYPROC, glEnableVertexAttribArray)                 \
    GLFUNC(PFNGLDISABLEVERTEXATTRIBARRAYPROC, glDisableVertexAttribArray)               \
    GLFUNC(PFNGLVERTEXATTRIBPOINTERPROC, glVertexAttribPointer)                         \
    GLFUNC(PFNGLVERTEXATTRIBIPOINTERPROC, glVertexAttribIPointer)                       \
    GLFUNC(PFNGLDRAWARRAYSPROC, glDrawArrays)                                           \
    GLFUNC(PFNGLDRAWELEMENTSPROC, glDrawElements)                                       \
    GLFUNC(PFNGLGENQUERIESPROC, glGenQueries)                                           \
    GLFUNC(PFNGLDELETEQUERIESPROC, glDeleteQueries)                                     \
    GLFUNC(PFNGLBEGINQUERYPROC, glBeginQuery)                                           \
    GLFUNC(PFNGLENDQUERYPROC, glEndQuery)                                               \
    GLFUNC(PFNGLQUERYCOUNTERPROC, glQueryCounter)                                       \
    GLFUNC(PFNGLGETQUERYIVPROC, glGetQueryiv)                                           \
    GLFUNC(PFNGLGETQUERYOBJECTIVPROC, glGetQueryObjectiv)                               \
    GLFUNC(PFNGLGETQUERYOBJECTI64VPROC, glGetQueryObjecti64v)                           \
    GLFUNC(PFNGLGETQUERYOBJECTUIVPROC, glGetQueryObjectuiv)                             \
    GLFUNC(PFNGLGETQUERYOBJECTUI64VPROC, glGetQueryObjectui64v)                         \
    GLFUNC(PFNGLGENFRAMEBUFFERSPROC, glGenFramebuffers)                                 \
    GLFUNC(PFNGLDELETEFRAMEBUFFERSPROC, glDeleteFramebuffers)                           \
    GLFUNC(PFNGLBINDFRAMEBUFFERPROC, glBindFramebuffer)                                 \
    GLFUNC(PFNGLCHECKFRAMEBUFFERSTATUSPROC, glCheckFramebufferStatus)                   \
    GLFUNC(PFNGLBLITFRAMEBUFFERPROC, glBlitFramebuffer)                                 \
    GLFUNC(PFNGLFRAMEBUFFERTEXTUREPROC, glFramebufferTexture)                           \
    GLFUNC(PFNGLFRAMEBUFFERTEXTURE2DPROC, glFramebufferTexture2D)                       \
    GLFUNC(PFNGLGETUNIFORMLOCATIONPROC, glGetUniformLocation)                           \
    GLFUNC(PFNGLUNIFORM1FPROC, glUniform1f)                                             \
    GLFUNC(PFNGLUNIFORM1FVPROC, glUniform1fv)                                           \
    GLFUNC(PFNGLUNIFORM1IPROC, glUniform1i)                                             \
    GLFUNC(PFNGLUNIFORM1IVPROC, glUniform1iv)                                           \
    GLFUNC(PFNGLUNIFORM2FPROC, glUniform2f)                                             \
    GLFUNC(PFNGLUNIFORM2FVPROC, glUniform2fv)                                           \
    GLFUNC(PFNGLUNIFORM2IPROC, glUniform2i)                                             \
    GLFUNC(PFNGLUNIFORM2IVPROC, glUniform2iv)                                           \
    GLFUNC(PFNGLUNIFORM3FPROC, glUniform3f)                                             \
    GLFUNC(PFNGLUNIFORM3FVPROC, glUniform3fv)                                           \
    GLFUNC(PFNGLUNIFORM3IPROC, glUniform3i)                                             \
    GLFUNC(PFNGLUNIFORM3IVPROC, glUniform3iv)                                           \
    GLFUNC(PFNGLUNIFORM4FPROC, glUniform4f)                                             \
    GLFUNC(PFNGLUNIFORM4FVPROC, glUniform4fv)                                           \
    GLFUNC(PFNGLUNIFORM4IPROC, glUniform4i)                                             \
    GLFUNC(PFNGLUNIFORM4IVPROC, glUniform4iv)                                           \
    GLFUNC(PFNGLUNIFORMMATRIX3FVPROC, glUniformMatrix3fv)                               \
    GLFUNC(PFNGLUNIFORMMATRIX4FVPROC, glUniformMatrix4fv)

#define GLFUNCLIST_ARB_DEBUG_OUTPUT                                                     \
    GLFUNC(PFNGLDEBUGMESSAGECONTROLARBPROC, glDebugMessageControlARB)                   \
    GLFUNC(PFNGLDEBUGMESSAGEINSERTARBPROC, glDebugMessageInsertARB)                     \
    GLFUNC(PFNGLDEBUGMESSAGECALLBACKARBPROC, glDebugMessageCallbackARB)                 \
    GLFUNC(PFNGLGETDEBUGMESSAGELOGARBPROC, glGetDebugMessageLogARB)

#define GLFUNC(type, name) extern type name;
GLFUNCLIST_CORE
GLFUNCLIST_ARB_DEBUG_OUTPUT
#undef GLFUNC

//

struct GLInfo
{
    bool ARB_debug_output;

    int major_version;
    int minor_version;
};

typedef void* GL_GetProcAddress_t(const char*);
void GL_Init(GL_GetProcAddress_t* glGetProcAddress);

extern GLInfo g_GLInfo;

//

#define MAX_NUM_SHADERS_PER_PROGRAM 8

struct GLShader
{
    const char*     name;
    GLenum          stage;

    GLuint          id;
    bool            valid;
};

struct GLProgram
{
    const char*     name;
    GLShader        shaders[MAX_NUM_SHADERS_PER_PROGRAM];

    GLuint          id;
    bool            valid;
};

void GL_InitShaderProgram(GLProgram* program);
void GL_FreeShaderProgram(GLProgram* program);

#endif
