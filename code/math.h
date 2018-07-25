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

#ifndef MATH_H
#define MATH_H

#include "common.h"

extern "C" double exp(double);
extern "C" double sin(double);
extern "C" double cos(double);
extern "C" double tan(double);
extern "C" double sqrt(double);
extern "C" double floor(double);
extern "C" double ceil(double);
extern "C" double asin(double);
extern "C" double acos(double);
extern "C" double atan2(double, double);

namespace Math
{
    static const float PI = 3.141592653589793f;
    static const float SQRT_2 = 1.414213562373095f;

    float Abs(float x);
    float Min(float a, float b);
    float Max(float a, float b);
    float Lerp(float x0, float x1, float t);
    float Clamp(float x, float min, float max);
}

/*
 * Vector3
 */

struct Vector3
{
    union
    {
        float data[3];

        struct { float x, y, z; };
        struct { float r, g, b; };
    };

    Vector3() {}
    Vector3(float x, float y, float z) : data {x, y, z} {}
};

Vector3 operator+(Vector3 v);
Vector3 operator-(Vector3 v);
Vector3 operator+(Vector3 v1, Vector3 v2);
Vector3 operator-(Vector3 v1, Vector3 v2);
Vector3& operator+=(Vector3& v1, Vector3 v2);
Vector3& operator-=(Vector3& v1, Vector3 v2);

Vector3 operator*(float s, Vector3 v);
Vector3 operator*(Vector3 v, float s);

namespace Math
{
    float Length(Vector3 v);
    float LengthSq(Vector3 v);
    Vector3 Normalize(Vector3 v);
    Vector3 Lerp(Vector3 v0, Vector3 v1, float t);
    float Dot(Vector3 v1, Vector3 v2);
    Vector3 Cross(Vector3 v1, Vector3 v2);
}

Vector3 Abs(Vector3 v);
float MaxElement(Vector3 v);
float MinElement(Vector3 v);

/*
 * Vector4
 */

struct Vector4
{
    union
    {
        float data[4];

        struct { float x, y, z, w; };
        Vector3 xyz;

        struct { float r, g, b, a; };
        Vector3 rgb;
    };

    Vector4() {}
    Vector4(float x, float y, float z, float w) : data {x, y, z, w} {}
    Vector4(Vector3 xyz, float w) : data {xyz.x, xyz.y, xyz.z, w} {}
};

Vector4 operator+(Vector4 v);
Vector4 operator-(Vector4 v);
Vector4 operator+(Vector4 v1, Vector4 v2);
Vector4 operator-(Vector4 v1, Vector4 v2);

Vector4 operator*(float s, Vector4 v);
Vector4 operator*(Vector4 v, float s);

namespace Math
{
    float Length(Vector4 v);
    float LengthSq(Vector4 v);
    Vector4 Normalize(Vector4 v);
    Vector4 Lerp(Vector4 v0, Vector4 v1, float t);
    float Dot(Vector4 v1, Vector4 v2);
}

/*
 * Matrix4
 */

struct Matrix4
{
    union
    {
        float data[16];

        struct
        {
            float _11, _21, _31, _41;
            float _12, _22, _32, _42;
            float _13, _23, _33, _43;
            float _14, _24, _34, _44;
        };
    };

    static Matrix4 MakeIdentity();
    static Matrix4 MakeRotation(Vector3 axis, float angle);
    static Matrix4 MakeTranslation(float x, float y, float z);

    static Matrix4 MakeOrtho(float l, float r, float b, float t, float znear, float zfar);

    static Matrix4 MakePerspective(float l, float r, float b, float t, float z_near, float z_far);
    static Matrix4 MakePerspectiveInverse(float l, float r, float b, float t, float z_near, float z_far);

    static Matrix4 MakePerspective(float fovy, float aspect, float znear, float zfar);
    static Matrix4 MakePerspectiveInverse(float fovy, float aspect, float znear, float zfar);
};

Matrix4 operator+(Matrix4 m);
Matrix4 operator-(Matrix4 m);
Matrix4 operator*(Matrix4 m1, Matrix4 m2);
Matrix4 operator+(Matrix4 m1, Matrix4 m2);
Matrix4 operator-(Matrix4 m1, Matrix4 m2);

Vector4 operator*(Matrix4 m, Vector4 v);

/*
 * Quaternion
 */

struct Quaternion
{
    float x, y, z, w;

    static Quaternion MakeIdentity();
};

Quaternion operator+(Quaternion q);
Quaternion operator-(Quaternion q);
Quaternion operator*(Quaternion q1, Quaternion q2);
Quaternion operator+(Quaternion q1, Quaternion q2);
Quaternion operator-(Quaternion q1, Quaternion q2);

Quaternion operator*(float s, Quaternion q);
Quaternion operator*(Quaternion q, float s);

namespace Math
{
    float Norm(Quaternion q);
    float NormSq(Quaternion q);
    Quaternion Conjugate(Quaternion q);
    Quaternion Inverse(Quaternion q);
    Quaternion Lerp(Quaternion q0, Quaternion q1, float t);
}

/*
 * Transform
 */

struct Transform
{
    Vector3     translation;
    Quaternion  rotation;
    float       scale;

    static Transform MakeIdentity();

    Matrix4 GetMatrix() const;
};

Transform operator*(Transform t1, Transform t2);

//
//
//

Vector3 TransformVector(Quaternion q, Vector3 v);

void ExtractClippingPlanes(Matrix4 m, Vector4 planes[6]);

Matrix4 MakeLocalToWorldMatrix(Vector3 position, Vector3 basis[3]);
Matrix4 MakeWorldToLocalMatrix(Vector3 position, Vector3 basis[3]);
Matrix4 MakeLocalToWorldMatrix(Vector3 position, Quaternion orientation, float scale);
Matrix4 MakeWorldToLocalMatrix(Vector3 position, Quaternion orientation, float scale);

Quaternion AxisAngleToQuaternion(Vector3 axis, float angle);

void GetFrustumVertices(float fovy, float aspect, float znear, float zfar, Vector3 vertices[8]);
void GetFrustumVertices(float l, float r, float b, float t, float znear, float zfar, Vector3 vertices[8]);

#endif
