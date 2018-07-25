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

#include "math.h"

float Math::Abs(float x)
{
    return x >= 0 ? x : -x;
}

float Math::Min(float a, float b)
{
    return a < b ? a : b;
}

float Math::Max(float a, float b)
{
    return a > b ? a : b;
}

float Math::Lerp(float x0, float x1, float t)
{
    return (1 - t) * x0 + t * x1;
}

float Math::Clamp(float x, float min, float max)
{
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

/*
 * Vector3
 */

Vector3 operator+(Vector3 v)
{
    return v;
}

Vector3 operator-(Vector3 v)
{
    return Vector3(-v.x, -v.y, -v.z);
}

Vector3 operator+(Vector3 v1, Vector3 v2)
{
    return Vector3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

Vector3 operator-(Vector3 v1, Vector3 v2)
{
    return Vector3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

Vector3& operator+=(Vector3& v1, Vector3 v2)
{
    v1 = v1 + v2;
    return v1;
}

Vector3& operator-=(Vector3& v1, Vector3 v2)
{
    v1 = v1 - v2;
    return v1;
}

Vector3 operator*(float s, Vector3 v)
{
    return Vector3(s * v.x, s * v.y, s * v.z);
}

Vector3 operator*(Vector3 v, float s)
{
    return Vector3(v.x * s, v.y * s, v.z * s);
}

float Math::Length(Vector3 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

float Math::LengthSq(Vector3 v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

Vector3 Math::Normalize(Vector3 v)
{
    const float epsilon = 1e-7;

    float len = Math::Length(v);
    if (len < epsilon)
        return Vector3(0, 0, 0);

    return 1.0f / len * v;
}

Vector3 Math::Lerp(Vector3 v0, Vector3 v1, float t)
{
    return (1 - t) * v0 + t * v1;
}

float Math::Dot(Vector3 v1, Vector3 v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vector3 Math::Cross(Vector3 v1, Vector3 v2)
{
    Vector3 res;

    res.x = v1.y * v2.z - v2.y * v1.z;
    res.y = - (v1.x * v2.z - v2.x * v1.z);
    res.z = v1.x * v2.y - v2.x * v1.y;

    return res;
}

Vector3 Abs(Vector3 v)
{
    Vector3 res;

    res.x = Math::Abs(v.x);
    res.y = Math::Abs(v.y);
    res.z = Math::Abs(v.z);

    return res;
}

float MaxElement(Vector3 v)
{
    return Math::Max(v.x, Math::Max(v.y, v.z));
}

float MinElement(Vector3 v)
{
    return Math::Min(v.x, Math::Min(v.y, v.z));
}

/*
 * Vector4
 */

Vector4 operator+(Vector4 v)
{
    return v;
}

Vector4 operator-(Vector4 v)
{
    return Vector4(-v.x, -v.y, -v.z, -v.w);
}

Vector4 operator+(Vector4 v1, Vector4 v2)
{
    return Vector4(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w + v2.w);
}

Vector4 operator-(Vector4 v1, Vector4 v2)
{
    return Vector4(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w);
}

Vector4 operator*(float s, Vector4 v)
{
    return Vector4(s * v.x, s * v.y, s * v.z, s * v.w);
}

Vector4 operator*(Vector4 v, float s)
{
    return Vector4(v.x * s, v.y * s, v.z * s, v.w * s);
}

float Math::Length(Vector4 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
}

float Math::LengthSq(Vector4 v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w;
}

Vector4 Math::Normalize(Vector4 v)
{
    const float epsilon = 1e-7;

    float len = Math::Length(v);
    if (len < epsilon)
        return Vector4(0, 0, 0, 0);

    return 1.0f / len * v;
}

Vector4 Math::Lerp(Vector4 v0, Vector4 v1, float t)
{
    return (1 - t) * v0 + t * v1;
}

float Math::Dot(Vector4 v1, Vector4 v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w;
}

/*
 * Matrix4
 */

Matrix4 Matrix4::MakeIdentity()
{
    Matrix4 res = {};

    res._11 = res._22 = res._33 = res._44 = 1.0f;

    return res;
}

Matrix4 Matrix4::MakeRotation(Vector3 axis, float angle)
{
    Matrix4 res = Matrix4::MakeIdentity();

    Vector3 r = Math::Normalize(axis);
    float s = sin(angle), c = cos(angle);

    res._11 = (1 - c) * r.x * r.x + c;
    res._12 = (1 - c) * r.x * r.y - s * r.z;
    res._13 = (1 - c) * r.x * r.z + s * r.y;

    res._21 = (1 - c) * r.y * r.x + s * r.z;
    res._22 = (1 - c) * r.y * r.y + c;
    res._23 = (1 - c) * r.y * r.z - s * r.x;

    res._31 = (1 - c) * r.z * r.x - s * r.y;
    res._32 = (1 - c) * r.z * r.y + s * r.x;
    res._33 = (1 - c) * r.z * r.z + c;

    return res;
}

Matrix4 Matrix4::MakeTranslation(float x, float y, float z)
{
    Matrix4 res = Matrix4::MakeIdentity();

    res._14 = x; res._24 = y; res._34 = z;

    return res;
}

Matrix4 Matrix4::MakeOrtho(float left, float right, float bottom, float top, float znear, float zfar)
{
    Matrix4 res = {};

    res._11 = 2.0f / (right - left);
    res._14 = - (right + left) / (right - left);

    res._22 = 2.0f / (top - bottom);
    res._24 = - (top + bottom) / (top - bottom);

    res._33 = - 2.0f / (zfar - znear);
    res._34 = - (zfar + znear) / (zfar - znear);

    res._44 = 1.0f;

    return res;
}

Matrix4 Matrix4::MakePerspective(float l, float r, float b, float t, float z_near, float z_far)
{
    Matrix4 res = {};

    res._11 = 2 * z_near / (r - l);
    res._13 = (r + l) / (r - l);

    res._22 = 2 * z_near / (t - b);
    res._23 = (t + b) / (t - b);

    res._33 = (z_near + z_far) / (z_near - z_far);
    res._34 = 2 * z_near * z_far / (z_near - z_far);

    res._43 = -1.0f;
    res._44 = 0.0f;

    return res;
}

Matrix4 Matrix4::MakePerspectiveInverse(float l, float r, float b, float t, float z_near, float z_far)
{
    Matrix4 res = {};

    res._11 = (r - l) / (2 * z_near);
    res._14 = (r + l) / (2 * z_near);

    res._22 = (t - b) / (2 * z_near);
    res._23 = (t + b) / (2 * z_near);

    res._33 = 0.0f;
    res._34 = -1.0f;

    res._43 = (z_near - z_far) / (2 * z_near * z_far);
    res._44 = (z_near + z_far) / (2 * z_near * z_far);

    return res;
}

Matrix4 Matrix4::MakePerspective(float fovy, float aspect, float znear, float zfar)
{
    float half_height = tan(fovy / 2) * znear;
    float half_width = aspect * half_height;

    return MakePerspective(-half_width, half_width, -half_height, half_height, znear, zfar);
}

Matrix4 Matrix4::MakePerspectiveInverse(float fovy, float aspect, float znear, float zfar)
{
    float half_height = tan(fovy / 2) * znear;
    float half_width = aspect * half_height;

    return MakePerspectiveInverse(-half_width, half_width, -half_height, half_height, znear, zfar);
}

Matrix4 operator+(Matrix4 m)
{
    return m;
}

Matrix4 operator-(Matrix4 m)
{
    Matrix4 res;

    res._11 = -m._11; res._12 = -m._12; res._13 = -m._13; res._14 = -m._14;
    res._21 = -m._21; res._22 = -m._22; res._23 = -m._23; res._24 = -m._24;
    res._31 = -m._31; res._32 = -m._32; res._33 = -m._33; res._34 = -m._34;
    res._41 = -m._41; res._42 = -m._42; res._43 = -m._43; res._44 = -m._44;

    return res;
}

Matrix4 operator*(Matrix4 m1, Matrix4 m2)
{
    Matrix4 res;

    res._11 = m1._11 * m2._11 + m1._12 * m2._21 + m1._13 * m2._31 + m1._14 * m2._41;
    res._12 = m1._11 * m2._12 + m1._12 * m2._22 + m1._13 * m2._32 + m1._14 * m2._42;
    res._13 = m1._11 * m2._13 + m1._12 * m2._23 + m1._13 * m2._33 + m1._14 * m2._43;
    res._14 = m1._11 * m2._14 + m1._12 * m2._24 + m1._13 * m2._34 + m1._14 * m2._44;

    res._21 = m1._21 * m2._11 + m1._22 * m2._21 + m1._23 * m2._31 + m1._24 * m2._41;
    res._22 = m1._21 * m2._12 + m1._22 * m2._22 + m1._23 * m2._32 + m1._24 * m2._42;
    res._23 = m1._21 * m2._13 + m1._22 * m2._23 + m1._23 * m2._33 + m1._24 * m2._43;
    res._24 = m1._21 * m2._14 + m1._22 * m2._24 + m1._23 * m2._34 + m1._24 * m2._44;

    res._31 = m1._31 * m2._11 + m1._32 * m2._21 + m1._33 * m2._31 + m1._34 * m2._41;
    res._32 = m1._31 * m2._12 + m1._32 * m2._22 + m1._33 * m2._32 + m1._34 * m2._42;
    res._33 = m1._31 * m2._13 + m1._32 * m2._23 + m1._33 * m2._33 + m1._34 * m2._43;
    res._34 = m1._31 * m2._14 + m1._32 * m2._24 + m1._33 * m2._34 + m1._34 * m2._44;

    res._41 = m1._41 * m2._11 + m1._42 * m2._21 + m1._43 * m2._31 + m1._44 * m2._41;
    res._42 = m1._41 * m2._12 + m1._42 * m2._22 + m1._43 * m2._32 + m1._44 * m2._42;
    res._43 = m1._41 * m2._13 + m1._42 * m2._23 + m1._43 * m2._33 + m1._44 * m2._43;
    res._44 = m1._41 * m2._14 + m1._42 * m2._24 + m1._43 * m2._34 + m1._44 * m2._44;

    return res;
}

Matrix4 operator+(Matrix4 m1, Matrix4 m2)
{
    Matrix4 res;

    res._11 = m1._11 + m2._11; res._12 = m1._12 + m2._12; res._13 = m1._13 + m2._13; res._14 = m1._14 + m2._14;
    res._21 = m1._21 + m2._21; res._22 = m1._22 + m2._22; res._23 = m1._23 + m2._23; res._24 = m1._24 + m2._24;
    res._31 = m1._31 + m2._31; res._32 = m1._32 + m2._32; res._33 = m1._33 + m2._33; res._34 = m1._34 + m2._34;
    res._41 = m1._41 + m2._41; res._42 = m1._42 + m2._42; res._43 = m1._43 + m2._43; res._44 = m1._44 + m2._44;

    return res;
}

Matrix4 operator-(Matrix4 m1, Matrix4 m2)
{
    Matrix4 res;

    res._11 = m1._11 - m2._11; res._12 = m1._12 - m2._12; res._13 = m1._13 - m2._13; res._14 = m1._14 - m2._14;
    res._21 = m1._21 - m2._21; res._22 = m1._22 - m2._22; res._23 = m1._23 - m2._23; res._24 = m1._24 - m2._24;
    res._31 = m1._31 - m2._31; res._32 = m1._32 - m2._32; res._33 = m1._33 - m2._33; res._34 = m1._34 - m2._34;
    res._41 = m1._41 - m2._41; res._42 = m1._42 - m2._42; res._43 = m1._43 - m2._43; res._44 = m1._44 - m2._44;

    return res;
}

Vector4 operator*(Matrix4 m, Vector4 v)
{
    Vector4 res;

    res.x = m._11 * v.x + m._12 * v.y + m._13 * v.z + m._14 * v.w;
    res.y = m._21 * v.x + m._22 * v.y + m._23 * v.z + m._24 * v.w;
    res.z = m._31 * v.x + m._32 * v.y + m._33 * v.z + m._34 * v.w;
    res.w = m._41 * v.x + m._42 * v.y + m._43 * v.z + m._44 * v.w;

    return res;
}

/*
 * Quaternion
 */

Quaternion Quaternion::MakeIdentity()
{
    return {0, 0, 0, 1};
}

Quaternion operator+(Quaternion q)
{
    return q;
}

Quaternion operator-(Quaternion q)
{
    return {-q.x, -q.y, -q.z, -q.w};
}

Quaternion operator*(Quaternion q1, Quaternion q2)
{
    Quaternion res = {};

    res.x = q1.y * q2.z - q1.z * q2.y + q2.w * q1.x + q1.w * q2.x;
    res.y = q1.z * q2.x - q1.x * q2.z + q2.w * q1.y + q1.w * q2.y;
    res.z = q1.x * q2.y - q1.y * q2.x + q2.w * q1.z + q1.w * q2.z;
    res.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;

    return res;
}

Quaternion operator+(Quaternion q1, Quaternion q2)
{
    return {q1.x + q2.x, q1.y + q2.y, q1.z + q2.z, q1.w + q2.w};
}

Quaternion operator-(Quaternion q1, Quaternion q2)
{
    return {q1.x - q2.x, q1.y - q2.y, q1.z - q2.z, q1.w - q2.w};
}

Quaternion operator*(float s, Quaternion q)
{
    return {s*q.x, s*q.y, s*q.z, s*q.w};
}

Quaternion operator*(Quaternion q, float s)
{
    return {q.x*s, q.y*s, q.z*s, q.w*s};
}

float Math::Norm(Quaternion q)
{
    return sqrt(q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w);
}

float Math::NormSq(Quaternion q)
{
    return q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w;
}

Quaternion Math::Conjugate(Quaternion q)
{
    return {-q.x, -q.y, -q.z, q.w};
}

Quaternion Math::Inverse(Quaternion q)
{
    return Math::Conjugate(q) * (1.0f / Math::NormSq(q));
}

Quaternion Math::Lerp(Quaternion q0, Quaternion q1, float t)
{
    Vector4 q0v = Vector4(q0.x, q0.y, q0.z, q0.w);
    Vector4 q1v = Vector4(q1.x, q1.y, q1.z, q1.w);
    float s = (Math::Dot(q0v, q1v) > 0) ? 1 : -1;
    return (1 - t) * q0 + s * t * q1;
}

/*
 * Transform
 */

Transform Transform::MakeIdentity()
{
    Transform res;

    res.translation = Vector3(0, 0, 0);
    res.rotation = Quaternion::MakeIdentity();
    res.scale = 1;

    return res;
}

Matrix4 Transform::GetMatrix() const
{
    // TODO: remove this function, its functionality was moved to MakeLocalToWorldMatrix
    return MakeLocalToWorldMatrix(translation, rotation, scale);
}

Transform operator*(Transform t1, Transform t2)
{
    Transform res;

    Quaternion tmp = {t2.translation.x, t2.translation.y, t2.translation.z, 0.0f};
    tmp = t1.rotation * tmp * Math::Conjugate(t1.rotation);
    res.translation = t1.translation + t1.scale * Vector3(tmp.x, tmp.y, tmp.z);

    res.rotation = t1.rotation * t2.rotation;
    res.scale = t1.scale * t2.scale;

    return res;
}

Vector3 TransformVector(Quaternion q, Vector3 v)
{
    Quaternion vq = {v.x, v.y, v.z, 0.0};
    vq = q * vq * Math::Inverse(q);
    return Vector3(vq.x, vq.y, vq.z);
}

void ExtractClippingPlanes(Matrix4 m, Vector4 planes[6])
{
    planes[0] = Vector4(m._41 + m._31, m._42 + m._32, m._43 + m._33, m._44 + m._34); // near
    planes[1] = Vector4(m._41 - m._31, m._42 - m._32, m._43 - m._33, m._44 - m._34); // far
    planes[2] = Vector4(m._41 + m._21, m._42 + m._22, m._43 + m._23, m._44 + m._24); // bottom
    planes[3] = Vector4(m._41 - m._21, m._42 - m._22, m._43 - m._23, m._44 - m._24); // top
    planes[4] = Vector4(m._41 + m._11, m._42 + m._12, m._43 + m._13, m._44 + m._14); // left
    planes[5] = Vector4(m._41 - m._11, m._42 - m._12, m._43 - m._13, m._44 - m._14); // right

    for (int i = 0; i < 6; ++i)
    {
        float normal_length = Math::Length(planes[i].xyz);
        planes[i] = (1.0f / normal_length) * planes[i];
    }
}

/*
 * MakeLocalToWorldMatrix - returns a matrix representing transformation from local frame to world frame
 *     position - position of the local frame expressed in world coordinates
 *     basis - basis vectors of the local frame expressed in world coordinates
 */
// TODO: add scale parameter
Matrix4 MakeLocalToWorldMatrix(Vector3 position, Vector3 basis[3])
{
    Matrix4 matrix = Matrix4::MakeIdentity();

    matrix._11 = basis[0].x;
    matrix._21 = basis[0].y;
    matrix._31 = basis[0].z;

    matrix._12 = basis[1].x;
    matrix._22 = basis[1].y;
    matrix._32 = basis[1].z;

    matrix._13 = basis[2].x;
    matrix._23 = basis[2].y;
    matrix._33 = basis[2].z;

    matrix._14 = position.x;
    matrix._24 = position.y;
    matrix._34 = position.z;

    return matrix;
}

/*
 * MakeWorldToLocalMatrix - returns a matrix representing transformation from world frame to local frame
 *     position - position of the local frame expressed in world coordinates
 *     basis - basis vectors of the local frame expressed in world coordinates
 */
// TODO: add scale parameter
Matrix4 MakeWorldToLocalMatrix(Vector3 position, Vector3 basis[3])
{
    Matrix4 matrix = Matrix4::MakeIdentity();

    matrix._11 = basis[0].x;
    matrix._21 = basis[1].x;
    matrix._31 = basis[2].x;

    matrix._12 = basis[0].y;
    matrix._22 = basis[1].y;
    matrix._32 = basis[2].y;

    matrix._13 = basis[0].z;
    matrix._23 = basis[1].z;
    matrix._33 = basis[2].z;

    matrix._14 = -Math::Dot(basis[0], position);
    matrix._24 = -Math::Dot(basis[1], position);
    matrix._34 = -Math::Dot(basis[2], position);

    return matrix;
}

// TODO: express this in terms of MakeLocalToWorldMatrix(position, basis)
Matrix4 MakeLocalToWorldMatrix(Vector3 position, Quaternion orientation, float scale)
{
    Matrix4 matrix = Matrix4::MakeIdentity();

    float qxqx = 2 * orientation.x * orientation.x;
    float qyqy = 2 * orientation.y * orientation.y;
    float qzqz = 2 * orientation.z * orientation.z;
    float qwqx = 2 * orientation.w * orientation.x;
    float qwqy = 2 * orientation.w * orientation.y;
    float qwqz = 2 * orientation.w * orientation.z;
    float qxqy = 2 * orientation.x * orientation.y;
    float qxqz = 2 * orientation.x * orientation.z;
    float qyqz = 2 * orientation.y * orientation.z;

    matrix._11 = scale * (1 - (qyqy + qzqz));
    matrix._12 = scale * (qxqy - qwqz);
    matrix._13 = scale * (qxqz + qwqy);
    matrix._14 = position.x;

    matrix._21 = scale * (qxqy + qwqz);
    matrix._22 = scale * (1 - (qxqx + qzqz));
    matrix._23 = scale * (qyqz - qwqx);
    matrix._24 = position.y;

    matrix._31 = scale * (qxqz - qwqy);
    matrix._32 = scale * (qyqz + qwqx);
    matrix._33 = scale * (1 - (qxqx + qyqy));
    matrix._34 = position.z;

    return matrix;
}

// TODO: express this in terms of MakeWorldToLocalMatrix(position, basis)
Matrix4 MakeWorldToLocalMatrix(Vector3 position, Quaternion orientation, float scale)
{
    Matrix4 matrix = Matrix4::MakeIdentity();

    float qxqx = 2 * orientation.x * orientation.x;
    float qyqy = 2 * orientation.y * orientation.y;
    float qzqz = 2 * orientation.z * orientation.z;
    float qwqx = 2 * orientation.w * orientation.x;
    float qwqy = 2 * orientation.w * orientation.y;
    float qwqz = 2 * orientation.w * orientation.z;
    float qxqy = 2 * orientation.x * orientation.y;
    float qxqz = 2 * orientation.x * orientation.z;
    float qyqz = 2 * orientation.y * orientation.z;

    Vector3 p = -position;
    float s = 1.0f / scale;

    matrix._11 = s * (1 - (qyqy + qzqz));
    matrix._12 = s * (qxqy + qwqz);
    matrix._13 = s * (qxqz - qwqy);
    matrix._14 = matrix._11 * p.x + matrix._12 * p.y + matrix._13 * p.z;

    matrix._21 = s * (qxqy - qwqz);
    matrix._22 = s * (1 - (qxqx + qzqz));
    matrix._23 = s * (qyqz + qwqx);
    matrix._24 = matrix._21 * p.x + matrix._22 * p.y + matrix._23 * p.z;

    matrix._31 = s * (qxqz + qwqy);
    matrix._32 = s * (qyqz - qwqx);
    matrix._33 = s * (1 - (qxqx + qyqy));
    matrix._34 = matrix._31 * p.x + matrix._32 * p.y + matrix._33 * p.z;

    return matrix;
}

Quaternion AxisAngleToQuaternion(Vector3 axis, float angle)
{
    Quaternion res;

    axis = Math::Normalize(axis) * sin(angle / 2);
    res.x = axis.x; res.y = axis.y; res.z = axis.z; res.w = cos(angle / 2);

    return res;
}

void GetFrustumVertices(float fovy, float aspect, float znear, float zfar, Vector3 vertices[8])
{
    float tan_half_fov = tan(fovy / 2);

    float half_height = tan_half_fov * znear;
    float half_width = half_height * aspect;

    GetFrustumVertices(-half_width, half_width, -half_height, half_height, znear, zfar, vertices);
}

void GetFrustumVertices(float l, float r, float b, float t, float znear, float zfar,  Vector3 vertices[8])
{
    float z_ratio = zfar / znear;

    vertices[0] = Vector3(l, b, -znear);
    vertices[1] = Vector3(r, b, -znear);
    vertices[2] = Vector3(l, t, -znear);
    vertices[3] = Vector3(r, t, -znear);

    vertices[4] = vertices[0] * z_ratio;
    vertices[5] = vertices[1] * z_ratio;
    vertices[6] = vertices[2] * z_ratio;
    vertices[7] = vertices[3] * z_ratio;
}
