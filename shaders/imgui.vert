#version 330 core

uniform mat4 u_ProjectionMatrix;

layout (location=0) in vec4 in_Position;
layout (location=1) in vec2 in_UV;
layout (location=2) in vec4 in_Color;

out Data
{
    vec2 UV;
    vec4 Color;
} vs_out;

void main()
{
    gl_Position = u_ProjectionMatrix * in_Position;
    vs_out.UV = in_UV;
    vs_out.Color = in_Color;
}
