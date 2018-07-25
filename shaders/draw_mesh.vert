#version 330 core

uniform mat4 u_WorldToClipMatrix;
uniform mat4 u_ObjectToWorldMatrix;

layout (location=0) in vec4 in_Position;

out vec4 Position;
out vec4 ClipPosition;

void main()
{
    Position = u_ObjectToWorldMatrix * in_Position;
    ClipPosition = u_WorldToClipMatrix * Position;
    gl_Position = ClipPosition;
}
