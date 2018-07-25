#version 330 core

uniform sampler2D u_Font;

in Data
{
    vec2 UV;
    vec4 Color;
} fs_in;

out vec4 out_Color;

void main()
{
    out_Color = fs_in.Color * texture(u_Font, fs_in.UV);
}
