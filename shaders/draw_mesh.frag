#version 330 core

uniform vec3 u_CameraPosition;

uniform int u_DisplayMode;

uniform samplerCube u_EnvMap;

uniform float u_MipLevel;
uniform float u_Exposure;

in vec4 Position;
in vec4 ClipPosition;

out vec4 out_Color;

void main()
{
    vec3 Normal = normalize(Position.xyz / Position.w);

    vec3 TexCoord = Normal;
    if (u_DisplayMode == 1)
    {
        vec3 V = u_CameraPosition - normalize(Position.xyz / Position.w);
        TexCoord = -reflect(V, Normal);
    }

    out_Color = textureLod(u_EnvMap, TexCoord, u_MipLevel);

    // exposure
    out_Color.rgb = out_Color.rgb * exp2(u_Exposure);

    // tone mapping
    out_Color.rgb = out_Color.rgb / (1 + out_Color.rgb);
}
