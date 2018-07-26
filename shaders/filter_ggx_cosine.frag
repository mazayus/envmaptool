#version 330 core

#define PI 3.14159265

uniform samplerCube u_EnvMap;

uniform int u_NumSamples;
uniform int u_StartSample;
uniform int u_EndSample;

uniform float u_Roughness;

uniform mat3 u_FaceBasis;

in vec3 Position;

out vec4 out_Integral;

float RadicalInverse(uint bits)
{
    bits = (bits << 16) | (bits >> 16);
    bits = ((bits & 0x55555555u) << 1) | ((bits & 0xAAAAAAAAu) >> 1);
    bits = ((bits & 0x33333333u) << 2) | ((bits & 0xCCCCCCCCu) >> 2);
    bits = ((bits & 0x0F0F0F0Fu) << 4) | ((bits & 0xF0F0F0F0u) >> 4);
    bits = ((bits & 0x00FF00FFu) << 8) | ((bits & 0xFF00FF00u) >> 8);

    const float InverseMaxU32 = 2.3283064365386963e-10;
    return float(bits) * InverseMaxU32;
}

vec2 Hammersley(uint i, uint N)
{
    return vec2(float(i) / float(N), RadicalInverse(i));
}

vec3 FilterEnvMap(in mat3 TBN)
{
    vec3 Integral = vec3(0, 0, 0);

    uint NumSamples = uint(u_NumSamples);
    uint StartSample = uint(u_StartSample);
    uint EndSample = uint(u_EndSample);

    for (uint i = StartSample; i < EndSample; ++i)
    {
        vec2 Xi = Hammersley(i, NumSamples);

        float CosTheta = sqrt((1 - Xi.x) / (1 + Xi.x * (u_Roughness * u_Roughness - 1)));
        float SinTheta = sqrt(1 - CosTheta * CosTheta);

        float CosPhi = cos(2 * PI * Xi.y);
        float SinPhi = sin(2 * PI * Xi.y);

        vec3 l = TBN * vec3(CosPhi * SinTheta, SinPhi * SinTheta, CosTheta);
        vec3 Radiance = texture(u_EnvMap, l).rgb;

        Integral += Radiance;
    }

    return Integral / float(NumSamples);
}

void main()
{
    vec3 Normal = u_FaceBasis * normalize(Position);
    vec3 Tangent = normalize(cross(vec3(0, 1, 0), Normal));
    vec3 Bitangent = cross(Normal, Tangent);

    mat3 TBN = mat3(Tangent, Bitangent, Normal);

    out_Integral = vec4(FilterEnvMap(TBN), 1);
}
