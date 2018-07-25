#version 330 core

out vec3 Position;

void main()
{
    const vec2 Vertices[4] = vec2[4](vec2(-1, -1), vec2(1, -1), vec2(-1, 1), vec2(1, 1));

    Position = vec3(Vertices[gl_VertexID], 1);
    gl_Position = vec4(Vertices[gl_VertexID], 0, 1);
}
