#version 330 core
layout (location = 0) in vec2 Position_array;
layout (location = 1) in vec2 Texture_array;

out vec2 Texture_coordinates;

uniform float scale;

void main()
{
   gl_Position = vec4(Position_array.x * scale, Position_array.y * scale, 0 * scale, 1.0);

   Texture_coordinates = Texture_array;

}