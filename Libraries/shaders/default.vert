#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aTex;

out vec2 texCoord;

uniform float scale;

void main()
{
   gl_Position = vec4(aPos.x*scale, aPos.y*scale, 0*scale, 1.0);

   texCoord = aTex;

}