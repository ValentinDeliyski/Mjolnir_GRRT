#pragma once
#ifndef RENDERING_ENGINE
	
	#define RENDERING_ENGINE

	#include <glad/glad.h>
	#include <GLFW/glfw3.h>
	#include <string>
	#include <fstream>
	#include <sstream>
	#include "Constants.h"

	GLFWwindow* OpenGL_init();

	void set_pixel_color(double intensity, int pixel_count);

	/***************************************
	|									   |
	| Classes that abstract the OpenGL api |
	|                                      |
	***************************************/


	class Vertex_Buffer {

	public:

		GLuint ID;
		Vertex_Buffer(const GLfloat* verticies, GLsizeiptr size);

		void Bind();
		void Unbind();
		void Delete();

	};

	class Vertex_array {

	public:

		GLuint ID;
		Vertex_array();

		void Linkattrib(Vertex_Buffer Vertex_Buffer, GLuint layout, GLuint numComponents, GLenum type, GLsizei stride, void* offset);
		void Bind();
		void Unbind();
		void Delete();

	};

	class Element_Buffer {

	public:
		GLuint ID;
		Element_Buffer(const GLuint* verticies, GLsizeiptr size);

		void Bind();
		void Unbind();
		void Delete();

	};

	class Shader {

	public:

		GLuint ID;
		Shader(const char* vertexFile, const char* fragmentFile);

		void Activate();
		void Delete();

	};

	std::string get_file_contents(const char* filename);

	GLuint init_texture();

	const GLfloat vertices[] =
	{ 
   // Vertex coordinates  Texture Coordinates
		-0.5f, -0.5f,		  0.0f, 0.0f,     // Lower left corner
		-0.5f,  0.5f,		  0.0f, 1.0f,     // Upper left corner
		 0.5f,  0.5f,		  1.0f, 1.0f,     // Upper right corner
		 0.5f, -0.5f,		  1.0f, 0.0f,     // Lower right corner
	};

	const GLuint Vertex_order[] = {0, 2, 1, 0, 3, 2};

#endif