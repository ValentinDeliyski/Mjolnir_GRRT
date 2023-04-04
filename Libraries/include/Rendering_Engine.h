#pragma once
#ifndef RENDERING_ENGINE
    
    #define RENDERING_ENGINE

    #include <glad/glad.h>
    #include <GLFW/glfw3.h>
    #include <string>
    #include <fstream>
    #include <sstream>
    #include "Constants.h"

    GLFWwindow* OpenGL_init(double aspect_ratio);

    void update_rendering_window(GLFWwindow* window, double aspect_ratio);

    void set_pixel_color(double intensity, int pixel_count);

    void set_background_pattern_color(double State_vector[], double old_state[], int texture_indexer, double J);

    void normalize_colormap(Initial_conditions_type* s_Initial_Conditions);

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

    class Window_Callbacks {

    public:

        static void define_button_callbacks(GLFWwindow* window, int key, int scancode, int action, int mods);

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