#pragma once
#ifndef RENDERING_ENGINE
    
    #define RENDERING_ENGINE

    #include <glad/glad.h>
    #include <GLFW/glfw3.h>
    #include <string>
    #include <fstream>
    #include <sstream>
    #include "Constants.h"

    class Rendering_engine {

        public:

            GLFWwindow* window{};

            bool renormalize_colormap_flag = false;

            int texture_indexer{};
            
            float Max_Intensity{};
            float Intensity_buffer[NUM_RAYS]{};
            float texture_buffer[NUM_RAYS * 3]{};

            float aspect_ratio = H_angle_max / V_angle_max;
            

            const GLuint Vertex_order[6] = { 0, 2, 1, 0, 3, 2 };
            const GLfloat vertices[16] =
            {
                // Vertex coordinates  Texture Coordinates
                     -0.5f, -0.5f,		  0.0f, 0.0f,     // Lower left corner
                     -0.5f,  0.5f,		  0.0f, 1.0f,     // Upper left corner
                      0.5f,  0.5f,		  1.0f, 1.0f,     // Upper right corner
                      0.5f, -0.5f,		  1.0f, 0.0f,     // Lower right corner
            };

            static std::string get_file_contents(const char* filename);

            void OpenGL_init();

            void update_rendering_window();

            void set_pixel_color(float intensity, int pixel_count);

            void renormalize_colormap();

            void update_max_intensity(float Intensity);

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

            GLuint init_texture();

    };

#endif