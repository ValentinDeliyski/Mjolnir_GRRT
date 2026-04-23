#define _USE_MATH_DEFINES

#include "Rendering_Engine.h"
#include "lensing.h"
#include "General_GR_functions.h"
#include "Sim_Modes.h"
#include <iostream>

void Rendering_engine::OpenGL_init(Initial_conditions_type* p_Init_Conditions) {

    // Allocate the intensity and texture buffers

    this->ray_number_x = p_Init_Conditions->Observer_params.resolution_x;
    this->ray_number_y = p_Init_Conditions->Observer_params.resolution_y;

    this->Intensity_buffer = (float*)calloc(static_cast<size_t>(this->ray_number_x) * static_cast<size_t>(this->ray_number_y), sizeof(float));
    this->texture_buffer   = (float*)calloc(static_cast<size_t>(this->ray_number_x) * static_cast<size_t>(this->ray_number_y) * 3, sizeof(float));

    // Calculate the aspect ratio of the rendering window

    float Y_angle_max = float(p_Init_Conditions->Observer_params.y_angle_max);
    float Y_angle_min = float(p_Init_Conditions->Observer_params.y_angle_min);
    float X_angle_max = float(p_Init_Conditions->Observer_params.x_angle_max);
    float X_angle_min = float(p_Init_Conditions->Observer_params.x_angle_min);

    float aspect_ratio = (Y_angle_max - Y_angle_min) / (X_angle_max - X_angle_min);
    
    if (!p_Init_Conditions->Observer_params.Use_angular_coords) {

        Y_angle_max = float(atan2(p_Init_Conditions->Observer_params.y_max, p_Init_Conditions->Observer_params.distance));
        Y_angle_min = float(atan2(p_Init_Conditions->Observer_params.y_min, p_Init_Conditions->Observer_params.distance));
        X_angle_max = float(atan2(p_Init_Conditions->Observer_params.x_max, p_Init_Conditions->Observer_params.distance));
        X_angle_min = float(atan2(p_Init_Conditions->Observer_params.x_min, p_Init_Conditions->Observer_params.distance));

        aspect_ratio = (Y_angle_max - Y_angle_min) / (X_angle_max - X_angle_min);

    }
    
    // Initialize GLFW
    glfwInit();   

    // Tell GLFW what version of OpenGL I am using -> OpenGL 4.6
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);

    // Tell GLFW we are using the CORE profile -> we only have the modern functions
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    this->window = glfwCreateWindow(1024, int(aspect_ratio * 1024), "Mjolnir GRRT", NULL, NULL);

    // Introduce the window into the current context
    glfwMakeContextCurrent(this->window);

    // Turn off vsync
    glfwSwapInterval(0);

    // Load GLAD so it configures OpenGL
    gladLoadGL();

    // Specify the viewport of OpenGL in the Window -> x = [0, aspect_ratio * 1024], y = [0, 1024]
    glViewport(0, 0, 1024, int(aspect_ratio * 1024));

    // The simulation image is interpreted as a texture
    init_texture();

    // This thing (after linkning) combines the bottom two things into one object
    // NEEDS TO BE BEFORE THE VERTEX BUFFER AND ELEMENT BUFFER CALLS
    Vertex_array Vertex_array;
    Vertex_array.Bind();

    // This thing holds the edges of the triangles that the renderer draws
    Vertex_Buffer Vertex_buffer(this->vertices, sizeof(this->vertices));
    // This thing holds the sequence in which the edges should be connected
    Element_Buffer Element_buffer(this->Vertex_order, sizeof(this->Vertex_order));

    /* Tell openGL that the data inside the vertex buffer at index "0" has "2" components of size "4 * sizeof(float)" and are offset in the array by 0.
       These "indecies" are specified in the vertex shader - I tell it there that index (called position there) 0 are the vertex coordinates. */
    Vertex_array.Linkattrib(0, 2, 4 * sizeof(float), (const void*)0);

    /* Tell openGL that the data inside the vertex buffer at index "1" has "2" components of size "4 * sizeof(float)" and are offset in the array by 2 floats.
       These "indecies" are specified in the vertex shader - I tell it there that index (called position there) 1 are the texture coordinates. */
    Vertex_array.Linkattrib(1, 2, 4 * sizeof(float), (const void*)(2 * sizeof(float)));

    // Generates a Shader object using the shaders defualt.vert and default.frag
    Shader shaderProgram(static_cast<const char*>(p_Init_Conditions->File_manager_params.Vert_shader_path.c_str()),
                         static_cast<const char*>(p_Init_Conditions->File_manager_params.Frag_shader_path.c_str()));

    shaderProgram.Activate();

    // Returns an integer handle for the variable "scale" inside the shader program.
    GLuint Scale_factor_handle = glGetUniformLocation(shaderProgram.ID, "scale");

    // Returns an integer handle for the variable "u_Texture" inside the shader program.
    GLuint Texture_uniform_handle = glGetUniformLocation(shaderProgram.ID, "u_Texture");

    /* Sets the value of the texture uniform to the index of the binded texture (0 in this case). 
       This tells openGL which texture to sample. */
    glUniform1i(Texture_uniform_handle, 0);

    // Set the value of the the scaler to 1.5f
    glUniform1f(Scale_factor_handle, 1.5f);

}

void Rendering_engine::update_rendering_window() const {

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, this->ray_number_x, this->ray_number_y, 0, GL_RGB, GL_FLOAT, this->texture_buffer);
    // Specify the color of the background
    glClearColor(0.07f, 0.13f, 0.17f, 1.0f);
    // Clean the back buffer and assign the new color to it
    glClear(GL_COLOR_BUFFER_BIT);
    // Draw primitives, number of indices, datatype of indices, index of indices
    glDrawElements(GL_TRIANGLES, sizeof(this->Vertex_order) / sizeof(float), GL_UNSIGNED_INT, 0);
    // Swap the back buffer with the front buffer
    glfwSwapBuffers(this->window);
    // Take care of all GLFW events
    glfwPollEvents();

}

void Rendering_engine::set_pixel_color(float Intensity, int texture_idx) {

    float x = Intensity / this->Max_Intensity;

    // Red Channel

    float R = 1.0f / 0.5f * x;

    if (R > 1.0f) { R = 1.0f; }

    // Blue Channel

    float G{};

    if (x > 0.5f) { G = 1.0f / 0.5f * x - 1.0f; }

    if (G > 1.0f) { G = 1.0f; }

    // Green Channel

    float B{};

    if (x > 0.75f) { B = 1.0f / 0.25f * x - 0.75f / 0.25f; }

    if (B > 1.0f) { B = 1.0f; }

    texture_buffer[texture_idx + 0] = R;
    texture_buffer[texture_idx + 1] = G;
    texture_buffer[texture_idx + 2] = B;

}

void Rendering_engine::renormalize_colormap() {

    // TODO: think about simplifying this

    float Current_max{};

    for (int index = 0; index <= this->texture_indexer; index += 3) {

        if (this->Intensity_buffer[int(index / 3)] > Current_max) {

            Current_max = this->Intensity_buffer[int(index / 3)];

        }

    }

    this->Max_Intensity = Current_max;

    for (int index = 0; index <= this->texture_indexer - 1; index += 3) {

        set_pixel_color(Intensity_buffer[int(index / 3)], index);

    }

}

/*******************************************
|                                          |
| Vertex Buffer Class Function Definitions |
|                                          |
*******************************************/

Rendering_engine::Vertex_Buffer::Vertex_Buffer(const GLfloat* Vertex_attribute_data, GLsizeiptr size) {

    /* Calls openGL to (internally) allocate one buffer and store its handle in this->Vertex_buffer_ID. */
    glGenBuffers(1, &this->Vertex_buffer_ID);

    /* Calls openGL to set the above buffer as an active ARRAY_BUFFER (a.e. one that holds vertex attirubtes). */
    glBindBuffer(GL_ARRAY_BUFFER, this->Vertex_buffer_ID);

    /* Calls openGL to store the data that the pointer "Vertex_attribute_data" points to, inside the now active buffer with ID, this->Vertex_buffer_ID. */
    glBufferData(GL_ARRAY_BUFFER, size, Vertex_attribute_data, GL_DYNAMIC_DRAW);

}

void Rendering_engine::Vertex_Buffer::Bind() const { glBindBuffer(GL_ARRAY_BUFFER, this->Vertex_buffer_ID); }

void Rendering_engine::Vertex_Buffer::Unbind() const { glBindBuffer(GL_ARRAY_BUFFER, 0); }

void Rendering_engine::Vertex_Buffer::Delete() const { glDeleteBuffers(1, &this->Vertex_buffer_ID); }

/******************************************
|                                         |
| Vertex Array Class Function Definitions |
|                                         |
******************************************/

Rendering_engine::Vertex_array::Vertex_array()
{
    glGenVertexArrays(1, &this->ID);
}

void Rendering_engine::Vertex_array::Linkattrib(GLuint index, GLuint numComponents, GLsizei stride, const void* offset) {

    /* Tell openGL that the Vertex attirbutes stored in a buffer with ID "index" and have "cumComponents" many components of type "GL_FLOAT".
       Further tell openGL that consecutive atributes are speperated by "stride" array elements, and the offset of the first component is "offset". */
    glVertexAttribPointer(index, numComponents, GL_FLOAT, GL_FALSE, stride, offset);
    // Enable the Vertex Attribute so that OpenGL knows to use it
    glEnableVertexAttribArray(index);

}

void Rendering_engine::Vertex_array::Bind() const {

    glBindVertexArray(this->ID);

}

void Rendering_engine::Vertex_array::Unbind() const {

    glBindVertexArray(0);

}

void Rendering_engine::Vertex_array::Delete() const {


    glDeleteVertexArrays(1, &this->ID);

}

/********************************************
|                                           |
| Element Buffer Class Function Definitions |
|                                           |
********************************************/

Rendering_engine::Element_Buffer::Element_Buffer(const GLuint* vertices, GLsizeiptr size) {

    glGenBuffers(1, &this->Element_buffer_ID);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->Element_buffer_ID);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, vertices, GL_DYNAMIC_DRAW);

}

void Rendering_engine::Element_Buffer::Bind() const { glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->Element_buffer_ID); }

void Rendering_engine::Element_Buffer::Unbind() const { glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); }

void Rendering_engine::Element_Buffer::Delete() const { glDeleteBuffers(1, &this->Element_buffer_ID); }

GLuint Rendering_engine::init_texture() {

    /* Declares a texture hangle (uint) that numbers the different textures. */
    GLuint texture_handle;

    /* Calls openGL to (internally) generate a texture and store its label in the "texture_handle" variable. */
    glGenTextures(1, &texture_handle);

    /* Calls openGL to set the current ative texture to be the one in its internal "slot 0". */
    glActiveTexture(GL_TEXTURE0);

    /* Calls openGL to allocate the active "slot 0" texture to the one generated above. */
    glBindTexture(GL_TEXTURE_2D, texture_handle);

    /* Calls openGL to set the pixel interpolation to be linear. */
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    /* Calls openGL to set the image to clamp to the edge, rather than repeat.
       In theory this shouldn't matter for my use case? */
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    
    return texture_handle;
}

std::string Rendering_engine::get_file_contents(const char* filename, std::string file_type)
{
    std::ifstream in(filename, std::ios::binary);

    if (in)
    {
        std::string contents;
        in.seekg(0, std::ios::end);
        contents.resize(in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(&contents[0], contents.size());
        in.close();

        return contents;

    }
    else {

        std::cout << "ERROR! Could not parse the " << file_type << ". Check the file path.\n";

        exit(ERROR);

    }
 
}

Rendering_engine::Shader::Shader(const char* vertexFile, const char* fragmentFile)
{

    std::string vertexCode = get_file_contents(vertexFile, "Vertex Shader");
    std::string fragmentCode = get_file_contents(fragmentFile, "Fragment Shader");

    const char* vertexSource = vertexCode.c_str();
    const char* fragmentSource = fragmentCode.c_str();

    // Create Vertex Shader Object and get its reference
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    // Attach Vertex Shader source to the Vertex Shader Object
    glShaderSource(vertexShader, 1, &vertexSource, NULL);
    // Compile the Vertex Shader into machine code
    glCompileShader(vertexShader);

    // Create Fragment Shader Object and get its reference
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    // Attach Fragment Shader source to the Fragment Shader Object
    glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
    // Compile the Vertex Shader into machine code
    glCompileShader(fragmentShader);

    // Create Shader Program Object and get its reference
    this->ID = glCreateProgram();
    // Attach the Vertex and Fragment Shaders to the Shader Program
    glAttachShader(this->ID, vertexShader);
    glAttachShader(this->ID, fragmentShader);
    // Wrap-up/Link all the shaders together into the Shader Program
    glLinkProgram(this->ID);

    // Delete the now useless Vertex and Fragment Shader objects
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

}

void Rendering_engine::Shader::Activate() const {

    glUseProgram(this->ID);

}

void Rendering_engine::Shader::Delete() const {

    glDeleteProgram(this->ID);

}

// GSL expects the first argument to be a window handle, the second argument to be an int "scancode" and the final to be an int "mods". I dont use them, so I do not label them
void Rendering_engine::Window_Callbacks::define_button_callbacks(GLFWwindow*, int key, int, int action, int) {
    
    if (key == GLFW_KEY_UP && action == GLFW_PRESS) {

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    }

    if (key == GLFW_KEY_DOWN && action == GLFW_PRESS) {

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    }

}

void Rendering_engine::Free_memory() const {

    free(this->Intensity_buffer);
    free(this->texture_buffer);

}
