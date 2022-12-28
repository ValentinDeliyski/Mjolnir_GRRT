#define _USE_MATH_DEFINES

#include "Rendering_Engine.h"
#include "Constants.h"
#include "Enumerations.h"
#include <cmath>
#include <iostream>


extern float Max_Intensity;
extern float texture_buffer[];

GLFWwindow* OpenGL_init(double aspect_ratio) {

    // Initialize GLFW
    glfwInit();
    // Tell GLFW what version of OpenGL we are using 
    // In this case we are using OpenGL 3.3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    // Tell GLFW we are using the CORE profile
    // So that means we only have the modern functions
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // Create a GLFWwindow object of 800 by 800 pixels, naming it "Gravitational Ray Tracer"
    GLFWwindow* window = glfwCreateWindow(1200, 1200, "Gravitational Ray Tracer", NULL, NULL);
    // Introduce the window into the current context
    glfwMakeContextCurrent(window);
    // Turn off vsync because if slows down the simulation A LOT
    glfwSwapInterval(0);
    //Load GLAD so it configures OpenGL
    gladLoadGL();
    // Specify the viewport of OpenGL in the Window
    // In this case the viewport goes from x = 0, y = 0, to x = 800, y = 800
    glViewport(0, 0, 1200, aspect_ratio * 1200);

    return window;

}

void set_pixel_color(double Intensity, int texture_indexer) {

    float x = Intensity / Max_Intensity;

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

    texture_buffer[texture_indexer + 0] = R;
    texture_buffer[texture_indexer + 1] = G;
    texture_buffer[texture_indexer + 2] = B;

}

void set_background_pattern_color(double State_vector[], double old_state[], int texture_indexer, double J) {

    double theta = (State_vector[e_theta] + old_state[e_theta]) / 2;
    double phi = (State_vector[e_phi] + old_state[e_phi]) / 2;

    if (J*J < 1e-5) {

        phi = phi + M_PI_2;

    }

    double grayscale_value = pow((1 + sin(10 * phi) * sin(10 * theta)) / 2, 1.0 / 5);



    texture_buffer[texture_indexer] = grayscale_value;
    texture_buffer[texture_indexer + 1] = grayscale_value;
    texture_buffer[texture_indexer + 2] = grayscale_value;

}

/*******************************************
|                                          |
| Vertex Buffer Class Function Definitions |
|                                          |
*******************************************/

Vertex_Buffer::Vertex_Buffer(const GLfloat* vertices, GLsizeiptr size) {

    glGenBuffers(1, &ID);
    glBindBuffer(GL_ARRAY_BUFFER, ID);
    glBufferData(GL_ARRAY_BUFFER, size, vertices, GL_DYNAMIC_DRAW);

}

void Vertex_Buffer::Bind()
{

    glBindBuffer(GL_ARRAY_BUFFER, ID);

}

void Vertex_Buffer::Unbind()
{

    glBindBuffer(GL_ARRAY_BUFFER, 0);

}

void Vertex_Buffer::Delete()
{

    glDeleteBuffers(1, &ID);

}

/******************************************
|                                         |
| Vertex Array Class Function Definitions |
|                                         |
******************************************/

Vertex_array::Vertex_array()
{
    glGenVertexArrays(1, &ID);
}

void Vertex_array::Linkattrib(Vertex_Buffer Vertex_Buffer, GLuint index, GLuint numComponents, GLenum type, GLsizei stride, void* offset) {


    Vertex_Buffer.Bind();
    glVertexAttribPointer(index, numComponents, type, GL_FALSE, stride, offset);
    // Enable the Vertex Attribute so that OpenGL knows to use it
    glEnableVertexAttribArray(index);
    Vertex_Buffer.Unbind();

}

void Vertex_array::Bind() {

    glBindVertexArray(ID);

}

void Vertex_array::Unbind() {

    glBindVertexArray(0);

}

void Vertex_array::Delete() {


    glDeleteVertexArrays(1, &ID);

}

/********************************************
|                                           |
| Element Buffer Class Function Definitions |
|                                           |
********************************************/

Element_Buffer::Element_Buffer(const GLuint* vertices, GLsizeiptr size) {

    glGenBuffers(1, &ID);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ID);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, vertices, GL_DYNAMIC_DRAW);

}

void Element_Buffer::Bind()
{

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ID);

}

void Element_Buffer::Unbind()
{

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

}

void Element_Buffer::Delete()
{

    glDeleteBuffers(1, &ID);

}

GLuint init_texture() {

    // Texture

    GLuint texture;
    glGenTextures(1, &texture);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texture);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_LINEAR);
    
    return texture;
}

std::string get_file_contents(const char* filename)
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
        return(contents);
    }
    throw(errno);
}

Shader::Shader(const char* vertexFile, const char* fragmentFile)
{

    std::string vertexCode = get_file_contents(vertexFile);
    std::string fragmentCode = get_file_contents(fragmentFile);

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
    ID = glCreateProgram();
    // Attach the Vertex and Fragment Shaders to the Shader Program
    glAttachShader(ID, vertexShader);
    glAttachShader(ID, fragmentShader);
    // Wrap-up/Link all the shaders together into the Shader Program
    glLinkProgram(ID);

    // Delete the now useless Vertex and Fragment Shader objects
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

}

void Shader::Activate() {

    glUseProgram(ID);
    

}

void Shader::Delete() {

    glDeleteProgram(ID);

}

void Window_Callbacks::define_button_callbacks(GLFWwindow* window, int key, int scancode, int action, int mods) {
    
    if (key == GLFW_KEY_UP && action == GLFW_PRESS) {

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_NEAREST);

    }

    if (key == GLFW_KEY_DOWN && action == GLFW_PRESS) {

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_LINEAR);

    }


}
