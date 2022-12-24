#include "Rendering_Engine.h"
#include "Constants.h"
#include <iostream>

extern float Max_Intensity;
extern float texture_buffer[];

GLFWwindow* OpenGL_init() {

    // Initialize GLFW
    glfwInit();
    // Tell GLFW what version of OpenGL we are using 
    // In this case we are using OpenGL 3.3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    // Tell GLFW we are using the CORE profile
    // So that means we only have the modern functions
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // Create a GLFWwindow object of 800 by 800 pixels, naming it "Gravitational Lenser"
    GLFWwindow* window = glfwCreateWindow(800, 800, "Gravitational Lenser", NULL, NULL);
    // Introduce the window into the current context
    glfwMakeContextCurrent(window);
    // Turn off vsync because if slows down the simulation A LOT
    glfwSwapInterval(0);
    //Load GLAD so it configures OpenGL
    gladLoadGL();
    // Specify the viewport of OpenGL in the Window
    // In this case the viewport goes from x = 0, y = 0, to x = 800, y = 800
    glViewport(0, 0, 800, 800);

    return window;

}

void set_pixel_color(double Intensity, int texture_indexer) {

    Max_Intensity = 2.8051e-05;

    float x = Intensity / Max_Intensity;

    // Red Channel

    float R = 1.0 / 0.5 * x;

    if (R > 1) { R = 1; }

    // Blue Channel

    float G{};

    if (x > 0.5) { G = 1.0 / 0.5 * x - 1.0; }

    if (G > 1) { G = 1; }

    // Green Channel

    float B{};

    if (x > 0.75) { B = 1.0 / 0.25 * x - 0.75 / 0.25; }

    if (B > 1) { B = 1; }

    texture_buffer[texture_indexer] = R;
    texture_buffer[texture_indexer + 1] = G;
    texture_buffer[texture_indexer + 2] = B;

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
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, vertices, GL_STATIC_DRAW);

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

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_NEAREST);
    
    return  texture;
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