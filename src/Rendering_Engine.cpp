#define _USE_MATH_DEFINES

#include "Rendering_Engine.h"
#include "lensing.h"
#include "General_GR_functions.h"
#include "General_math_functions.h"
#include "Sim_Modes.h"
#include <iostream>

void Rendering_engine::OpenGL_init() {

    // Initialize GLFW
    glfwInit();

    // Tell GLFW what version of OpenGL we are using -> OpenGL 3.3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

    // Tell GLFW we are using the CORE profile -> we only have the modern functions
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow(1200, 1200, "Gravitational Ray Tracer", NULL, NULL);

    // Introduce the window into the current context
    glfwMakeContextCurrent(window);

    // Turn off vsync because if slows down the simulation A LOT
    glfwSwapInterval(0);

    // Load GLAD so it configures OpenGL
    gladLoadGL();

    // Specify the viewport of OpenGL in the Window -> x = [0, aspect_ratio * 1200], y = [0, 1200]
    glViewport(0, 0, aspect_ratio * 1200, 1200);

    // The simulation image is interpreted as a texture
    GLuint texture = init_texture();

    // This thing (after linkning) combines the bottom two things into one object
    // NEEDS TO BE BEFORE THE VERTEX BUFFER AND ELEMENT BUFFER CALLS
    Vertex_array Vertex_array;
    Vertex_array.Bind();

    // This thing holds the edges of the triangles that the renderer draws
    Vertex_Buffer Vertex_buffer(vertices, sizeof(vertices));
    // This thing holds the sequence in which the edges should be connected
    Element_Buffer Element_buffer(Vertex_order, sizeof(Vertex_order));

    Vertex_array.Linkattrib(Vertex_buffer, 0, 2, GL_FLOAT, 4 * sizeof(float), (void*)0);
    Vertex_array.Linkattrib(Vertex_buffer, 1, 2, GL_FLOAT, 4 * sizeof(float), (void*)(2 * sizeof(float)));

    // Generates a Shader object using the shaders defualt.vert and default.frag
    Shader shaderProgram(vert_shader_path, frag_shader_path);
    shaderProgram.Activate();

    // Generates a float (with an int ID), that scales the output image
    GLuint uniID = glGetUniformLocation(shaderProgram.ID, "scale");

    // Generates an int (with an int ID), that tells the shader *insert what it tells it here*
    GLuint tex0Uni = glGetUniformLocation(shaderProgram.ID, "tex0");
    glUniform1i(tex0Uni, 0);

    // Activates the scaler with a value of 1.5f
    glUniform1f(uniID, 1.5f);
    // Binds the texture array (RGB values 
    glBindTexture(GL_TEXTURE_2D, texture);

}

void Rendering_engine::update_rendering_window() {

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, RESOLUTION, RESOLUTION * aspect_ratio, 0, GL_RGB, GL_FLOAT, texture_buffer);
    // Specify the color of the background
    glClearColor(0.07f, 0.13f, 0.17f, 1.0f);
    // Clean the back buffer and assign the new color to it
    glClear(GL_COLOR_BUFFER_BIT);
    // Draw primitives, number of indices, datatype of indices, index of indices
    glDrawElements(GL_TRIANGLES, sizeof(Vertex_order) / sizeof(float), GL_UNSIGNED_INT, 0);
    // Swap the back buffer with the front buffer
    glfwSwapBuffers(window);
    // Take care of all GLFW events
    glfwPollEvents();

}

void Rendering_engine::update_max_intensity(float Intensity) {

    if (Intensity > Max_Intensity) {

        Max_Intensity = Intensity;
        renormalize_colormap_flag = true;

    }

}

void Rendering_engine::set_pixel_color(float Intensity, int texture_indexer) {

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

void Rendering_engine::set_background_pattern_color(double State_vector[], double old_state[], int texture_indexer, double J) {

    double theta = (State_vector[e_theta] + old_state[e_theta]) / 2;
    double phi = (State_vector[e_phi] + old_state[e_phi]) / 2;

    if (J*J < 1e-5) {

        phi = phi + M_PI_2;

    }

    double grayscale_value = pow((1 + sin(10 * phi) * sin(10 * theta)) / 2, 1.0 / 5);



    texture_buffer[texture_indexer]     = grayscale_value;
    texture_buffer[texture_indexer + 1] = grayscale_value;
    texture_buffer[texture_indexer + 2] = grayscale_value;

}

void Rendering_engine::normalize_colormap(Initial_conditions_type* s_Initial_Conditions) {

    /*

    Do one scan line in the middle of the image to find the maximum intensity for use in the colormap

    */

    std::cout << "Initial y = 0 scan to normalize the colormap..." << '\n';

    Results_type Colormap_normalization_results{}; 

    for (int pixel_num = 0; pixel_num <= RESOLUTION - 1; pixel_num++) {

        get_intitial_conditions_from_angles(s_Initial_Conditions, 0, H_angle_max - pixel_num * Scan_Step);

        Colormap_normalization_results = Propagate_ray(s_Initial_Conditions);

        double Intensity_current = Colormap_normalization_results.Intensity[direct] +
                                   Colormap_normalization_results.Intensity[first]  +
                                   Colormap_normalization_results.Intensity[second] +
                                   Colormap_normalization_results.Intensity[third];

        if (Intensity_current > Max_Intensity) {

            Max_Intensity = Intensity_current;

        }

        print_progress(pixel_num, RESOLUTION - 1, false, true);

    }

}

void Rendering_engine::renormalize_colormap() {

    if (renormalize_colormap_flag == true) {

        for (int index = 0; index <= texture_indexer - 1; index += 3) {


            set_pixel_color(Intensity_buffer[int(index / 3)], index);


        }

        renormalize_colormap_flag = false;

    }
    else {

        for (int index = texture_indexer - (RESOLUTION) * 3; index <= texture_indexer - 1; index += 3) {

            set_pixel_color(Intensity_buffer[int(index / 3)], index);

        }

    }

}

/*******************************************
|                                          |
| Vertex Buffer Class Function Definitions |
|                                          |
*******************************************/

Rendering_engine::Vertex_Buffer::Vertex_Buffer(const GLfloat* vertices, GLsizeiptr size) {

    glGenBuffers(1, &ID);
    glBindBuffer(GL_ARRAY_BUFFER, ID);
    glBufferData(GL_ARRAY_BUFFER, size, vertices, GL_DYNAMIC_DRAW);

}

void Rendering_engine::Vertex_Buffer::Bind()
{

    glBindBuffer(GL_ARRAY_BUFFER, ID);

}

void Rendering_engine::Vertex_Buffer::Unbind()
{

    glBindBuffer(GL_ARRAY_BUFFER, 0);

}

void Rendering_engine::Vertex_Buffer::Delete()
{

    glDeleteBuffers(1, &ID);

}

/******************************************
|                                         |
| Vertex Array Class Function Definitions |
|                                         |
******************************************/

Rendering_engine::Vertex_array::Vertex_array()
{
    glGenVertexArrays(1, &ID);
}

void Rendering_engine::Vertex_array::Linkattrib(Vertex_Buffer Vertex_Buffer, GLuint index, GLuint numComponents, GLenum type, GLsizei stride, void* offset) {


    Vertex_Buffer.Bind();
    glVertexAttribPointer(index, numComponents, type, GL_FALSE, stride, offset);
    // Enable the Vertex Attribute so that OpenGL knows to use it
    glEnableVertexAttribArray(index);
    Vertex_Buffer.Unbind();

}

void Rendering_engine::Vertex_array::Bind() {

    glBindVertexArray(ID);

}

void Rendering_engine::Vertex_array::Unbind() {

    glBindVertexArray(0);

}

void Rendering_engine::Vertex_array::Delete() {


    glDeleteVertexArrays(1, &ID);

}

/********************************************
|                                           |
| Element Buffer Class Function Definitions |
|                                           |
********************************************/

Rendering_engine::Element_Buffer::Element_Buffer(const GLuint* vertices, GLsizeiptr size) {

    glGenBuffers(1, &ID);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ID);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, vertices, GL_DYNAMIC_DRAW);

}

void Rendering_engine::Element_Buffer::Bind()
{

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ID);

}

void Rendering_engine::Element_Buffer::Unbind()
{

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

}

void Rendering_engine::Element_Buffer::Delete()
{

    glDeleteBuffers(1, &ID);

}

GLuint Rendering_engine::init_texture() {

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

std::string Rendering_engine::get_file_contents(const char* filename)
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

Rendering_engine::Shader::Shader(const char* vertexFile, const char* fragmentFile)
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

void Rendering_engine::Shader::Activate() {

    glUseProgram(ID);
    

}

void Rendering_engine::Shader::Delete() {

    glDeleteProgram(ID);

}

void Rendering_engine::Window_Callbacks::define_button_callbacks(GLFWwindow* window, int key, int scancode, int action, int mods) {
    
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
