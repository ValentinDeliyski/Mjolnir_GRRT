/****************************************************************************************************
|                                                                                                   |
|                          ---------  Gravitational Ray Tracer  ---------                           | 
|                                                                                                   |
|    * @Version: 3.5                                                                                |
|    * @Author: Valentin Deliyski                                                                   |
|    * @Description: This program numeriaclly integrates the equations of motion                    |
|    for null geodesics and radiative transfer in a curved spacetime,then projects                  |
|    them onto an observer's screen to construct relativistic images of accretion disks             |
|                                                                                                   |
|    * @Supported Spacetimes:                                                                       |
|        ** Kerr Black Holes                                                                        |
|        ** Static Regular Black Holes                                                              |
|        ** Rotating Traversable Wormholes                                                          |
|        ** Janis - Newman - Winicour Naked Singularities                                           |
|                                                                                                   |
|    * @Supported Disk Models                                                                       |
|        ** Novikov-Thorne                                                                          |
|        ** Generic Optically Thin Disk With Arbitrary Density, Emission and Absorbtion Profiles    |
|                                                                                                   |
****************************************************************************************************/

#define _USE_MATH_DEFINES

#include <string>
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "Spacetimes.h"
#include "Constants.h"
#include "Enumerations.h"
#include "IO_files.h"

#include "Disk_Models.h"
#include "General_GR_functions.h"
#include "Lensing.h"

#include "Rendering_Engine.h"

e_Spacetimes e_metric = Kerr;
e_Emission_model e_emission = Synchotron_phenomenological;

/* 

Define classes that hold the spacetime properites

*/

std::vector<c_Spacetime_Base*> Spacetimes = {

    new derived_Kerr_class(),
    new derived_RBH_class(),
    new derived_Wormhole_class(),
    new derived_JNW_class()
};

/*

Define the Observer class

*/

extern Real r_obs = 1e3;
extern Real theta_obs = 60. / 180 * M_PI;
Real phi_obs = 0;

c_Observer Observer_class(r_obs, theta_obs, phi_obs);

/*

Define the Optically Thin Disk Class

*/

Optically_Thin_Toroidal_Model OTT_Model;

/*

Define the Novikov-Thorne Disk Class

*/

Real r_in = Spacetimes[e_metric]->get_ISCO(Prograde);
Real r_out = 25;

Novikov_Thorne_Model NT_Model(r_in, r_out);

/*

Define some global boolians

*/

extern Const_bool lens_from_file = false;
extern Const_bool truncate = true;

/*

Rendering Engine variables

*/

int pixel_count{}; // This is sometimes off by one which is why Im literally counting them - prolly to due with the loops iterating over floats...
float Max_Intensity{};
int texture_indexer{};
bool Normalizing_colormap{};
float texture_buffer[RESOLUTION * 3]{};

void print_ASCII_art() {

    std::cout <<

        " ######   ########     ###    ##     ## #### ########    ###    ######## ####  #######  ##    ##    ###    ##          ########     ###    ##    ##    ######## ########     ###     ######  ######## ########  \n"
        "##    ##  ##     ##   ## ##   ##     ##  ##     ##      ## ##      ##     ##  ##     ## ###   ##   ## ##   ##          ##     ##   ## ##    ##  ##        ##    ##     ##   ## ##   ##    ## ##       ##     ## \n"
        "##        ##     ##  ##   ##  ##     ##  ##     ##     ##   ##     ##     ##  ##     ## ####  ##  ##   ##  ##          ##     ##  ##   ##    ####         ##    ##     ##  ##   ##  ##       ##       ##     ## \n"
        "##   #### ########  ##     ## ##     ##  ##     ##    ##     ##    ##     ##  ##     ## ## ## ## ##     ## ##          ########  ##     ##    ##          ##    ########  ##     ## ##       ######   ########  \n"
        "##    ##  ##   ##   #########  ##   ##   ##     ##    #########    ##     ##  ##     ## ##  #### ######### ##          ##   ##   #########    ##          ##    ##   ##   ######### ##       ##       ##   ##   \n"
        "##    ##  ##    ##  ##     ##   ## ##    ##     ##    ##     ##    ##     ##  ##     ## ##   ### ##     ## ##          ##    ##  ##     ##    ##          ##    ##    ##  ##     ## ##    ## ##       ##    ##  \n"
        " ######   ##     ## ##     ##    ###    ####    ##    ##     ##    ##    ####  #######  ##    ## ##     ## ########    ##     ## ##     ##    ##          ##    ##     ## ##     ##  ######  ######## ##     ## \n";

    std::cout << '\n' << '\n';

}

void print_progress(int current, int max, bool lens_from_file) {

    int current_digits = 1;

    if (current != 0) {

        current_digits = floor(log10f(current) + 1);

    }

    int max_digits = floor(log10f(max) + 1);

    if (current == 0) {

        if (lens_from_file) {

            std::cout << "Number Of Rays Cast: ";

        }
        else {

            std::cout << "Number Of Lines Scanned: ";

        }

        for (int i = 0; i <= max_digits + current_digits; i += 1) {

            std::cout << "0";

        }

    }

    for (int i = 0; i <= max_digits + current_digits; i += 1) {

        std::cout << "\b";

    }

    std::cout << current + 1 << "/" << max + 1;

}

double V_angle_min = -atan(15 / r_obs);
double V_angle_max = atan(15 / r_obs);

double H_angle_min = -atan(15 / r_obs);
double H_angle_max = atan(15 / r_obs);

void main() {
    
    auto start_time = std::chrono::high_resolution_clock::now();

    double Scan_Step = (H_angle_max - H_angle_min) / 128;

    GLFWwindow* window = OpenGL_init(V_angle_max / H_angle_max);

    //The simulation image is interpreted as a texture
    GLuint texture = init_texture();

    // This thing (after linkning) combines the bottom two things into one object
    // NEEDS TO BEFORE THE VERTEX BUFFER AND ELEMENT BUFFER CALLS
    Vertex_array Vertex_array;
    Vertex_array.Bind();

    // This thing holds the edges of the triangles that the renderer draws
    Vertex_Buffer Vertex_buffer(vertices, sizeof(vertices));
    // This thing holds the sequence in which the edges should be connected
    Element_Buffer Element_buffer(Vertex_order, sizeof(Vertex_order));

    Vertex_array.Linkattrib(Vertex_buffer, 0, 2, GL_FLOAT, 4 * sizeof(float), (void*)0);
    Vertex_array.Linkattrib(Vertex_buffer, 1, 2, GL_FLOAT, 4 * sizeof(float), (void*)(2 * sizeof(float)));

    // Generates Shader object using the shaders defualt.vert and default.frag
    Shader shaderProgram(".\\Libraries\\shaders\\default.vert", ".\\Libraries\\shaders\\default.frag");
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

    glfwSetKeyCallback(window, Window_Callbacks::define_button_callbacks);

    /*
    
    Create/Open the logging files

    */

    std::ofstream data[4], momentum_data[4];

        open_output_files(data, momentum_data);
  
    /*
    
    Get the metric at the observer to feed into the initial conditions functions

    */

    Initial_conditions_type s_Initial_Conditions{};
    s_Initial_Conditions.init_Pos[e_r] = r_obs;
    s_Initial_Conditions.init_Pos[e_theta] = theta_obs;
    s_Initial_Conditions.init_Pos[e_phi] = phi_obs;

    double metric[4][4]{}, N_obs, omega_obs;

        Spacetimes[e_metric]->get_metric(s_Initial_Conditions.init_metric, &N_obs, &omega_obs, r_obs, theta_obs);
    
        s_Initial_Conditions.init_metric_Redshift_func = N_obs;
        s_Initial_Conditions.init_metric_Shitft_func = omega_obs;

    print_ASCII_art();

    std::cout << "Observer Radial Position [GM/c^2] = " << r_obs << '\n';
    std::cout << "Observer Inclination [deg]        = " << int(theta_obs / M_PI * 180) << '\n';

    if (lens_from_file) {

        /*
        
        Read the initial conditions from file

        */

        double J_data[500]{}, p_theta_data[500]{};
        int Data_number = 0;

        get_geodesic_data(J_data, p_theta_data, &Data_number);

        for (int photon = 0; photon <= Data_number; photon += 1) {

            /*

            This function polulates the initial momentum inside the s_Initial_Conditions struct

            */

            Spacetimes[e_metric]->get_initial_conditions_from_file(&s_Initial_Conditions, J_data, p_theta_data, photon);

            Lens(&s_Initial_Conditions, data, momentum_data);
    
            print_progress(photon, Data_number, lens_from_file);
        }

        std::cout << '\n';
    }
    else{

        /*
        
        Do one scan line in the middle of the image to find the maximum intensity for use in the colormap
        
        */
        
        Normalizing_colormap = true;
  
        for (double H_angle = H_angle_min; H_angle <= H_angle_max - Scan_Step; H_angle += Scan_Step) {

            get_intitial_conditions_from_angles(&s_Initial_Conditions, 0, H_angle);

            Lens(&s_Initial_Conditions, data, momentum_data);

            pixel_count += 1;

        }

        Normalizing_colormap = false;

        /*
        
        Setup a viewing window for the observer and loop trough it

        */
   
        int progress = 0;

        for (double V_angle = V_angle_min; V_angle <= V_angle_max - Scan_Step; V_angle += Scan_Step) {

            print_progress(progress, int((V_angle_max - V_angle_min) / Scan_Step), lens_from_file);

            progress += 1;

            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, pixel_count, pixel_count * V_angle_max / H_angle_max, 0, GL_RGB, GL_FLOAT, texture_buffer);
            // Specify the color of the background
            glClearColor(0.07f, 0.13f, 0.17f, 1.0f);
            // Clean the back buffer and assign the new color to it
            glClear(GL_COLOR_BUFFER_BIT);
            // Draw primitives, number of indices, datatype of indices, index of indices
            glDrawElements(GL_TRIANGLES, sizeof(Vertex_order) / sizeof(float), GL_UNSIGNED_INT, 0);
            // Swap the back buffer with the front buffer
            glfwSwapBuffers(window);
            // Take care of all GLFW events
            
            for (double H_angle = H_angle_max; H_angle >= H_angle_min + Scan_Step; H_angle -= Scan_Step) {

                /*
                
                This function polulates the initial momentum inside the s_Initial_Conditions struct
                
                */

                get_intitial_conditions_from_angles(&s_Initial_Conditions, V_angle, H_angle);

                Lens(&s_Initial_Conditions, data, momentum_data);

                glfwPollEvents();

            }

        }
    
    }            

    auto end_time = std::chrono::high_resolution_clock::now();

    close_output_files(data, momentum_data);

    std::cout << '\n' << "Simulation finished!" << '\n';

    std::cout << "Simulation time: " << std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time);

    while (!glfwWindowShouldClose(window)) {

        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, pixel_count, pixel_count* V_angle_max / H_angle_max, 0, GL_RGB, GL_FLOAT, texture_buffer);
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

}