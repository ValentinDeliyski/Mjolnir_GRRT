/****************************************************************************************************
|                                                                                                   |
|                          ---------  Gravitational Ray Tracer  ---------                           | 
|                                                                                                   |
|    @ Version: 3.5.3                                                                               |
|    @ Author: Valentin Deliyski                                                                    |
|    @ Description: This program numeriaclly integrates the equations of motion                     |
|    for null geodesics and radiative transfer in a curved spacetime,then projects                  |
|    them onto an observer's screen to construct relativistic images of accretion disks             |
|                                                                                                   |
|    @ Supported Spacetimes:                                                                        |
|      * Kerr Black Holes                                                                           |
|      * Static Regular Black Holes                                                                 |
|      * Rotating Traversable Wormholes                                                             |
|      * Janis - Newman - Winicour Naked Singularities                                              |
|                                                                                                   |
|    @ Supported Disk Models                                                                        |
|      * Novikov-Thorne                                                                             |
|      * Generic Optically Thin Disk With Arbitrary Density, Emission and Absorbtion Profiles       |
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

Spacetime_enums e_metric = Naked_Singularity;
Emission_model_enums e_emission = Synchotron_phenomenological;

/* 

Define classes that holds the spacetime properites

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

extern Real r_obs = 1e4;
extern Real theta_obs = 80. / 180 * M_PI;
Real phi_obs = 0;

c_Observer Observer_class(r_obs, theta_obs, phi_obs);

/*

Define the Optically Thin Disk Class

*/

Optically_Thin_Toroidal_Model OTT_Model;

/*

Define the Novikov-Thorne Disk Class

*/

Real r_in = 3;
Real r_out = 1500;

Novikov_Thorne_Model NT_Model(r_in, r_out);

/*

Define some global boolians

*/

extern Const_bool lens_from_file = false;
extern Const_bool truncate = true;

/*

Rendering Engine variables

*/

float Max_Intensity{};
int texture_indexer{};
bool Normalizing_colormap{};
float texture_buffer[TEXTURE_BUFFER_SIZE * 3]{};

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

void print_progress(int current, int max, bool lens_from_file, bool Normalizing_colormap) {

    int current_digits = 1;

    if (current != 0) {

        current_digits = floor(log10f(current) + 1);

    }

    int max_digits = floor(log10f(max) + 1);

    if (current == 0) {

        if (lens_from_file) {

            std::cout << "Number Of Rays Cast: ";

        }
        else if(!Normalizing_colormap){

            std::cout << "Number Of Lines Scanned: ";

        }
        else {

            std::cout << "Progress: ";

        }

        for (int i = 0; i <= max_digits + current_digits; i += 1) {

            std::cout << "0";

        }

    }

    for (int i = 0; i <= max_digits + current_digits + 1; i += 1) {

        std::cout << "\b";

    }

    std::cout << current + 1 << "/" << max + 1 << " ";

}

/*

Setup a viewing window for the observer

*/

double V_angle_min = -atan(30 / r_obs);
double V_angle_max = atan(30 / r_obs);

double H_angle_min = -atan(30 / r_obs);
double H_angle_max = atan(30 / r_obs);

void main() {
    
    auto start_time = std::chrono::high_resolution_clock::now();

    double Scan_Step = (H_angle_max - H_angle_min) / (RESOLUTION - 1);

    GLFWwindow* window = OpenGL_init(H_angle_max / V_angle_max);
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
    s_Initial_Conditions.init_Pos[e_r]     = r_obs;
    s_Initial_Conditions.init_Pos[e_theta] = theta_obs;
    s_Initial_Conditions.init_Pos[e_phi]   = phi_obs;

    double metric[4][4]{}, N_obs, omega_obs;

        Spacetimes[e_metric]->get_metric(s_Initial_Conditions.init_metric, &N_obs, &omega_obs, r_obs, theta_obs);
    
        s_Initial_Conditions.init_metric_Redshift_func = N_obs;
        s_Initial_Conditions.init_metric_Shitft_func   = omega_obs;

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
    
            print_progress(photon, Data_number, lens_from_file, Normalizing_colormap);
        }

        std::cout << '\n';
    }
    else{

        /*
        
        Do one scan line in the middle of the image to find the maximum intensity for use in the colormap
        
        */
        
        Normalizing_colormap = true;

        std::cout << "Initial y = 0 scan to normalize the colormap..." << '\n';
  
        for (int pixel_num = 0; pixel_num <= RESOLUTION - 1; pixel_num++) {

            get_intitial_conditions_from_angles(&s_Initial_Conditions, 0, H_angle_max - pixel_num * Scan_Step);

            Lens(&s_Initial_Conditions, data, momentum_data);

            print_progress(pixel_num, RESOLUTION - 1, lens_from_file, Normalizing_colormap);

        }

        Normalizing_colormap = false;

        /*
        
        Loop trough the viewing window

        */
   
        int progress = 0;

        std::cout << '\n' << "Simulation starts..." << '\n';

        for (int V_pixel_num = 0; V_pixel_num <= RESOLUTION - 1; V_pixel_num++) {

            update_rendering_window(window, V_angle_max / H_angle_max);
            print_progress(progress, RESOLUTION - 1, lens_from_file, Normalizing_colormap);

            progress += 1;

            for (int H_pixel_num = 0; H_pixel_num <= RESOLUTION - 1; H_pixel_num++) {

                /*
                
                This function polulates the initial momentum inside the s_Initial_Conditions struct
                
                */

                get_intitial_conditions_from_angles(&s_Initial_Conditions, V_angle_min + V_pixel_num * Scan_Step,
                                                                           H_angle_max - H_pixel_num * Scan_Step);

                Lens(&s_Initial_Conditions, data, momentum_data);

            }

        }
    
    }            

    auto end_time = std::chrono::high_resolution_clock::now();

    close_output_files(data, momentum_data);

    std::cout << '\n' << "Simulation finished!" << '\n';

    std::cout << "Simulation time: " << std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time);

    while (!glfwWindowShouldClose(window)) {

        update_rendering_window(window, V_angle_max / H_angle_max);

    }

}