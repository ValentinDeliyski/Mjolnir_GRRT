/****************************************************************************************************
|                                                                                                   |
|                          ---------  Gravitational Ray Tracer  ---------                           | 
|                                                                                                   |
|    @ Version: 3.8.1                                                                               |
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
|      * Gauss - Bonnet Black Holes / Naked Singularities                                           |
|      * Black Holes with a Dark Matter Halo                                                        |
|                                                                                                   |
|    @ Supported Disk Models                                                                        |
|      * Novikov-Thorne                                                                             |
|      * Generic Optically Thin Disk With Arbitrary Density, Emission and Absorbtion Profiles       |
|                                                                                                   |
****************************************************************************************************/

#define _USE_MATH_DEFINES

#include "Constants.h"
#include "Spacetimes.h"
#include "Enumerations.h"
#include "IO_files.h"

#include "Disk_Models.h"
#include "General_GR_functions.h"
#include "Lensing.h"

#include "Rendering_Engine.h"

#include "Sim_Modes.h"

/* 

Define classes that holds the spacetime properites

*/


Kerr_class Kerr_class_instance = Kerr_class();
Wormhole_class Wormhole_class_instance = Wormhole_class();
RBH_class RHB_class_instance = RBH_class();
JNW_class JNW_class_intance = JNW_class();
Gauss_Bonnet_class Gauss_Bonet_class_instance = Gauss_Bonnet_class();
Black_Hole_w_Dark_Matter_Halo_class BH_w_DM_class_instace = Black_Hole_w_Dark_Matter_Halo_class();

Spacetime_Base_Class* Spacetimes[] = {

    &Kerr_class_instance,
    &Wormhole_class_instance,
    &RHB_class_instance,
    &JNW_class_intance,
    &Gauss_Bonet_class_instance,
    &BH_w_DM_class_instace

};

/*

Define the Observer class

*/

c_Observer Observer_class(r_obs, theta_obs, phi_obs);

/*

Define the Optically Thin Disk Class

*/

Optically_Thin_Toroidal_Model OTT_Model(HOTSPOT_R_COORD, 0);

/*

Define the Novikov-Thorne Disk Class

*/

Novikov_Thorne_Model NT_Model(r_in, r_out);

/*

Rendering Engine variables

*/

float Max_Intensity{};
float texture_buffer[TEXTURE_BUFFER_SIZE * 3]{};

/*

Initialize the file manager

*/

File_manager_class File_manager(input_file_path, Truncate_files);

/*

Precomputed variables definitions

*/

double sin_electron_pitch_angles[NUM_SAMPLES_TO_AVG]{};
double one_over_sqrt_sin[NUM_SAMPLES_TO_AVG];
double one_over_cbrt_sin[NUM_SAMPLES_TO_AVG];
double one_over_sin_to_1_6[NUM_SAMPLES_TO_AVG];

void precompute_electron_pitch_angles() {

    for (int index = 0; index <= NUM_SAMPLES_TO_AVG - 1; index++) {

        double pitch_angle = double(index) / NUM_SAMPLES_TO_AVG * M_PI;
        sin_electron_pitch_angles[index] = sin(pitch_angle);

        if (sin_electron_pitch_angles[index] != 0) {

            one_over_sqrt_sin[index]   = 1. / sqrt(sin_electron_pitch_angles[index]);
            one_over_cbrt_sin[index]   = 1. / cbrt(sin_electron_pitch_angles[index]);
            one_over_sin_to_1_6[index] = 1. / sqrt(one_over_cbrt_sin[index]);

        }

    }

}

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

int main() {

    precompute_electron_pitch_angles();

    /*

    Get the metric at the observer to feed into the initial conditions struct

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

    switch (Active_Sim_Mode) {

        case 1:

            run_simulation_mode_1(&s_Initial_Conditions);

            break;

        case 2:

            run_simulation_mode_2(&s_Initial_Conditions);

            break;

        case 3:

            run_simulation_mode_3(&s_Initial_Conditions);

            break;

        case 4:

            run_simulation_mode_4(&s_Initial_Conditions);

            break;

        default:

            std::cout << "Unsuported simulation mode!";

            return ERROR;

    }

    return OK;
}