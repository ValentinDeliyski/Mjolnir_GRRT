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
#include <string>
#include <vector>

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

std::vector<c_Spacetime_Base*> Spacetimes = {

    new derived_Kerr_class(),
    new derived_RBH_class(),
    new derived_Wormhole_class(),
    new derived_JNW_class()
};

/*

Define the Observer class

*/

c_Observer Observer_class(r_obs, theta_obs, phi_obs);

/*

Define the Optically Thin Disk Class

*/

Optically_Thin_Toroidal_Model OTT_Model;

/*

Define the Novikov-Thorne Disk Class

*/

Novikov_Thorne_Model NT_Model(r_in, r_out);

/*

Define some global boolians

*/

extern Const_bool truncate = true;

/*

Rendering Engine variables

*/

float Max_Intensity{};
int texture_indexer{};
bool Normalizing_colormap = true;
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

int main() {

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

    switch (Active_Sim_Mode) {

        case 1:

            run_simulation_mode_1(&s_Initial_Conditions, data, momentum_data);

            break;

        case 2:

            run_simulation_mode_2(&s_Initial_Conditions, data, momentum_data);

            break;

        default:

            std::cout << "Unsuported simulation mode!";

            return ERROR;

    }

    close_output_files(data, momentum_data);

    return OK;
}