/****************************************************************************************************
|                                                                                                   |
|                          ---------  Gravitational Ray Tracer  ---------                           | 
|                                                                                                   |
|    @ Version: 3.9.2                                                                               |
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
#include <iostream>

#include "Disk_Models.h"
#include "General_GR_functions.h"
#include "Lensing.h"

#include "Rendering_Engine.h"
#include "Structs.h"
#include "Sim_Modes.h"
#include "Console_printing.h"

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

Observer_class Observer(r_obs, theta_obs, phi_obs);

/*

Initialize the file manager

*/

File_manager_class File_manager(input_file_path, Truncate_files);

int main() {

    /*

    Define the Optically Thin Disk Class

    */

    Optically_Thin_Toroidal_Model OTT_Model({}, {});

    OTT_Model.precompute_electron_pitch_angles();

    Disk_model_parameters Disk_params{ { HOTSPOT_R_COORD, 0.0, 0.0 },
                                         HOTSPOT_SCALE,
                                         HOTSPOT_REL_SCALE,
                                         N_ELECTRON_EXACT_CGS,
                                         DISK_OPENING_ANGLE,
                                         DISK_CUTOFF_SCALE,
                                         R_Cutoff,
                                         R_0,
                                         DISK_HEIGHT_SCALE,
                                         DISK_RADIAL_SCALE,
                                         T_ELECTRON_EXACT_CGS,
                                         DISK_MAGNETIZATION };

    Emission_law_parameters Emission_params{ EMISSION_SCALE_PHENOMENOLOGICAL,
                                             DISK_ABSORBTION_COEFF,
                                             EMISSION_POWER_LAW,
                                             SOURCE_F_POWER_LAW };

    int result = OTT_Model.load_parameters(&Disk_params, &Emission_params);

    /*

    Define the Novikov-Thorne Disk Class

    */

    Novikov_Thorne_Model NT_Model(r_in, r_out, Spacetimes);


    if (ERROR != result) {

        Console_Printer_class Console_Printer;
        Console_Printer.print_ASCII_art();
        Console_Printer.print_sim_parameters();

        /*

        Get the metric at the observer to feed into the initial conditions struct

        */

        Initial_conditions_type s_Initial_Conditions{};
        s_Initial_Conditions.init_Pos[e_r]     = r_obs;
        s_Initial_Conditions.init_Pos[e_theta] = theta_obs;
        s_Initial_Conditions.init_Pos[e_phi]   = phi_obs;

        Metric_type s_Metric = Spacetimes[e_metric]->get_metric(s_Initial_Conditions.init_Pos);
        
        memcpy(s_Initial_Conditions.init_metric, s_Metric.Metric, sizeof(s_Metric.Metric));
        s_Initial_Conditions.init_metric_Redshift_func = s_Metric.Lapse_function;
        s_Initial_Conditions.init_metric_Shitft_func   = s_Metric.Shift_function;

        memcpy(s_Initial_Conditions.Spacetimes, Spacetimes, sizeof(Spacetimes));

        s_Initial_Conditions.OTT_model = &OTT_Model;
        s_Initial_Conditions.NT_model  = &NT_Model;

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
}