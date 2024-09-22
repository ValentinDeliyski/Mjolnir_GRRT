/****************************************************************************************************
|                                                                                                   |
|                ---------  Mjølnir General Relativistic Ray Tracer  ---------                      | 
|                                                                                                   |
|    @ Version: 2.0                                                                                 |
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


void static Allocate_Spacetime_Class(Spacetime_enums e_metric, Spacetime_Base_Class** p_Spacetime) {

    // These do not ever get "delete" called on them, because they need to exist for the entire duration of the program

    switch (e_metric) {

    case Kerr:
        *p_Spacetime = new Kerr_class;
        break;

    case Wormhole:      
        *p_Spacetime = new Wormhole_class;
        break;

    case Reg_Black_Hole:       
        *p_Spacetime = new RBH_class;
        break;

    case Naked_Singularity:       
        *p_Spacetime = new JNW_class;
        break;

    case Gauss_Bonnet:       
        *p_Spacetime = new Gauss_Bonnet_class;
        break;

    case BH_w_Dark_Matter:      
        *p_Spacetime = new Black_Hole_w_Dark_Matter_Halo_class;
        break;

    }

}

void static Allocate_GOT_Model_class_instance(Simulation_Context_type* s_Sim_Context) {

    static Generic_Optically_Thin_Model GOT_Model_class_instance = Generic_Optically_Thin_Model();
    s_Sim_Context->p_GOT_Model = &(GOT_Model_class_instance);

    s_Sim_Context->p_GOT_Model->precompute_electron_pitch_angles();

    // TODO: These need to get loaded from a file...

    s_Sim_Context->p_Init_Conditions->Disk_params = { { HOTSPOT_R_COORD, M_PI_2, 0.0 },
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

    s_Sim_Context->p_Init_Conditions->Emission_params = { EMISSION_SCALE_PHENOMENOLOGICAL,
                                                          DISK_ABSORBTION_COEFF,
                                                          EMISSION_POWER_LAW,
                                                          SOURCE_F_POWER_LAW };

    if (ERROR == s_Sim_Context->p_GOT_Model->load_parameters(&s_Sim_Context->p_Init_Conditions->Disk_params, &s_Sim_Context->p_Init_Conditions->Emission_params)) {

        std::cout << "Could not load emission model parameters!" << "\n";
        exit(ERROR);
    }
}

int main() {
    
    Console_Printer_class Console_Printer;
    Console_Printer.print_ASCII_art();
    Console_Printer.print_sim_parameters();

    /*
    
    |============================== Define the Simulation Context struct ==============================|
    
    */

    Simulation_Context_type s_Sim_Context{};

    // Populate the Spacetime class instance 
    Allocate_Spacetime_Class(e_metric, &s_Sim_Context.p_Spacetime);

    // TODO: Read this from a file...

    Initial_conditions_type Init_Conditions{};
    s_Sim_Context.p_Init_Conditions = &Init_Conditions;

    s_Sim_Context.p_Init_Conditions->Metric_Parameters = { WH_REDSHIFT,
                                                            STOP_AT_THROAT,
                                                            JNW_GAMMA,
                                                            GAUSS_BONNET_GAMMA,
                                                            RBH_PARAM,
                                                            COMPACTNESS,
                                                            M_HALO,
                                                            SPIN };

    s_Sim_Context.p_Spacetime->load_parameters(s_Sim_Context.p_Init_Conditions->Metric_Parameters);
    s_Sim_Context.e_Spacetime = e_metric;

    // Populate the Initial Conditions
    s_Sim_Context.p_Init_Conditions->init_Pos[e_r]     = r_obs;
    s_Sim_Context.p_Init_Conditions->init_Pos[e_theta] = theta_obs;
    s_Sim_Context.p_Init_Conditions->init_Pos[e_phi]   = phi_obs;

    Metric_type s_init_Metric = s_Sim_Context.p_Spacetime->get_metric(s_Sim_Context.p_Init_Conditions->init_Pos);
    
    memcpy(s_Sim_Context.p_Init_Conditions->init_metric, s_init_Metric.Metric, sizeof(s_init_Metric.Metric));
    s_Sim_Context.p_Init_Conditions->init_metric_Redshift_func = s_init_Metric.Lapse_function;
    s_Sim_Context.p_Init_Conditions->init_metric_Shitft_func   = s_init_Metric.Shift_function;

    // Populate the Emission Model class instances
    Allocate_GOT_Model_class_instance(&s_Sim_Context);

    Novikov_Thorne_Model NT_class_instance(r_in, r_out, s_Sim_Context.p_Spacetime);
    s_Sim_Context.p_NT_model = &NT_class_instance;

    // Populate the Observer class instance
    Observer_class Observer_class_instance(s_Sim_Context.p_Init_Conditions);
    s_Sim_Context.p_Observer = &Observer_class_instance;

    // Populate the File Manager class instance
    File_manager_class File_manager_instance(s_Sim_Context.p_Init_Conditions, input_file_path, Truncate_files);
    s_Sim_Context.File_manager = &File_manager_instance;

    /*

    |============================== Run the simulation ==============================|

    */

    switch (Active_Sim_Mode) {
   
    case 1:
         run_simulation_mode_1(&s_Sim_Context);
         break;
   
    case 2:
         run_simulation_mode_2(&s_Sim_Context);
         break;
   
    case 3:
   
         // With the python wrapper this wont need to exist anymore (yay)

         //run_simulation_mode_3(s_Sim_Context);
   
         break;
   
    case 4:
         run_simulation_mode_4(&s_Sim_Context);
         break;
   
    default:
         std::cout << "Unsuported simulation mode!" << "\n";
         exit(ERROR);
   
    }

    return OK;

}