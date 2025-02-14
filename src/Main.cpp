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

#include "Input_parser.h"


void static Allocate_Spacetime_Class(Simulation_Context_type* p_Sim_context) {

    // These do not ever get "delete" called on them, because they need to exist for the entire duration of the program

    switch (p_Sim_context->p_Init_Conditions->Metric_parameters.e_Spacetime) {

    case Kerr:
        p_Sim_context->p_Spacetime = new Kerr_class;
        break;

    case Wormhole:      
        p_Sim_context->p_Spacetime = new Wormhole_class;
        break;

    case Reg_Black_Hole:       
        p_Sim_context->p_Spacetime = new RBH_class;
        break;

    case Janis_Newman_Winicour:       
        p_Sim_context->p_Spacetime = new JNW_class;
        break;

    case Einstein_Gauss_Bonnet:       
        p_Sim_context->p_Spacetime = new Gauss_Bonnet_class;
        break;

    case BH_w_Dark_Matter:      
        p_Sim_context->p_Spacetime = new Black_Hole_w_Dark_Matter_Halo_class;
        break;

    }

}

void static Allocate_GOT_Model_class_instance(Simulation_Context_type* p_Sim_Context) {

    p_Sim_Context->p_GOT_Model = new Generic_Optically_Thin_Model();
    p_Sim_Context->p_GOT_Model->precompute_electron_pitch_angles(p_Sim_Context->p_Init_Conditions);

    if (ERROR == p_Sim_Context->p_GOT_Model->load_parameters(p_Sim_Context)) {

        exit(ERROR);

    }
}

int main(int argument_count, char** cmd_line_args) {

    std::string Input_file_path{};
    bool print_to_console{};
    if (argument_count == 5 && 0 == strcmp(cmd_line_args[1], "-in") && 0 == strcmp(cmd_line_args[3], "-print_to_console")) {

        Input_file_path = cmd_line_args[2];
        print_to_console = std::stoi(cmd_line_args[4]);

    }
    else {

        std::cout << "To run Mjolnir, use the following call structure:" << "\n";
        std::cout << "Mjolnir_GRRT.exe -in __INPUT_FILE_PATH__ -print_to_console __1 FOR YES 0 FOR NO__" << "\n";

        exit(ERROR);

    }

    /*
    
    |============================== Define the Simulation Context struct ==============================|
    
    */

    Simulation_Context_type s_Sim_Context{};

    s_Sim_Context.p_Init_Conditions = new Initial_conditions_type();

    if (ERROR == parse_simulation_input_XML(Input_file_path, s_Sim_Context.p_Init_Conditions)){
    
        exit(ERROR);
    
    }

    s_Sim_Context.p_Init_Conditions->Print_to_console = print_to_console;
    s_Sim_Context.p_Init_Conditions->Hotspot_params.Coord_time_offset += s_Sim_Context.p_Init_Conditions->Observer_params.distance;

    // Populate the Spacetime class instance 
    Allocate_Spacetime_Class(&s_Sim_Context);

    s_Sim_Context.p_Spacetime->load_parameters(s_Sim_Context.p_Init_Conditions->Metric_parameters);

    // Get the observer position and populate the Observer class instance.
    s_Sim_Context.p_Observer = new Observer_class(&s_Sim_Context);

    double init_state[4] = {0,
                            s_Sim_Context.p_Init_Conditions->Observer_params.distance,
                            s_Sim_Context.p_Init_Conditions->Observer_params.inclination,
                            s_Sim_Context.p_Init_Conditions->Observer_params.azimuth };

    Metric_type s_init_Metric = s_Sim_Context.p_Spacetime->get_metric(init_state);
    
    memcpy(&s_Sim_Context.p_Init_Conditions->Init_metric, &s_init_Metric, sizeof(Metric_type));

    // Populate the Emission Model class instances
    Allocate_GOT_Model_class_instance(&s_Sim_Context);

    // Allocate the Novikov-Thorne Model class
     s_Sim_Context.p_NT_model = new Novikov_Thorne_Model(&s_Sim_Context);

    // Populate the File Manager class instance
    s_Sim_Context.File_manager = new File_manager_class(s_Sim_Context.p_Init_Conditions);

    // Initialize the struct that holds the ray results (as static in order to not blow up the stack -> this must always be passed around as a pointer!)
    static Results_type s_Ray_results{};

    s_Ray_results.Ray_log_struct.Ray_path_log = new double[s_Sim_Context.p_Init_Conditions->Integrator_params.Max_integration_count * e_State_Number]();

    for (int index = I; index <= STOKES_PARAM_NUM - 1; index++) {

        s_Ray_results.Ray_log_struct.Ray_emission_log[index] = new double[2 * s_Sim_Context.p_Init_Conditions->Integrator_params.Max_integration_count]();

    }

    Console_Printer_class Console_Printer;
    
    if (s_Sim_Context.p_Init_Conditions->Print_to_console) { 

        Console_Printer.print_ASCII_art();
        Console_Printer.print_sim_parameters(s_Sim_Context.p_Init_Conditions);

    }

    /*

    |============================== Run the simulation ==============================|

    */

    switch (s_Sim_Context.p_Init_Conditions->Simulation_mode) {
   
    case 1:
         run_simulation_mode_1(&s_Sim_Context, &s_Ray_results);
         break;
   
    case 2:
         run_simulation_mode_2(&s_Sim_Context, &s_Ray_results);
         break;

    case 3:
         run_simulation_mode_3(&s_Sim_Context, &s_Ray_results);
         break;
   
    default:
         std::cout << "Unsuported simulation mode!" << "\n";
         exit(ERROR);
   
    }

    return OK;

}