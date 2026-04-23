#define _USE_MATH_DEFINES

#include "Constants.h"
#include "Spacetimes.h"
#include "Enumerations.h"
#include "IO_files.h"
#include <iostream>

#include "Emission_Models.h"
#include "Novikov_Thorne_model.h"
#include "General_GR_functions.h"
#include "Lensing.h"

#include "Rendering_Engine.h"
#include "Structs.h"
#include "Sim_Modes.h"

#include "Input_parser.h"

void static print_ASCII_art() {

  std::cout << " __       __                    __            __                   ______   _______   _______   ________   \n"
            << "/  \\     /  |                  /  |          /  |                 /      \\ /       \\ /      \\ /        |  \n"
            << "$$  \\   /$$ |     __   ______  $$ | _______  $$/   ______        /$$$$$$  |$$$$$$$  |$$$$$$$  |$$$$$$$$/   \n"
            << "$$$  \\ /$$$ |    /  | /      \\ $$ |/       \\ /  | /      \\       $$ | _$$/ $$ |__$$ |$$ |__$$ |   $$ |     \n"
            << "$$$$  /$$$$ |    $$/ /$$$$$$  |$$ |$$$$$$$  |$$ |/$$$$$$  |      $$ |/    |$$    $$< $$    $$<    $$ |     \n"
            << "$$ $$ $$/$$ |    /  |$$ |  $$ |$$ |$$ |  $$ |$$ |$$ |  $$/       $$ |$$$$ |$$$$$$$  |$$$$$$$  |   $$ |     \n"
            << "$$ |$$$/ $$ |    $$ |$$ \\__$$ |$$ |$$ |  $$ |$$ |$$ |            $$ \\__$$ |$$ |  $$ |$$ |  $$ |   $$ |     \n"
            << "$$ | $/  $$ |    $$ |$$    $$/ $$ |$$ |  $$ |$$ |$$ |            $$    $$/ $$ |  $$ |$$ |  $$ |   $$ |     \n"
            << "$$/      $$/__   $$ | $$$$$$/  $$/ $$/   $$/ $$/ $$/              $$$$$$/  $$/   $$/ $$/   $$/    $$/      \n"
            << "           /  \\__$$ |                                                                                      \n"
            << "           $$    $$/                                                                                       \n"
            << "            $$$$$$/                                                                                        \n";

        std::cout << '\n';

}

void static Allocate_Spacetime_Class(Simulation_Context_type* p_Sim_context) {

    // These do not ever get "delete" called on them, because they need to exist for the entire duration of the program

    switch (p_Sim_context->p_Init_Conditions->Metric_parameters.e_Spacetime) {

    case Kerr:
        p_Sim_context->p_Spacetime = new Kerr_class(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case Wormhole:      
        p_Sim_context->p_Spacetime = new Wormhole_class(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case Reg_Black_Hole:       
        p_Sim_context->p_Spacetime = new RBH_class(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case Janis_Newman_Winicour:       
        p_Sim_context->p_Spacetime = new JNW_class(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case Einstein_Gauss_Bonnet:       
        p_Sim_context->p_Spacetime = new Gauss_Bonnet_class(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case BH_w_Dark_Matter:      
        p_Sim_context->p_Spacetime = new Black_Hole_w_Dark_Matter_Halo_class(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case Numerical:
        p_Sim_context->p_Spacetime = new Numerical_metric(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case Minkowski:
        p_Sim_context->p_Spacetime = new Minkowski_class(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    }

}

int main(int argument_count, char** cmd_line_args) {

    try {

        std::string Input_file_path{};
        bool print_to_console{};

        if (argument_count == 5 and 0 == strcmp(cmd_line_args[1], "-in") and 0 == strcmp(cmd_line_args[3], "-print_to_console")) {

            Input_file_path = cmd_line_args[2];
            print_to_console = std::stoi(cmd_line_args[4]);

        }
        else {

            throw std::runtime_error("To run Mjolnir, use the following call structure:\n Mjolnir_GRRT.exe -in __INPUT_FILE_PATH__ -print_to_console __1 FOR YES 0 FOR NO__ \n");

        }

        /*

        |============================== Define the Simulation Context struct ==============================|

        */

        Simulation_Context_type s_Sim_Context{};

        s_Sim_Context.p_Init_Conditions = new Initial_conditions_type();

        if (ERROR == parse_simulation_input_XML(Input_file_path, s_Sim_Context.p_Init_Conditions)) { throw std::runtime_error("Could not parse input file!"); }

        s_Sim_Context.p_Init_Conditions->Print_to_console = print_to_console;
        s_Sim_Context.p_Init_Conditions->Hotspot_params.Profile_params.Coord_time_offset += s_Sim_Context.p_Init_Conditions->Observer_params.distance;

        // Populate the Spacetime class instance 
        Allocate_Spacetime_Class(&s_Sim_Context);

        // Get the observer position and populate the Observer class instance.
        s_Sim_Context.p_Observer = new Observer_class(&s_Sim_Context);

        double init_state[4] = { s_Sim_Context.p_Init_Conditions->Observer_params.init_time,
                                 s_Sim_Context.p_Init_Conditions->Observer_params.distance,
                                 s_Sim_Context.p_Init_Conditions->Observer_params.inclination,
                                 s_Sim_Context.p_Init_Conditions->Observer_params.azimuth };

        Metric_type s_init_Metric = s_Sim_Context.p_Spacetime->get_global_metric(init_state);

        memcpy(&s_Sim_Context.p_Init_Conditions->Init_metric, &s_init_Metric, sizeof(Metric_type));

        // Populate the Emission Model class instances
        s_Sim_Context.p_Emission_Model = new Emission_models_class(&s_Sim_Context);
        s_Sim_Context.p_Emission_Model->precompute_electron_pitch_angles(s_Sim_Context.p_Init_Conditions);


        if (s_Sim_Context.p_Init_Conditions->Disk_params.e_Disk_model == e_Novikov_Thorne) {

            // Allocate the Novikov-Thorne Model class
            s_Sim_Context.p_NT_model = new Novikov_Thorne_Model_class(&s_Sim_Context);
            
        }
        else {

            s_Sim_Context.p_NT_model = nullptr;

        }

        // Populate the File Manager class instance
        s_Sim_Context.File_manager = new File_manager_class(s_Sim_Context.p_Init_Conditions);

        // Initialize the struct that holds the ray results (as static in order to not blow up the stack -> this must always be passed around as a pointer!)
        static Results_type s_Ray_results{};

        s_Ray_results.Ray_log_struct.Ray_path_log_local = new double[s_Sim_Context.p_Init_Conditions->Integrator_params.Max_integration_count * e_Full_state_size] {};
        s_Ray_results.Ray_log_struct.Ray_path_log_global = new double[s_Sim_Context.p_Init_Conditions->Integrator_params.Max_integration_count * e_Full_state_size] {};

        s_Ray_results.RK_integrator_debug_log.N_steps_rejected = new double[s_Sim_Context.p_Init_Conditions->Integrator_params.Max_integration_count] {};
        s_Ray_results.RK_integrator_debug_log.State_error_history = new double[s_Sim_Context.p_Init_Conditions->Integrator_params.Max_integration_count] {};

        s_Ray_results.Polarization_debug_log.PW_constant[0] = new double[s_Sim_Context.p_Init_Conditions->Integrator_params.Max_integration_count] {};
        s_Ray_results.Polarization_debug_log.PW_constant[1] = new double[s_Sim_Context.p_Init_Conditions->Integrator_params.Max_integration_count] {};

        for (int index = I; index < e_Stokes_param_num; index++) {

            s_Ray_results.Ray_log_struct.Ray_emission_log[index] = new double[s_Sim_Context.p_Init_Conditions->Integrator_params.Max_integration_count]();

        }

        for (int index = e_x; index <= e_y; index++) {

            s_Ray_results.Ray_log_struct.Ray_polarization_log[index] = new double[s_Sim_Context.p_Init_Conditions->Integrator_params.Max_integration_count]();

        }

        if (s_Sim_Context.p_Init_Conditions->Print_to_console) {

            print_ASCII_art();

        }

        /*

        |============================== Run the simulation ==============================|

        */

        switch (s_Sim_Context.p_Init_Conditions->Simulation_mode) {

        case Image_generation:
            run_image_generation(&s_Sim_Context, &s_Ray_results);
            break;

        case Make_geodesic_sweep:
            run_geodesic_sweep(&s_Sim_Context, &s_Ray_results);
            break;

        case Make_geodesic_log:
            make_geodesic_log(&s_Sim_Context, &s_Ray_results);
            break;

        case Debug_mode:
            run_debug_simulation(&s_Sim_Context);
            break;

        default:

            throw std::runtime_error("Unsupported simulation mode! \n");

        }

        return OK;

    }
    catch (const std::exception& error) {

        std::cout << std::endl << "Mjolnir ERROR: " << error.what() << std::endl << std::endl;
        return ERROR;

    }

}