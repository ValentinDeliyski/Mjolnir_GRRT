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
        p_Sim_context->p_Spacetime = std::make_shared<Kerr_class>(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case Wormhole:      
        p_Sim_context->p_Spacetime = std::make_shared<Wormhole_class>(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case Reg_Black_Hole:       
        p_Sim_context->p_Spacetime = std::make_shared<RBH_class>(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case Janis_Newman_Winicour:       
        p_Sim_context->p_Spacetime = std::make_shared<JNW_class>(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case Einstein_Gauss_Bonnet:       
        p_Sim_context->p_Spacetime = std::make_shared<Gauss_Bonnet_class>(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case BH_w_Dark_Matter:      
        p_Sim_context->p_Spacetime = std::make_shared<Black_Hole_w_Dark_Matter_Halo_class>(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case Numerical:
        p_Sim_context->p_Spacetime = std::make_shared<Numerical_metric>(&p_Sim_context->p_Init_Conditions->Metric_parameters);
        break;

    case Minkowski:
        p_Sim_context->p_Spacetime = std::make_shared<Minkowski_class>(&p_Sim_context->p_Init_Conditions->Metric_parameters);
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

        s_Sim_Context.p_Init_Conditions = std::make_shared<Initial_conditions_type>(Initial_conditions_type());

        if (Return_Values::ERROR == parse_simulation_input_XML(Input_file_path, s_Sim_Context.p_Init_Conditions)) { throw std::runtime_error("Could not parse input file!"); }

        s_Sim_Context.p_Init_Conditions->Print_to_console = print_to_console;
        s_Sim_Context.p_Init_Conditions->Hotspot_params.Profile_params.Coord_time_offset += s_Sim_Context.p_Init_Conditions->Observer_params.distance;

        // Populate the Spacetime class instance 
        Allocate_Spacetime_Class(&s_Sim_Context);

        // Get the observer position and populate the Observer class instance.
        s_Sim_Context.p_Observer = std::make_shared<Observer_class>(&s_Sim_Context);

        double init_state[4] = { s_Sim_Context.p_Init_Conditions->Observer_params.init_time,
                                 s_Sim_Context.p_Init_Conditions->Observer_params.distance,
                                 s_Sim_Context.p_Init_Conditions->Observer_params.inclination,
                                 s_Sim_Context.p_Init_Conditions->Observer_params.azimuth };

        Metric_type s_init_Metric = s_Sim_Context.p_Spacetime->get_global_metric(init_state);
        memcpy(&s_Sim_Context.p_Init_Conditions->Init_metric, &s_init_Metric, sizeof(Metric_type));

        // Populate the Emission Model class instances
        s_Sim_Context.p_Emission_Model = std::make_shared<Emission_models_class>(&s_Sim_Context);
        s_Sim_Context.p_Emission_Model->precompute_electron_pitch_angles(s_Sim_Context.p_Init_Conditions);

        if (Disk_model_enums::e_Novikov_Thorne == s_Sim_Context.p_Init_Conditions->Disk_params.e_Disk_model) {

            /* ---- Allocate the Novikov - Thorne Model class
                    NOTE: This is a unique pointer because it does not get copied anywhere ---- */
            s_Sim_Context.p_NT_model = std::make_unique<Novikov_Thorne_Model_class>(&s_Sim_Context);
            
        }
        else {

            s_Sim_Context.p_NT_model = nullptr;

        }

        /* --- Initialize the struct that holds the ray results --- */
        std::unique_ptr<Results_type> s_Ray_results = std::make_unique<Results_type>();

        /* --- The ray log struct will need to be accessed by the geodesic and emission integrators, so I create it as a shared pointer --- */
        s_Ray_results->Ray_log_struct = std::make_shared<Ray_log_type>();

        /* -------------------------- Reference for the sake of readability -------------------------- */
        size_t& Max_log_size = s_Sim_Context.p_Init_Conditions->Integrator_params.Max_integration_count;

        /* --- The logs themselves should only exist in one point in memory, so I create them as unique pointers --- */
        s_Ray_results->Ray_log_struct->Ray_path_log_local = std::make_unique<double[]>(Max_log_size * e_Full_state_size);
        s_Ray_results->Ray_log_struct->Ray_path_log_global = std::make_unique<double[]>(Max_log_size * e_Full_state_size);

        s_Ray_results->RK_integrator_debug_log.N_steps_rejected = std::make_unique<double[]>(Max_log_size * e_Full_state_size);;
        s_Ray_results->RK_integrator_debug_log.State_error_history = std::make_unique<double[]>(Max_log_size * e_Full_state_size);;

        for (int index = I; index < e_Stokes_param_num; index++) {

            s_Ray_results->Ray_log_struct->Ray_emission_log[index] = std::make_unique<double[]>(Max_log_size);

        }

        for (int index = e_x; index <= e_y; index++) {

            s_Ray_results->Ray_log_struct->Ray_polarization_log[index] = std::make_unique<double[]>(Max_log_size);
            s_Ray_results->Polarization_debug_log.PW_constant[index] = new double[Max_log_size] {};

        }

        if (s_Sim_Context.p_Init_Conditions->Print_to_console) {

            print_ASCII_art();

        }

        /* ============================== Run the simulation ============================== */

        switch (s_Sim_Context.p_Init_Conditions->Simulation_mode) {

        case Image_generation:
            run_image_generation(&s_Sim_Context, s_Ray_results.get());
            break;

        case Make_geodesic_sweep:
            run_geodesic_sweep(&s_Sim_Context, s_Ray_results.get());
            break;

        case Make_geodesic_log:
            make_geodesic_log(&s_Sim_Context, s_Ray_results.get());
            break;

        case Debug_mode:
            run_debug_simulation(&s_Sim_Context);
            break;

        default:

            throw std::runtime_error("Unsupported simulation mode! \n");

        }

        return Return_Values::OK;

    }
    catch (const std::exception& error) {

        std::cout << std::endl << "Mjolnir ERROR: " << error.what() << std::endl << std::endl;
        return Return_Values::ERROR;

    }

}