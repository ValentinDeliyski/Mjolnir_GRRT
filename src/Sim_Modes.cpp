#pragma once

#ifndef SIM_MODES

    #define SIM_MODES

    #include "IO_files.h"
    #include "Constants.h"
    #include "General_GR_functions.h"
    #include "Rendering_Engine.h"
    #include "Disk_models.h"
    #include "Lensing.h"
    #include <thread>
    #include "Spacetimes.h"

    #include <iostream>

    void static print_progress(int current, int max, bool lens_from_file) {

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

        for (int i = 0; i <= max_digits + current_digits + 1; i += 1) {

            std::cout << "\b";

        }

        std::cout << current + 1 << "/" << max + 1 << " ";

    }

    void static Rendering_function(Rendering_engine* Renderer) {

        Renderer->OpenGL_init();
        glfwSetKeyCallback(Renderer->window, Rendering_engine::Window_Callbacks::define_button_callbacks);

        while (!glfwWindowShouldClose(Renderer->window)) {

            Renderer->renormalize_colormap();

            Renderer->update_rendering_window();

            std::this_thread::sleep_for(std::chrono::milliseconds(20));

        }

        Renderer->Free_memory();

    }

    void static Generate_Image(Simulation_Context_type* s_Sim_Context, Rendering_engine* Renderer) {

        /*

        Create/Open the logging files

        */

        s_Sim_Context->File_manager->open_image_output_files(int(0));

        /*

        Loop trough the viewing window

        */

        auto start_time = std::chrono::high_resolution_clock::now();
        int progress = 0;

        std::cout << '\n' << "Generating image..." << '\n';

        for (int V_pixel_num = 0; V_pixel_num <= RESOLUTION - 1; V_pixel_num++) {

            print_progress(progress, RESOLUTION - 1, false);

            progress += 1;

            for (int H_pixel_num = 0; H_pixel_num <= RESOLUTION * H_angle_max / V_angle_max - 1; H_pixel_num++) {

                /*

                This function polulates the initial momentum inside the s_Initial_Conditions struct

                */

                get_intitial_conditions_from_angles(s_Sim_Context->p_Init_Conditions,
                                                    V_angle_min + V_pixel_num * Scan_Step,
                                                    H_angle_max - H_pixel_num * Scan_Step);

               /*

               Ray propagation happens here

               */

               Results_type* s_Ray_results = Propagate_ray(s_Sim_Context);

               /*

               Updating the visualization happens here

               */

               Renderer->Intensity_buffer[int(Renderer->texture_indexer / 3)] = s_Ray_results->Intensity[direct][I] +
                                                                              s_Ray_results->Intensity[first][I] +
                                                                              s_Ray_results->Intensity[second][I] +
                                                                              s_Ray_results->Intensity[third][I];

               Renderer->texture_indexer += 3;

               s_Sim_Context->File_manager->write_image_data_to_file(s_Ray_results);

            }

        }

        auto end_time = std::chrono::high_resolution_clock::now();

        std::cout << '\n' << "Image Generation Finished!";
        std::cout << '\n' << "Simulation time: " << std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time) << "\n";

        s_Sim_Context->File_manager->close_image_output_files();

    }

    void run_simulation_mode_1(Simulation_Context_type* p_Sim_Context) {


        /*

        Initialize the rendering engine (the Renderer instance must be static to not blow up the stack - the texture and intensity buffer are inside of it)

        */

        static Rendering_engine Renderer = Rendering_engine();

        std::jthread GUI_Thread(Rendering_function, &Renderer);

        /*
        
        Run the main ray-tracer "kernel"
        
        */

        Generate_Image(p_Sim_Context, &Renderer);


    }

    void run_simulation_mode_2(Simulation_Context_type* p_Sim_Context) {

        /*

        Read the initial conditions from file

        */

        double J_data[500]{}, p_theta_data[500]{};
        int Data_number{};

        p_Sim_Context->File_manager->get_geodesic_data(J_data, p_theta_data, &Data_number);
        /*

        Create/Open the logging files

        */

        p_Sim_Context->File_manager->open_image_output_files(Data_number);

        for (int Param_Sweep = 0; Param_Sweep <= PARAM_SWEEP_NUMBER - 1; Param_Sweep++) {

            double Current_Param_Value = INIT_PARAM_VALUE;

            if (PARAM_SWEEP_NUMBER > 1) {

                Current_Param_Value = INIT_PARAM_VALUE * (1.0f - double(Param_Sweep) / (PARAM_SWEEP_NUMBER - 1)) + FINAL_PARAM_VALUE * Param_Sweep / (PARAM_SWEEP_NUMBER - 1);
            }

            p_Sim_Context->p_Spacetime->update_parameters(Current_Param_Value, PARAM_TYPE);

            for (int photon = 0; photon <= Data_number - 1; photon += 1) {

                /*

                This function polulates the initial momentum inside the s_Initial_Conditions struct

                */

                p_Sim_Context->p_Spacetime->get_initial_conditions_from_file(p_Sim_Context->p_Init_Conditions, J_data, p_theta_data, photon);

                /*

                Ray propagation happens here

                */

                Results_type* s_Ray_results = Propagate_ray(p_Sim_Context);

                p_Sim_Context->File_manager->write_image_data_to_file(s_Ray_results);

                print_progress(photon, Data_number - 1, true);
            }

            std::cout << '\n';
        }

        p_Sim_Context->File_manager->close_image_output_files();
    }

    void run_simulation_mode_4(Simulation_Context_type* p_Sim_Context) {

        p_Sim_Context->p_Spacetime->get_initial_conditions_from_file(p_Sim_Context->p_Init_Conditions, (double*) &X_INIT, (double*) &Y_INIT, 0);

        Results_type* s_Ray_results = Propagate_ray(p_Sim_Context);

        p_Sim_Context->File_manager->open_log_output_file();
        p_Sim_Context->File_manager->log_photon_path(s_Ray_results);
        p_Sim_Context->File_manager->close_log_output_file();
    }

#endif

