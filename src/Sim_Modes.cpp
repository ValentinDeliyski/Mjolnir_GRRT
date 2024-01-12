#pragma once

#ifndef SIM_MODES

    #define SIM_MODES

    #include "IO_files.h"
    #include "Constants.h"
    #include "General_GR_functions.h"
    #include "Rendering_Engine.h"
    #include "Disk_models.h"
    #include "Lensing.h"
    #include "Spacetimes.h"

    #include <iostream>

    extern File_manager_class File_manager;

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

        for (int i = 0; i <= max_digits + current_digits + 1; i += 1) {

            std::cout << "\b";

        }

        std::cout << current + 1 << "/" << max + 1 << " ";

    }

    void Rendering_function(Rendering_engine* Renderer) {

        Renderer->OpenGL_init();
        glfwSetKeyCallback(Renderer->window, Rendering_engine::Window_Callbacks::define_button_callbacks);

        while (!glfwWindowShouldClose(Renderer->window)) {

            Renderer->renormalize_colormap();

            Renderer->update_rendering_window();

            std::this_thread::sleep_for(std::chrono::milliseconds(20));

        }

    }

    void run_simulation_mode_1(Initial_conditions_type* s_Initial_Conditions) {

        /*

        Create/Open the logging files

        */

        File_manager.open_image_output_files();

        /*

        Initialize the rendering engine (the Renderer instance must be static to not blow up the stack - the texture and intensity buffer are inside of it)

        */

        static Rendering_engine Renderer = Rendering_engine();

        std::jthread GUI_Thread(Rendering_function, &Renderer);

        /*

        Loop trough the viewing window

        */

        auto start_time = std::chrono::high_resolution_clock::now();
        int progress = 0;

        std::cout << '\n' << "Simulation starts..." << '\n';

        for (int V_pixel_num = 0; V_pixel_num <= RESOLUTION - 1; V_pixel_num++) {

            print_progress(progress, RESOLUTION - 1, false);

            progress += 1;

            for (int H_pixel_num = 0; H_pixel_num <= RESOLUTION - 1; H_pixel_num++) {

                /*

                This function polulates the initial momentum inside the s_Initial_Conditions struct

                */

                get_intitial_conditions_from_angles(s_Initial_Conditions, 
                                                    V_angle_min + V_pixel_num * Scan_Step,
                                                    H_angle_max - H_pixel_num * Scan_Step);

               /*

               Ray propagation happens here

               */

               Results_type* s_Ray_results = Propagate_ray(s_Initial_Conditions);

               /*

               Updating the visualization happens here

               */

               Renderer.Intensity_buffer[int(Renderer.texture_indexer / 3)] = s_Ray_results->Intensity[direct][I] +
                                                                              s_Ray_results->Intensity[first][I] +
                                                                              s_Ray_results->Intensity[second][I] +
                                                                              s_Ray_results->Intensity[third][I];

               Renderer.texture_indexer += 3;

               File_manager.write_image_data_to_file(s_Ray_results);

            }

        }

        auto end_time = std::chrono::high_resolution_clock::now();

        std::cout << '\n' << "Simulation finished!";
        std::cout << '\n' << "Simulation time: " << std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time);

        File_manager.close_image_output_files();

    }

    void run_simulation_mode_2(Initial_conditions_type* s_Initial_Conditions) {
        
        /*

        Create/Open the logging files

        */

        File_manager.open_image_output_files();

        /*

        Read the initial conditions from file

        */

        double J_data[500]{}, p_theta_data[500]{};
        int Data_number{};

        File_manager.get_geodesic_data(J_data, p_theta_data, &Data_number);

        for (int photon = 0; photon <= Data_number; photon += 1) {

            /*

            This function polulates the initial momentum inside the s_Initial_Conditions struct

            */

            s_Initial_Conditions->Spacetimes[e_metric]->get_initial_conditions_from_file(s_Initial_Conditions, J_data, p_theta_data, photon);

            /*
            
            Ray propagation happens here
            
            */

            Results_type* s_Ray_results = Propagate_ray(s_Initial_Conditions);

            File_manager.write_image_data_to_file(s_Ray_results);

            print_progress(photon, Data_number, true);
        }

        std::cout << '\n';

        File_manager.close_image_output_files();
    }
    
    void run_simulation_mode_3(Initial_conditions_type* s_Initial_Conditions) {

        /*

        Initialize the rendering engine (the Renderer instance must be static to not blow up the stack - the texture and intensity buffer are inside of it)

        */

        static Rendering_engine Renderer = Rendering_engine();

        std::jthread GUI_Thread(Rendering_function, &Renderer);

        auto start_time = std::chrono::high_resolution_clock::now();

        /*
        
        Perform HOTSPOT_ANIMATION_NUMBER number of simulations in order to make an animation of the hotspot
        
        */

        std::cout << '\n' << "Simulation starts..." << '\n';

        for (int hotspot_number = 0; hotspot_number <= HOTSPOT_ANIMATION_NUMBER - 1; hotspot_number++) {

            //s_Initial_Conditions->OTT_model->HOTSPOT_PHI_COORD = double(hotspot_number) / HOTSPOT_ANIMATION_NUMBER * 2 * M_PI;

            /*

            Create/Open the logging files

            */

            File_manager.open_image_output_files();

            /*

            Loop trough the viewing window

            */

            int progress = 0;

            std::cout << "Hotspot position at: " << hotspot_number << " / 4 PI" << "\n";

            for (int V_pixel_num = 0; V_pixel_num <= RESOLUTION - 1; V_pixel_num++) {

                print_progress(progress, RESOLUTION - 1, false);

                progress += 1;

                for (int H_pixel_num = 0; H_pixel_num <= RESOLUTION - 1; H_pixel_num++) {

                    /*

                    This function polulates the initial momentum inside the s_Initial_Conditions struct

                    */

                    get_intitial_conditions_from_angles(s_Initial_Conditions, V_angle_min + V_pixel_num * Scan_Step,
                                                        H_angle_max - H_pixel_num * Scan_Step);
                    /*
            
                    Ray propagation happens here
            
                    */

                    Results_type* s_Ray_results = Propagate_ray(s_Initial_Conditions);

                    Renderer.Intensity_buffer[int(Renderer.texture_indexer / 3)] = s_Ray_results->Intensity[direct][I] +
                                                                                   s_Ray_results->Intensity[first][I] +
                                                                                   s_Ray_results->Intensity[second][I] +
                                                                                   s_Ray_results->Intensity[third][I];
                    
                    Renderer.texture_indexer += 3;

                    File_manager.write_image_data_to_file(s_Ray_results);

                }

               
            }

            Renderer.texture_indexer = 0;

            std::cout << "\n";

            File_manager.close_image_output_files();
        }

        auto end_time = std::chrono::high_resolution_clock::now();

        std::cout << "Simulation finished!" << '\n';
        std::cout << "Simulation time: " << std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time);

    }

    void run_simulation_mode_4(Initial_conditions_type* s_Initial_Conditions) {

        s_Initial_Conditions->Spacetimes[e_metric]->get_initial_conditions_from_file(s_Initial_Conditions, (double*) &X_INIT, (double*) &Y_INIT, 0);

        Results_type* s_Ray_results = Propagate_ray(s_Initial_Conditions);

        File_manager.open_log_output_file();
        File_manager.log_photon_path(s_Ray_results);
        File_manager.close_log_output_file();
    }

#endif

