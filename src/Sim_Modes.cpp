#pragma once

#ifndef SIM_MODES

    #define SIM_MODES

    #include "IO_files.h"
    #include "Constants.h"
    #include "General_GR_functions.h"
    #include "Rendering_Engine.h"

    extern std::vector<c_Spacetime_Base*> Spacetimes;
    extern bool Normalizing_colormap;

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
            else if (!Normalizing_colormap) {

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

	void run_simulation_mode_2(Initial_conditions_type* s_Initial_Conditions, std::ofstream data[], std::ofstream momentum_data[]) {
        
        Normalizing_colormap = false;

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

            Spacetimes[e_metric]->get_initial_conditions_from_file(s_Initial_Conditions, J_data, p_theta_data, photon);

            Lens(s_Initial_Conditions, data, momentum_data);

            print_progress(photon, Data_number, true, false);
        }

        std::cout << '\n';


	}

    void run_simulation_mode_1(Initial_conditions_type* s_Initial_Conditions, std::ofstream data[], std::ofstream momentum_data[]) {

        GLFWwindow* window = OpenGL_init(H_angle_max / V_angle_max);
        glfwSetKeyCallback(window, Window_Callbacks::define_button_callbacks);

        auto start_time = std::chrono::high_resolution_clock::now();

        /*

        Do one scan line in the middle of the image to find the maximum intensity for use in the colormap

        */

        std::cout << "Initial y = 0 scan to normalize the colormap..." << '\n';

        for (int pixel_num = 0; pixel_num <= RESOLUTION - 1; pixel_num++) {

            get_intitial_conditions_from_angles(s_Initial_Conditions, 0, H_angle_max - pixel_num * Scan_Step);

            Lens(s_Initial_Conditions, data, momentum_data);

            print_progress(pixel_num, RESOLUTION - 1, false, Normalizing_colormap);

        }

        Normalizing_colormap = false;

        /*

        Loop trough the viewing window

        */

        int progress = 0;

        std::cout << '\n' << "Simulation starts..." << '\n';

        for (int V_pixel_num = 0; V_pixel_num <= RESOLUTION - 1; V_pixel_num++) {

            update_rendering_window(window, V_angle_max / H_angle_max);
            print_progress(progress, RESOLUTION - 1, false, Normalizing_colormap);

            progress += 1;

            for (int H_pixel_num = 0; H_pixel_num <= RESOLUTION - 1; H_pixel_num++) {

                /*

                This function polulates the initial momentum inside the s_Initial_Conditions struct

                */

                get_intitial_conditions_from_angles(s_Initial_Conditions, V_angle_min + V_pixel_num * Scan_Step,
                                                    H_angle_max - H_pixel_num * Scan_Step);

                Lens(s_Initial_Conditions, data, momentum_data);

            }

        }

        auto end_time = std::chrono::high_resolution_clock::now();

        std::cout << '\n' << "Simulation finished!" << '\n';

        std::cout << "Simulation time: " << std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time);

        while (!glfwWindowShouldClose(window)) {

            update_rendering_window(window, V_angle_max / H_angle_max);

        }

    }

#endif

