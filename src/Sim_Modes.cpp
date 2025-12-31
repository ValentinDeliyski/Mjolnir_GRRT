#include "Sim_Modes.h"

void static print_progress(int current, int max, bool lens_from_file) {

    int current_digits = 1;

    if (current != 0) {

        current_digits = int(floor(log10(current) + 1));

    }

    int max_digits = int(floor(log10(max) + 1));

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

void static Rendering_function(std::stop_token stop_token, Rendering_engine* Renderer, Initial_conditions_type* p_Init_conditions) {

    Renderer->OpenGL_init(p_Init_conditions);
    glfwSetKeyCallback(Renderer->window, Rendering_engine::Window_Callbacks::define_button_callbacks);

    while (!glfwWindowShouldClose(Renderer->window) and !stop_token.stop_requested()) {

        Renderer->renormalize_colormap();

        Renderer->update_rendering_window();

        std::this_thread::sleep_for(std::chrono::milliseconds(20));

    }

}

void static Update_render(Disk_model_enums Disk_model, Results_type* const p_Ray_results, Rendering_engine* const Renderer) {

    if (e_Novikov_Thorne == Disk_model) {

        Renderer->Intensity_buffer[int(Renderer->texture_indexer / 3)] = float(p_Ray_results->Flux_NT * pow(p_Ray_results->Redshift_NT, 4));

    }
    else {

        Renderer->Intensity_buffer[int(Renderer->texture_indexer / 3)] = float(p_Ray_results->Intensity[I]);

    }

    Renderer->texture_indexer += 3;

}

void static Zero_results_struct(Results_type* const p_Ray_results) {

    memset(p_Ray_results->Intensity, 0, e_Stokes_param_num * sizeof(double));

    memset(p_Ray_results->Final_State_Vector, 0, e_Full_state_size * sizeof(double));
    memset(p_Ray_results->Thin_Disk_State_Vector, 0, e_Dynamic_state_size * sizeof(double));

    memset(&p_Ray_results->Flux_NT, 0, sizeof(double));
    memset(&p_Ray_results->Redshift_NT, 0, sizeof(double));
    memset(&p_Ray_results->Projected_polarization_vector, 0, 2 * sizeof(double));

    memset(&p_Ray_results->Thin_Disk_State_Vector, 0, e_Dynamic_state_size * sizeof(double));

    p_Ray_results->Ray_log_struct.Log_length = 0;
    p_Ray_results->Ray_log_struct.Log_offset = 0;

    memset(&p_Ray_results->Image_Coords, 0, 2 * sizeof(double));

    memset(&p_Ray_results->Celestial_sphere_crossing_coords, 0, 4 * sizeof(double));

    p_Ray_results->NT_Disk_found = false;

}

void static Generate_Image(const Simulation_Context_type* const p_Sim_Context, Rendering_engine* const Renderer, Results_type* const p_Ray_results) {

        int& X_resolution = p_Sim_Context->p_Init_Conditions->Observer_params.resolution_x;
        int& Y_resolution = p_Sim_Context->p_Init_Conditions->Observer_params.resolution_y;

        double Y_angle_max = p_Sim_Context->p_Init_Conditions->Observer_params.y_angle_max;
        double Y_angle_min = p_Sim_Context->p_Init_Conditions->Observer_params.y_angle_min;
        double X_angle_max = p_Sim_Context->p_Init_Conditions->Observer_params.x_angle_max;
        double X_angle_min = p_Sim_Context->p_Init_Conditions->Observer_params.x_angle_min;

        if (!p_Sim_Context->p_Init_Conditions->Observer_params.Use_angular_coords) {

            Y_angle_max = atan2(p_Sim_Context->p_Init_Conditions->Observer_params.y_max, p_Sim_Context->p_Init_Conditions->Observer_params.distance);
            Y_angle_min = atan2(p_Sim_Context->p_Init_Conditions->Observer_params.y_min, p_Sim_Context->p_Init_Conditions->Observer_params.distance);
            X_angle_max = atan2(p_Sim_Context->p_Init_Conditions->Observer_params.x_max, p_Sim_Context->p_Init_Conditions->Observer_params.distance);
            X_angle_min = atan2(p_Sim_Context->p_Init_Conditions->Observer_params.x_min, p_Sim_Context->p_Init_Conditions->Observer_params.distance);

        }
        else {

            p_Sim_Context->p_Init_Conditions->Observer_params.y_max = p_Sim_Context->p_Init_Conditions->Observer_params.distance * tan(Y_angle_max);
            p_Sim_Context->p_Init_Conditions->Observer_params.y_min = p_Sim_Context->p_Init_Conditions->Observer_params.distance * tan(Y_angle_min);
            p_Sim_Context->p_Init_Conditions->Observer_params.x_max = p_Sim_Context->p_Init_Conditions->Observer_params.distance * tan(X_angle_max);
            p_Sim_Context->p_Init_Conditions->Observer_params.x_min = p_Sim_Context->p_Init_Conditions->Observer_params.distance * tan(X_angle_min);

        }

        // Having a non-even resolution means that for a symmetric observation window, there will be a vertical line of pixels that coorespond 
        // to photons with zero azimuthal angular momentum. In that case the behavior of the theta and phi coordinates swap, and theta becomes unbounded.
        // This causes sin(theta) to take on negative values, which breaks the radiative transfer. In such cases I offset the observation window slightly
        // to cirmumvent this little hickup.
        if (p_Sim_Context->p_Init_Conditions->Observer_params.resolution_x % 2 != 0) {

            X_angle_max += X_angle_max * 1e-5;

        }

        double X_scan_step = (X_angle_max - X_angle_min) / (X_resolution - 1);
        double Y_scan_step = (Y_angle_max - Y_angle_min) / (Y_resolution - 1);

        /*

        Create/Open the logging files

        */

        p_Sim_Context->File_manager->create_output_file();

        /*

        Loop trough the viewing window

        */

        auto start_time = std::chrono::high_resolution_clock::now();
        int progress = 0;

        std::cout << '\n' << "Generating image for " << p_Sim_Context->p_Init_Conditions->File_manager_params.Simulation_name << "...\n";

        for (int V_pixel_num = 0; V_pixel_num < Y_resolution; V_pixel_num++) {

            if (p_Sim_Context->p_Init_Conditions->Print_to_console) { print_progress(progress, Y_resolution - 1, false); }

            progress += 1;

            for (int H_pixel_num = 0; H_pixel_num < X_resolution; H_pixel_num++) {

                /*  ------------------ This function polulates the initial momentum inside the s_Initial_Conditions struct ------------------ */
                get_intitial_conditions_from_angles(p_Sim_Context->p_Init_Conditions,
                                                    Y_angle_min + V_pixel_num * Y_scan_step,
                                                    X_angle_max - H_pixel_num * X_scan_step);
                
                /* ------------------  Ray propagation happens here ------------------ */
                Propagate_ray(p_Sim_Context, p_Ray_results);
                
                /* ------------------ Updating the visualization happens here ------------------ */
                Update_render(p_Sim_Context->p_Init_Conditions->Disk_params.e_Disk_model, p_Ray_results, Renderer);

                /* ------------------ Results logging happens here ------------------ */

                p_Sim_Context->File_manager->open_image_output_file();
                p_Sim_Context->File_manager->write_image_data_to_file(p_Ray_results);
                p_Sim_Context->File_manager->close_image_output_file();

                /* The final results must be manually set to 0s because the Ray_results struct is STATIC (and in an outer scope), 
                   and therefore not automatically reinitialized to 0s. I have to manually do it. */
                Zero_results_struct(p_Ray_results);

            }

        }

        auto end_time = std::chrono::high_resolution_clock::now();

        std::cout << '\n' << "Image Generation for " << p_Sim_Context->p_Init_Conditions->File_manager_params.Simulation_name << " Finished!" << 
                     '\n' << "Simulation time: " << std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time) << "\n";

    }

void run_image_generation(const Simulation_Context_type* const p_Sim_Context, Results_type* const p_Ray_results) {
       
        /*

        Initialize the rendering engine (the Renderer instance must be static to not blow up the stack - the texture and intensity buffers are inside of it)

        */

        static Rendering_engine Renderer = Rendering_engine();
        std::jthread GUI_Thread(Rendering_function, &Renderer, p_Sim_Context->p_Init_Conditions);

        /*
        
        Run the main ray-tracer "kernel"
        
        */

        Generate_Image(p_Sim_Context, &Renderer, p_Ray_results);

        GUI_Thread.request_stop();

        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        Renderer.Free_memory();

    }

void run_geodesic_sweep(const Simulation_Context_type* const p_Sim_Context, Results_type* const p_Ray_results) {

    /*

    Read the initial conditions from file

    */

    double p_phi_data[500]{}, p_theta_data[500]{};

    p_Sim_Context->File_manager->get_geodesic_data(p_phi_data, p_theta_data);

    /*

    Create/Open the logging files

    */

    p_Sim_Context->File_manager->create_output_file();
    p_Sim_Context->File_manager->open_image_output_file();

    for (int photon_idx = 0; photon_idx <= p_Sim_Context->File_manager->sim_mode_2_ray_number - 1; photon_idx += 1) {

        /*

        This function polulates the initial momentum inside the s_Initial_Conditions struct

        */

        get_initial_conditions_from_image_coords(p_Sim_Context->p_Init_Conditions, p_phi_data[photon_idx], p_theta_data[photon_idx]);

        /*

        Ray propagation happens here

        */

        Propagate_ray(p_Sim_Context, p_Ray_results);

        /*
        
        Results logging happens here
        
        */

        p_Sim_Context->File_manager->write_image_data_to_file(p_Ray_results);

        /* The final results must be manually set to 0s because the Ray_results struct is STATIC (and in an outer scope),
           and therefore not automatically reinitialized to 0s. I have to manually do it. */
        Zero_results_struct(p_Ray_results);

        print_progress(photon_idx, p_Sim_Context->File_manager->sim_mode_2_ray_number - 1, true);

    }

    p_Sim_Context->File_manager->close_image_output_file();

    std::cout << '\n';


}

void make_geodesic_log(const Simulation_Context_type* const p_Sim_Context, Results_type* const p_Ray_results) {

    double& X_init = p_Sim_Context->p_Init_Conditions->Sim_mode_3_X_init;
    double& Y_init = p_Sim_Context->p_Init_Conditions->Sim_mode_3_Y_init;

    get_initial_conditions_from_image_coords(p_Sim_Context->p_Init_Conditions, X_init, Y_init);

    Propagate_ray(p_Sim_Context, p_Ray_results);

    p_Sim_Context->File_manager->open_image_output_file();
    p_Sim_Context->File_manager->log_photon_path(p_Ray_results);
    p_Sim_Context->File_manager->close_log_output_file();
}


