#pragma once

Return_Value_enums Lens(double initial_conditions[], double M, double metric_parameter, double a, double r_throat, double r_in, double r_out, bool lens_from_file,
                        std::ofstream data[], std::ofstream momentum_data[], Spacetimes e_metric, c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double r_obs     = initial_conditions[e_r];
    double theta_obs = initial_conditions[e_theta];
    double phi_obs   = initial_conditions[e_phi];
    double J         = initial_conditions[3];
    double p_theta_0 = initial_conditions[e_p_theta];
    double p_r_0     = initial_conditions[e_p_r];

    double State_vector[6] = { r_obs, theta_obs, phi_obs, 0 , p_theta_0, p_r_0 };
    double State_vector_test[6]{}, Old_state[6]{};

    double Flux{}, redshift{}, Image_coordiantes[3]{};

    int const Vector_size = sizeof(State_vector) / sizeof(double);

    for (int vector_indexer = e_r; vector_indexer <= e_p_r; vector_indexer += 1) {

        State_vector_test[vector_indexer] = State_vector[vector_indexer];
        Old_state[vector_indexer] = State_vector[vector_indexer];

    }

    double r2[3]{}, r1[3]{
                          State_vector[e_r] * cos(State_vector[e_phi] + State_vector[e_phi_FD]) * sin(State_vector[e_theta]),
                          State_vector[e_r] * sin(State_vector[e_phi] + State_vector[e_phi_FD]) * sin(State_vector[e_theta]),
                          State_vector[e_r] * cos(State_vector[e_theta])
                         };

    double inter_State_vector[RK45_size * Vector_size]{}, Derivatives[RK45_size * Vector_size]{}, error[Vector_size]{};

    int integration_count{}, Image_Order{};

    double step = INIT_STEPSIZE;
    double affine_parameter{};

    bool continue_integration = false;
    bool found_disc[ORDER_NUM]{};

    Return_Value_enums RK45_Status = OK, Disc_model_status = OK;

    while (RK45_Status == OK && integration_count < MAX_INTEGRATION_COUNT) {

        RK45_Status = RK45_EOM(State_vector, Derivatives, &step, r_throat, a, metric_parameter, M, J, Kerr_class, e_metric, RBH_class, Wormhole_class, &continue_integration);

        if (continue_integration == true) {

            if (integration_count == 1) {

                for (int index = 0; index <= ORDER_NUM - 1; index++) {

                    found_disc[index] = false;

                }

                redshift = 0;
                Flux     = 0;

                Image_Order = direct;

                r2[x] = State_vector[e_r] * cos(State_vector[e_phi] + State_vector[e_phi_FD]) * sin(State_vector[e_theta]);
                r2[y] = State_vector[e_r] * sin(State_vector[e_phi] + State_vector[e_phi_FD]) * sin(State_vector[e_theta]);
                r2[z] = State_vector[e_r] * cos(State_vector[e_theta]);

                double photon_tangent[3]         = { r1[x] - r2[x], r1[y] - r2[y], r1[z] - r2[z] };
                double photon_LOS_parameter      = -dot_product(r1, r1) / dot_product(r1, photon_tangent);
                double obs_plane_intersection[3] = { r1[x] + photon_LOS_parameter * photon_tangent[x],
                                                     r1[y] + photon_LOS_parameter * photon_tangent[y],
                                                     r1[z] + photon_LOS_parameter * photon_tangent[z] };

                Rorate_to_obs_plane(theta_obs, phi_obs, obs_plane_intersection, Image_coordiantes);

            }
                                    
            if (Inside_disc(State_vector, Old_state, r_in, r_out)) {

                Disc_model_status = Evaluate_Disk_Model(Novikov_Thorne, State_vector, Kerr, J, M, r_throat, a, metric_parameter, r_obs, theta_obs,
                                                        r_in, r_out, &redshift, &Flux, Kerr_class, RBH_class, Wormhole_class);

                if (Disc_model_status != OK) {

                    std::cout << "Error Evaluating Disc Model!" << '\n';

                    return ERROR;

                }

               write_to_file(Image_coordiantes, redshift, Flux, State_vector, metric_parameter, J,
                             Image_Order, lens_from_file, data, momentum_data);

               found_disc[Image_Order] = true;

            }

            if (crossed_equatior(State_vector, Old_state)) {

                Image_Order += 1;

            }

            for (int vector_indexer = 0; vector_indexer <= Vector_size - 1; vector_indexer += 1) {

                Old_state[vector_indexer] = State_vector[vector_indexer];

            }

            bool scatter            = State_vector[e_r] > r_out && Derivatives[e_r] < 0;
            bool scatter_other_side = State_vector[e_r] < -sqrt(r_out * r_out + r_throat * r_throat);

            bool hit_horizon_kerr = State_vector[e_r] - Kerr_class.get_r_horizon() < 0.05;
            bool hit_horizon_RBH  = State_vector[e_r] - RBH_class.get_r_horizon() < 0.05;

            bool pass_trough_throat = State_vector[e_r] < 0;

            bool terminate_integration = scatter;

            switch (e_metric) {

            case Wormhole:

                scatter = State_vector[e_r] > sqrt(r_out * r_out + r_throat * r_throat) && Derivatives[0] < 0;

                terminate_integration = scatter || scatter_other_side;

                break;

            case Kerr:
                
                terminate_integration = scatter || hit_horizon_kerr;

                break;

            case Reg_Black_Hole:

                terminate_integration = scatter || hit_horizon_RBH;

                break;

            default:

                std::cout << "Wrong metric!" << '\n';

                return ERROR;

            }

            if (terminate_integration) {

                for (int Image_Order_Scan = direct; Image_Order_Scan <= third; Image_Order_Scan +=1) {

                    if (found_disc[Image_Order_Scan] == false && lens_from_file == false) {

                        write_to_file(Image_coordiantes, 0., 0., State_vector, metric_parameter, J,
                                      Image_Order_Scan, lens_from_file, data, momentum_data);

                    }
                }

                integration_count = 0;

                break;

            }

            integration_count += 1;

        }

    }

    return RK45_Status;

}
