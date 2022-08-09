#pragma once

int Lens(double initial_conditions[], double M, double alpha_metric, double a, double r_throat, RK45 Coeff_deriv[][6],
         RK45 Coeff_sol[], RK45 Coeff_test_sol[], double max_error, double safety_factor, double r_in, double r_out,
         bool lens_from_file, std::ofstream data[], std::ofstream momentum_data[], Spacetimes e_metric, c_Kerr Kerr_class,
         c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double r_obs     = initial_conditions[e_r];
    double theta_obs = initial_conditions[e_theta];
    double phi_obs   = initial_conditions[e_phi];
    double J         = initial_conditions[3];
    double p_theta_0 = initial_conditions[e_p_theta];
    double p_r_0     = initial_conditions[e_p_r];

    double State_vector[6] = { r_obs, theta_obs, phi_obs, 0 , p_theta_0, p_r_0 };
    double State_vector_test[6]{}, Old_state[6]{};

    double Flux{}, redshift{}, Delta_phi{}, Image_coordiantes[3]{};

    int const Vector_size = sizeof(State_vector) / sizeof(double);
    int const RK45_size = 7;

    for (int vector_indexer = 0; vector_indexer <= Vector_size - 1; vector_indexer += 1) {

        State_vector_test[vector_indexer] = State_vector[vector_indexer];
        Old_state[vector_indexer] = State_vector[vector_indexer];

    }

    double r2[3]{}, r1[3]{
                          State_vector[e_r] * cos(State_vector[e_phi] + State_vector[e_phi_FD]) * sin(State_vector[e_theta]),
                          State_vector[e_r] * sin(State_vector[e_phi] + State_vector[e_phi_FD]) * sin(State_vector[e_theta]),
                          State_vector[e_r] * cos(State_vector[e_theta])
                         };

    double inter_State_vector[RK45_size * Vector_size]{}, Derivatives[RK45_size * Vector_size]{}, error[Vector_size]{};

    int integration_count, max_integration_count;

    integration_count = 0;
    max_integration_count = 7600000;

    double step = INIT_STEPSIZE;
    double affine_parameter = 0;

    int equator_crossings = 0;

    bool continue_integration = false;
    bool outside_disc = true;
    bool found_disc = false;


    while (integration_count < max_integration_count) {

        int iteration = 0;

        if (continue_integration == false) {

            for (int vector_indexer = 0; vector_indexer <= Vector_size - 1; vector_indexer += 1) {

                State_vector[vector_indexer] = Old_state[vector_indexer];

            }
        }

        while (iteration <= RK45_size - 1) { //runs trough the EOM evaluations in-between t and t + step

            for (int vector_indexer = 0; vector_indexer <= Vector_size - 1; vector_indexer += 1) { //runs trough the state vector components

                inter_State_vector[vector_indexer + iteration * Vector_size] = State_vector[vector_indexer];

                for (int derivative_indexer = 0; derivative_indexer <= iteration - 1; derivative_indexer += 1) { //runs trough tough the Dormand-Prince coeficients matrix and adds on the contributions from the derivatives at the points between t and t + step;

                    inter_State_vector[vector_indexer + iteration * Vector_size] += -step * Coeff_deriv[iteration][derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * Vector_size];

                }

            }

            get_EOM(e_metric, inter_State_vector, J, Derivatives, iteration, r_throat, a, alpha_metric, M,
                    Kerr_class, RBH_class, Wormhole_class);

            iteration += 1;

        }

        for (int vector_indexer = 0; vector_indexer <= Vector_size - 1; vector_indexer += 1) {

            State_vector_test[vector_indexer] = State_vector[vector_indexer];

            for (int derivative_indexer = 0; derivative_indexer <= iteration - 1; derivative_indexer += 1) {

                State_vector[vector_indexer] += -step * Coeff_sol[derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * Vector_size];
                State_vector_test[vector_indexer] += -step * Coeff_test_sol[derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * Vector_size];

            }

            error[vector_indexer] = State_vector[vector_indexer] - State_vector_test[vector_indexer];

        }

        if (my_max(error) < max_error)
        {
            step = 0.8 * step * pow(max_error / (my_max(error) + safety_factor), 0.2);

            double z = State_vector[e_r] * cos(State_vector[e_theta]);

            bool near_disk = z * z < 0.1 * 0.1 &&
                             State_vector[e_r] * State_vector[e_r] < 1.1 * 1.1 * r_out * r_out;

            if (near_disk)
            {
                step = 0.06 * step * pow(max_error / (my_max(error) + safety_factor), 0.2);
            }

            integration_count += 1;
            affine_parameter += step;

            continue_integration = true;

        }
        else
        {

            step = 0.8 * step * pow(max_error / (my_max(error) + safety_factor), 0.25);
            continue_integration = false;

        }

        if (continue_integration == true) {

            if (integration_count == 1) {

                found_disc = false;

                equator_crossings = 1;

                r2[x] = State_vector[e_r] * cos(State_vector[e_phi] + State_vector[e_phi_FD]) * sin(State_vector[e_theta]);
                r2[y] = State_vector[e_r] * sin(State_vector[e_phi] + State_vector[e_phi_FD]) * sin(State_vector[e_theta]);
                r2[z] = State_vector[e_r] * cos(State_vector[e_theta]);

                double photon_tangent[3] = { r1[x] - r2[x], r1[y] - r2[y], r1[z] - r2[z] };
                double photon_LOS_parameter = -dot_product(r1, r1) / dot_product(r1, photon_tangent);
                double obs_plane_intersection[3] = { r1[x] + photon_LOS_parameter * photon_tangent[x],
                                                     r1[y] + photon_LOS_parameter * photon_tangent[y],
                                                     r1[z] + photon_LOS_parameter * photon_tangent[z] };

                Rorate_to_obs_plane(theta_obs, phi_obs, obs_plane_intersection, Image_coordiantes);

                Flux = 0;
                redshift = 1e-10;
                Delta_phi = 0;

            }

            if (cos(State_vector[e_theta]) * cos(Old_state[e_theta]) < 0) {

                equator_crossings += 1;

                bool inside_disk = State_vector[e_r] * State_vector[e_r] > r_in * r_in &&
                                   State_vector[e_r] * State_vector[e_r] < r_out * r_out &&
                                   outside_disc == true;

                if (inside_disk) {

                    outside_disc = false;

                    redshift = Redshift(e_metric, J, M, r_throat, a, alpha_metric, State_vector[e_r], State_vector[e_theta],
                                               r_obs, theta_obs, Kerr_class, RBH_class, Wormhole_class);

                    Flux = get_flux(e_metric, M, r_throat, a, alpha_metric, State_vector[e_r], r_in, State_vector[e_theta],
                                           Kerr_class, RBH_class, Wormhole_class);
                    if (e_metric != Kerr) {

                        Delta_phi = fabs(State_vector[e_phi]) + J_0_correction_function(State_vector[e_theta], theta_obs);

                    }
                    else {

                        Delta_phi = (equator_crossings - 1) * M_PI;

                    }

                    write_to_file(Image_coordiantes, redshift, Flux, State_vector, alpha_metric, J,
                                  Delta_phi, lens_from_file, data, momentum_data);

                    found_disc = true;
                }

            }
            else {

                outside_disc = true;

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

                return -1;

            }

            if (terminate_integration) {
 
                if (found_disc == false) {

                    write_to_file(Image_coordiantes, redshift, Flux, State_vector, alpha_metric, J,
                                  Delta_phi, lens_from_file, data, momentum_data);

                }

                integration_count = 0;

                break;

            }

        }

    }

    return 0;
}