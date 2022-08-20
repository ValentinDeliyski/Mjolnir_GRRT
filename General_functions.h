#pragma once

int Rorate_to_obs_plane(double theta_obs, double phi_obs, double Image_point[3], double rotated_Image_point[3]) {

    const double R_theta[3][3] =
    {
      {1,             0,                       0           },
      {0,  cos(M_PI_2 - theta_obs), sin(M_PI_2 - theta_obs)},
      {0, -sin(M_PI_2 - theta_obs), cos(M_PI_2 - theta_obs)}
    };

    const double R_phi[3][3] =
    {
      {cos(M_PI_2 - phi_obs), -sin(M_PI_2 - phi_obs), 0},
      {sin(M_PI_2 - phi_obs),  cos(M_PI_2 - phi_obs), 0},
      {           0,                      0,          1}
    };

    double Rotation_matrix[3][3]{};

    for (int row = 0; row <= 3 - 1; row += 1) {

        for (int column = 0; column <= 3 - 1; column += 1) {

            for (int k = 0; k <= 3 - 1; k += 1) {

                Rotation_matrix[row][column] += R_theta[row][k] * R_phi[k][column];

            }

        }
    }                
    
    for (int k = x; k <= z; k += 1) {

        rotated_Image_point[k] = 0;

    }

    for (int vector_index = 0; vector_index <= 3 - 1; vector_index += 1) {

        for (int k = 0; k <= 3 - 1; k += 1) {

            rotated_Image_point[vector_index] += Rotation_matrix[vector_index][k] * Image_point[k];

        }
    }

    return OK;
}

double my_max(double vector[]) {

    int index_max = e_Coord_Number - 1;

    for (int index = 0; index <= index_max; index += 1) {

        if (vector[index] < 0) {

            vector[index] = -1.0 * vector[index];
        }

    }

    double max = vector[0];

    for (int index = 1; index <= index_max; index += 1) {

        if (vector[index] > max) {

            max = vector[index];

        }
    }

    return max;
}

bool crossed_equatior(double State_vector[], double Old_State_Vector[]) {

    return cos(State_vector[e_theta]) * cos(Old_State_Vector[e_theta]) < 0;

}

bool Inside_disc(double State_Vector[], double Old_State_Vector[], double r_in, double r_out) {

    bool inside_disc = State_Vector[e_r] * State_Vector[e_r] > r_in * r_in &&
                       State_Vector[e_r] * State_Vector[e_r] < r_out * r_out &&
                       crossed_equatior(State_Vector, Old_State_Vector);


    return inside_disc;
}

double dot_product(double vector_1[3], double vector_2[3]) {

    return vector_1[0] * vector_2[0] + vector_1[1] * vector_2[1] + vector_1[2] * vector_2[2];

}

int get_metric(Spacetimes e_metric, double metric[4][4], double* N_metric, double* omega, double M, 
               double r_throat, double a, double alpha_metric, double r, double theta,
               c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    switch (e_metric) {

    case Kerr:

        Kerr_class.metric(metric, N_metric, omega, M, a, r, theta);

        break;

    case Wormhole:

        Wormhole_class.metric(metric, N_metric, omega, M, r_throat, a, alpha_metric, r, theta);

        break;

    case Reg_Black_Hole:

        RBH_class.metric(metric, N_metric, omega, M, a, r, theta);

        break;

    default:

        std::cout << "Wrong metric!" << '\n';

        return -1;

    }

    return OK;
}

int get_metric_fist_derivatives(Spacetimes e_metric, double dr_metric[4][4], double* dr_N_metric, double* dr_omega,
                                double M, double r_throat, double a, double alpha_metric, double r, double theta,
                                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    switch (e_metric) {

    case Kerr:

        Kerr_class.metric_first_derivatives(Kerr_class, dr_metric, dr_N_metric, dr_omega, M, a, r, theta);

        break;

    case Wormhole:

        Wormhole_class.metric_first_derivatives(Wormhole_class, dr_metric, dr_N_metric, dr_omega, M, r_throat, a, alpha_metric, r, theta);

        break;

    case Reg_Black_Hole:

        RBH_class.metric_first_derivatives(RBH_class, dr_metric, dr_N_metric, dr_omega, M, a, r, theta);

        break;

    default:

        std::cout << "Wrong metric!" << '\n';

        return -1;

    }

    return OK;

}

int get_metric_second_derivatives(Spacetimes e_metric, double d2r_metric[4][4], double* d2r_N_metric, double* d2r_omega,
                                  double M, double r_throat, double a, double alpha_metric, double r, double theta,
                                  c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    switch (e_metric) {

    case Kerr:

        Kerr_class.metric_second_derivatives(Kerr_class, d2r_metric, d2r_N_metric, d2r_omega, M, a, r, theta);

        break;

    case Wormhole:

        Wormhole_class.metric_second_derivatives(Wormhole_class, d2r_metric, d2r_N_metric, d2r_omega, M, r_throat, a, alpha_metric, r, theta);

        break;

    case Reg_Black_Hole:

        RBH_class.metric_second_derivatives(RBH_class, d2r_metric, d2r_N_metric, d2r_omega, M, a, r, theta);

        break;

    default:

        std::cout << "Wrong metric!" << '\n';

        return -1;

    }

    return OK;

}

int get_intitial_conditions_from_angles(double* J, double* p_theta, double* p_r, double metric[4][4],
                                        double V_angle, double H_angle) {

    double g2, gamma, ksi, L_z, E;

    g2 = pow(metric[0][3], 2) - metric[0][0] * metric[3][3];
    ksi = sqrt(metric[3][3] / g2);
    gamma = -metric[0][3] / metric[3][3] * ksi;

    L_z = sqrt(metric[3][3]) * sin(H_angle + 2 * M_PI) * cos(V_angle);
    E = (1 + gamma * L_z) / ksi;

    *J = L_z / E;
    *p_theta = sqrt(metric[2][2]) * sin(V_angle) / E;
    *p_r = sqrt(metric[1][1]) * cos(H_angle + 2 * M_PI) * cos(V_angle) / E;

    return OK;
}

int get_initial_conditions_from_file(Spacetimes e_metric, double* J, double J_data[], double* p_theta, double p_theta_data[],
                                     double* p_r, int photon, double r_obs, double theta_obs, double metric[4][4],
                                     double N_metric, double omega_metric, double M, double a, double alpha_metric,
                                     c_Kerr Kerr_class,c_RBH RBH_class, c_Wormhole Wormhole_class) {

    switch (e_metric) {

    case Kerr:

        Kerr_class.intitial_conditions_from_file(J, J_data, p_theta, p_theta_data, p_r, photon, r_obs, theta_obs, metric, M, a);

        break;

    case Wormhole:

        Wormhole_class.intitial_conditions_from_file(J, J_data, p_theta, p_theta_data, p_r, photon, r_obs, theta_obs, metric, N_metric, omega_metric);

        break;

    case Reg_Black_Hole:

        RBH_class.intitial_conditions_from_file(J, J_data, p_theta, p_theta_data, p_r, photon, r_obs, theta_obs, metric, M, a);

        break;
    }

    return OK;

}

double Keplerian_angular_velocity(Spacetimes e_metric, double M, double r_throat, double a,
                                  double alpha_metric, double r, double theta,
                                  c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double dr_metric[4][4], dr_N, dr_omega;

    get_metric_fist_derivatives(e_metric, dr_metric, &dr_N, &dr_omega, M, r_throat, a, alpha_metric, r, theta,
                                Kerr_class, RBH_class,Wormhole_class);

    return   (-dr_metric[0][3] + sqrt(dr_metric[0][3] * dr_metric[0][3] - dr_metric[0][0] * dr_metric[3][3])) / dr_metric[3][3];

}

double dr_Keplerian_angular_velocity(Spacetimes e_metric, double M, double r_throat, double a,
                                     double alpha_metric, double r, double theta, double Kepler,
                                     c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double dr_metric[4][4], dr_N, dr_omega;

    double d2r_metric[4][4], d2r_N, d2r_omega;

    get_metric_fist_derivatives(e_metric, dr_metric, &dr_N, &dr_omega, M, r_throat, a, alpha_metric, r, theta,
                                Kerr_class, RBH_class, Wormhole_class);

    get_metric_second_derivatives(e_metric, d2r_metric, &d2r_N, &d2r_omega, M, r_throat, a, alpha_metric, r, theta,
                                  Kerr_class, RBH_class, Wormhole_class);

    double root = sqrt(dr_metric[0][3] * dr_metric[0][3] - dr_metric[0][0] * dr_metric[3][3]);

    return  -Kepler / dr_metric[3][3] * d2r_metric[3][3] + (-d2r_metric[0][3] + 1.0 / root / 2 * (2 * dr_metric[0][3] * d2r_metric[0][3] - dr_metric[0][0] * d2r_metric[3][3] - d2r_metric[0][0] * dr_metric[3][3])) / dr_metric[3][3];

}

double Redshift(Spacetimes e_metric,double J, double M, double r_throat, double a, double alpha_metric,
                double r_source, double theta_source, double r_obs, double theta_obs,
                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double metric_source[4][4], N_source, omega_source;

    double metric_obs[4][4], N_obs, omega_obs;

    get_metric(e_metric, metric_source, &N_source, &omega_source, M, r_throat, a, alpha_metric, r_source, theta_source,
               Kerr_class, RBH_class, Wormhole_class);

    get_metric(e_metric, metric_obs, &N_obs, &omega_obs, M, r_throat, a, alpha_metric, r_obs, theta_obs,
               Kerr_class, RBH_class, Wormhole_class);

    double Kepler = Keplerian_angular_velocity(e_metric, M, r_throat, a, alpha_metric, r_source, theta_source,
                                               Kerr_class, RBH_class, Wormhole_class);

    double gamma = 1 / sqrt(-metric_source[0][0] - 2 * metric_source[0][3] * Kepler - metric_source[3][3] * Kepler * Kepler);

    double U_source[4] = { gamma, 0, 0, gamma * Kepler };

    double U_obs[4] = { 1.0 / N_obs, 0 ,0 , omega_obs / N_obs };

    return  (-U_obs[0] + U_obs[3] * J) / (-U_source[0] + U_source[3] * J);

}

double Flux_integrand(Spacetimes e_metric, double M, double r_throat, double a, double alpha_metric,
                      double r, double theta, double* metric_det, double* E_disk, double* L_disk,
                      c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double metric[4][4], N, omega;

    double dr_metric[4][4], dr_N, dr_omega;

    get_metric(e_metric, metric, &N, &omega, M, r_throat, a, alpha_metric, r, theta,
               Kerr_class, RBH_class, Wormhole_class);

    get_metric_fist_derivatives(e_metric, dr_metric, &dr_N, &dr_omega, M, r_throat, a, alpha_metric, r, theta,
                                Kerr_class, RBH_class, Wormhole_class);

    double Kepler    = Keplerian_angular_velocity(e_metric, M, r_throat, a, alpha_metric, r, theta,
                                                  Kerr_class, RBH_class, Wormhole_class);
    double dr_Kepler = dr_Keplerian_angular_velocity(e_metric, M, r_throat, a, alpha_metric, r, theta, Kepler,
                                                     Kerr_class, RBH_class, Wormhole_class);

    *metric_det = metric[1][1] * (metric[0][3] * metric[0][3] - metric[0][0] * metric[3][3]);

    double root = sqrt(-metric[0][0] - 2 * metric[0][3] * Kepler - metric[3][3] * Kepler * Kepler);
    double dr_root = (-dr_metric[0][0] - 2 * (dr_metric[0][3] * Kepler + metric[0][3] * dr_Kepler) - dr_metric[3][3] * Kepler * Kepler - 2 * metric[3][3] * Kepler * dr_Kepler);
    *E_disk = -(metric[0][0] + metric[0][3] * Kepler) / root;
    *L_disk = (metric[3][3] * Kepler + metric[0][3]) / root;
    double dr_L = (dr_metric[3][3] * Kepler + metric[3][3] * dr_Kepler + dr_metric[0][3]) / root - *L_disk / root / root / 2 * dr_root;

    return (*E_disk - Kepler * *L_disk) * dr_L;

}

double solve_integral(Spacetimes e_metric, double M, double r_throat, double a, double alpha_metric,
                      double lower_bound, double upper_bound, double theta, double tolerance,
                      c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double metric_det, E_disk, L_disk;

    double mid_point = (lower_bound + upper_bound) / 2;
    double left_mid_point = (lower_bound + mid_point) / 2;
    double right_mid_point = (mid_point + upper_bound) / 2;

    double F_lower_bound = Flux_integrand(e_metric, M, r_throat, a, alpha_metric, lower_bound, theta, &metric_det, &E_disk, &L_disk,
                                          Kerr_class, RBH_class, Wormhole_class);
    double F_mid_point = Flux_integrand(e_metric, M, r_throat, a, alpha_metric, mid_point, theta, &metric_det, &E_disk, &L_disk,
                                        Kerr_class, RBH_class, Wormhole_class);
    double F_upper_bound = Flux_integrand(e_metric, M, r_throat, a, alpha_metric, upper_bound, theta, &metric_det, &E_disk, &L_disk,
                                          Kerr_class, RBH_class, Wormhole_class);
 
    double F_left_mid = Flux_integrand(e_metric, M, r_throat, a, alpha_metric, left_mid_point, theta, &metric_det, &E_disk, &L_disk,
                                       Kerr_class, RBH_class, Wormhole_class);
    double F_right_mid = Flux_integrand(e_metric, M, r_throat, a, alpha_metric, right_mid_point, theta, &metric_det, &E_disk, &L_disk,
                                        Kerr_class, RBH_class, Wormhole_class);

    double S_left = (mid_point - lower_bound) / 6 * (F_lower_bound + 4 * F_left_mid + F_mid_point);
    double S_right = (upper_bound - mid_point) / 6 * (F_mid_point + 4 * F_right_mid + F_upper_bound);

    double S_2 = S_left + S_right;
    double S_1 = (upper_bound - lower_bound) / 6 * (F_lower_bound + 4 * F_mid_point + F_upper_bound);

    double error;

    if (S_2 >= S_1) {

        error = S_2 - S_1;

    }
    else {

        error = S_1 - S_2;

    }

    double integral;

    if (error < 15 * tolerance) {

        integral = S_2 + (S_2 - S_1) / 15;


    }
    else {

        double L_value = solve_integral(e_metric, M, r_throat, a, alpha_metric, lower_bound, mid_point, theta, tolerance / 2,
                                        Kerr_class, RBH_class, Wormhole_class);
        double R_value = solve_integral(e_metric, M, r_throat, a, alpha_metric, mid_point, upper_bound, theta, tolerance / 2,
                                        Kerr_class, RBH_class, Wormhole_class);

        integral = L_value + R_value;

    }

    return integral;
}

double get_flux(Spacetimes e_metric, double M, double r_throat, double a,
                double alpha_metric, double r, double r_in, double theta,
                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double metric_det{}, E_disk{}, L_disk{};

    double Kepler = Keplerian_angular_velocity(e_metric, M, r_throat, a, alpha_metric, r, theta,
                                               Kerr_class, RBH_class, Wormhole_class);
    double dr_Kepler = dr_Keplerian_angular_velocity(e_metric, M, r_throat, a, alpha_metric, r, theta, Kepler,
                                                     Kerr_class, RBH_class, Wormhole_class);

    Flux_integrand(e_metric, M, r_throat, a, alpha_metric, r, theta, &metric_det, &E_disk, &L_disk, 
                  Kerr_class, RBH_class, Wormhole_class);

    double Flux_coeff = -dr_Kepler / ((E_disk - Kepler * L_disk) * (E_disk - Kepler * L_disk)) / (4 * M_PI * sqrt(metric_det));

    if (e_metric == Wormhole) {

        r = sqrt(r * r + r_throat * r_throat);

    }
        
    double Flux_integral = solve_integral(e_metric, M, r_throat, a, alpha_metric, r_in, r, theta, INTEGRAL_ACCURACY,
                                          Kerr_class, RBH_class, Wormhole_class);

    return Flux_coeff * Flux_integral;

}

int get_EOM(Spacetimes e_metric, double inter_State_vector[7 * 6], double J, double Derivatives[7 * 6],
            int iteration, double r_throat, double a, double metric_parameter, double M,
            c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    switch (e_metric) {

    case Kerr:

        Kerr_class.EOM(inter_State_vector, J, Derivatives, iteration, a, M);

        break;

    case Wormhole:

        Wormhole_class.EOM(inter_State_vector, J, Derivatives, iteration, r_throat, a, metric_parameter);

        break;

    case Reg_Black_Hole:

        RBH_class.EOM(inter_State_vector, J, Derivatives, iteration, metric_parameter, M);

        break;

    
    default:

        std::cout << "Wrong metric!" << '\n';

        return ERROR;
    }

    return OK;
}

Return_Value_enums Evaluate_Disk_Model(Disk_Models e_Disk_Model, double State_vector[], Spacetimes e_metric, double J, double M, double r_throat,
                        double a, double metric_parameter, double r_obs, double theta_obs, double r_in, double r_out,
                        double* redshift, double* Flux, c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    switch (e_Disk_Model) {

        case Novikov_Thorne:

            *redshift = Redshift(e_metric, J, M, r_throat, a, metric_parameter, State_vector[e_r], State_vector[e_theta],
                                 r_obs, theta_obs, Kerr_class, RBH_class, Wormhole_class);

            *Flux = get_flux(e_metric, M, r_throat, a, metric_parameter, State_vector[e_r], r_in, State_vector[e_theta],
                             Kerr_class, RBH_class, Wormhole_class);

            break;

        case Optically_Thin_Toroidal:

            break;

        default:

            std::cout << "Wrong Disc Model!" << '\n';

            return ERROR;

    }

    return OK;

}