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

    int index_max = e_State_Number - 1;

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

double dot_product(double vector_1[3], double vector_2[3]) {

    return vector_1[0] * vector_2[0] + vector_1[1] * vector_2[1] + vector_1[2] * vector_2[2];

}

int invert_metric(double inv_metric[4][4], double metric[4][4]) {

    double g2 = metric[0][3] * metric[0][3] - metric[3][3] * metric[0][0];

    inv_metric[0][0] = -metric[3][3] / g2;
    inv_metric[0][3] = metric[0][3] / g2;
    inv_metric[3][0] = inv_metric[0][3];
    inv_metric[1][1] = 1. / metric[1][1];
    inv_metric[2][2] = 1. / metric[2][2];
    inv_metric[3][3] = -metric[0][0] / g2;

    return OK;

}

double get_ISCO(Spacetimes e_metric, c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    switch (e_metric) {

    case Kerr:

        return Kerr_class.get_ISCO();

    case Wormhole:

        return Wormhole_class.get_ISCO();

    case Reg_Black_Hole:

        return RBH_class.get_ISCO();


    }

}

int get_metric(Spacetimes e_metric, double metric[4][4], double* N_metric, double* omega, double r, double theta,
               c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    switch (e_metric) {

    case Kerr:

        Kerr_class.metric(metric, N_metric, omega, r, theta);

        break;

    case Wormhole:

        Wormhole_class.metric(metric, N_metric, omega, r, theta);

        break;

    case Reg_Black_Hole:

        RBH_class.metric(metric, N_metric, omega, r, theta);

        break;

    default:

        std::cout << "Wrong metric!" << '\n';

        return -1;

    }

    return OK;
}

int get_metric_fist_derivatives(Spacetimes e_metric, double dr_metric[4][4], double* dr_N_metric, double* dr_omega,
                                double r, double theta, c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {
    
    switch (e_metric) {

    case Kerr:

        Kerr_class.metric_first_derivatives(Kerr_class, dr_metric, dr_N_metric, dr_omega, r, theta);

        break;

    case Wormhole:

        Wormhole_class.metric_first_derivatives(Wormhole_class, dr_metric, dr_N_metric, dr_omega, r, theta);

        break;

    case Reg_Black_Hole:

        RBH_class.metric_first_derivatives(RBH_class, dr_metric, dr_N_metric, dr_omega, r, theta);

        break;

    default:

        std::cout << "Wrong metric!" << '\n';

        return ERROR;

    }

    return OK;

}

int get_metric_second_derivatives(Spacetimes e_metric, double d2r_metric[4][4], double* d2r_N_metric, double* d2r_omega,
                                  double r, double theta, c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    switch (e_metric) {

    case Kerr:

        Kerr_class.metric_second_derivatives(Kerr_class, d2r_metric, d2r_N_metric, d2r_omega, r, theta);

        break;

    case Wormhole:

        Wormhole_class.metric_second_derivatives(Wormhole_class, d2r_metric, d2r_N_metric, d2r_omega, r, theta);

        break;

    case Reg_Black_Hole:

        RBH_class.metric_second_derivatives(RBH_class, d2r_metric, d2r_N_metric, d2r_omega, r, theta);

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
                                     double N_metric, double omega_metric,
                                     c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    switch (e_metric) {

    case Kerr:

        Kerr_class.intitial_conditions_from_file(J, J_data, p_theta, p_theta_data, p_r, photon, r_obs, theta_obs, metric);

        break;

    case Wormhole:

        Wormhole_class.intitial_conditions_from_file(J, J_data, p_theta, p_theta_data, p_r, photon, r_obs, theta_obs, metric, N_metric, omega_metric);

        break;

    case Reg_Black_Hole:

        RBH_class.intitial_conditions_from_file(J, J_data, p_theta, p_theta_data, p_r, photon, r_obs, theta_obs, metric);

        break;
    }

    return OK;

}

double Keplerian_angular_velocity(Spacetimes e_metric, double r, double theta,
                                  c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double dr_metric[4][4], dr_N, dr_omega;

    if (r < get_ISCO(e_metric, Kerr_class, RBH_class, Wormhole_class)) {

        r = get_ISCO(e_metric, Kerr_class, RBH_class, Wormhole_class);
        theta = M_PI_2;

    }

    get_metric_fist_derivatives(e_metric, dr_metric, &dr_N, &dr_omega, r, theta,
                                Kerr_class, RBH_class, Wormhole_class);

    return (-dr_metric[0][3] + sqrt(dr_metric[0][3] * dr_metric[0][3] - dr_metric[0][0] * dr_metric[3][3])) / dr_metric[3][3];

}

double dr_Keplerian_angular_velocity(Spacetimes e_metric, double r, double theta, double Kepler,
                                     c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double dr_metric[4][4], dr_N, dr_omega;

    double d2r_metric[4][4], d2r_N, d2r_omega;

    get_metric_fist_derivatives(e_metric, dr_metric, &dr_N, &dr_omega, r, theta,
                                Kerr_class, RBH_class, Wormhole_class);

    get_metric_second_derivatives(e_metric, d2r_metric, &d2r_N, &d2r_omega, r, theta,
                                  Kerr_class, RBH_class, Wormhole_class);

    double root = sqrt(dr_metric[0][3] * dr_metric[0][3] - dr_metric[0][0] * dr_metric[3][3]);

    return  -Kepler / dr_metric[3][3] * d2r_metric[3][3] + (-d2r_metric[0][3] + 1.0 / root / 2 * (2 * dr_metric[0][3] * d2r_metric[0][3] - dr_metric[0][0] * d2r_metric[3][3] - d2r_metric[0][0] * dr_metric[3][3])) / dr_metric[3][3];

}

double Redshift(Spacetimes e_metric, Disk_Models Disk_Model, double J, double State_Vector[], double r_obs, double theta_obs,
                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double r_source = State_Vector[e_r];
    double theta_source = State_Vector[e_theta];

    double metric_source[4][4], N_source, omega_source;

    get_metric(e_metric, metric_source, &N_source, &omega_source, r_source, theta_source,
               Kerr_class, RBH_class, Wormhole_class);

    double metric_obs[4][4], N_obs, omega_obs;

    get_metric(e_metric, metric_obs, &N_obs, &omega_obs, r_obs, theta_obs,
               Kerr_class, RBH_class, Wormhole_class);

    double metric_ISCO[4][4]{}, N_ISCO{}, omega_ISCO{};

    get_metric(e_metric, metric_ISCO, &N_ISCO, &omega_ISCO, r_source, theta_source,
               Kerr_class, RBH_class, Wormhole_class);

    double U_obs[4] = { 1.0 / N_obs, 0 ,0 , omega_obs / N_obs };

    if (Disk_Model == Novikov_Thorne) {

        double Kepler = Keplerian_angular_velocity(e_metric, r_source, theta_source,
                                                   Kerr_class, RBH_class, Wormhole_class);

        double gamma = 1 / sqrt(-metric_source[0][0] - 2 * metric_source[0][3] * Kepler - metric_source[3][3] * Kepler * Kepler);

        double U_source[4] = { gamma, 0, 0, gamma * Kepler };

        return  (-U_obs[0] + U_obs[3] * J) / (-U_source[0] + U_source[3] * J);

    }
    else if (Disk_Model == Optically_Thin_Toroidal){

        /*
        Here u_t, u_r, u_theta, u_phi are covariant components 
        */

        double rho = r_source*fabs(sin(State_Vector[e_theta]));
        double ell = pow(rho, 3. / 2) / (1 + rho);

        double u_t{}, u_r{}, u_phi{};

        double inv_metric[4][4];

        invert_metric(inv_metric, metric_source);

        u_t = -1.0 / sqrt(-(inv_metric[0][0] - 2 * inv_metric[0][3] * ell + inv_metric[3][3] * ell * ell));
        u_phi = -u_t * ell;

        /*
        If the disk extends below ISCO, keep u_t and u_phi equal to their ISCO values
        and add a radial component to norm the 4-velocity
        */

        if (r_source < get_ISCO(e_metric, Kerr_class, RBH_class, Wormhole_class)) {

            double inv_metric_ISCO[4][4]{};

            invert_metric(inv_metric_ISCO, metric_ISCO);

            rho = get_ISCO(e_metric, Kerr_class, RBH_class, Wormhole_class) * fabs(sin(State_Vector[e_theta]));
            ell = pow(rho, 3. / 2) / (1 + rho);

            u_t = 1.0 / sqrt(-(inv_metric_ISCO[0][0] - 2 * inv_metric_ISCO[0][3] * ell + inv_metric_ISCO[3][3] * ell * ell));
            u_phi = -u_t * ell;
            u_r = sqrt(metric_source[1][1] * (-1 - inv_metric[0][0] * u_t * u_t - inv_metric[3][3] * u_phi * u_phi));

        }

        /*
        U_source is given in contravariant components
        */
        
        double U_source[4] = { inv_metric[0][0] * u_t + inv_metric[0][3] * u_phi,
                               inv_metric[1][1] * u_r                           ,
                               0                                                ,
                               inv_metric[3][3] * u_phi + inv_metric[3][0] * u_t };

        return  (-U_obs[0] + U_obs[3] * J) / (-U_source[0] + U_source[1] * State_Vector[e_p_r] + U_source[3] * J);

    }

}

double Flux_integrand(Spacetimes e_metric, double r, double theta, double* metric_det, double* E_disk, double* L_disk,
                      c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double metric[4][4], N, omega;

    double dr_metric[4][4], dr_N, dr_omega;

    get_metric(e_metric, metric, &N, &omega, r, theta,
               Kerr_class, RBH_class, Wormhole_class);

    get_metric_fist_derivatives(e_metric, dr_metric, &dr_N, &dr_omega, r, theta,
                                Kerr_class, RBH_class, Wormhole_class);

    double Kepler    = Keplerian_angular_velocity(e_metric, r, theta,
                                                  Kerr_class, RBH_class, Wormhole_class);
    double dr_Kepler = dr_Keplerian_angular_velocity(e_metric, r, theta, Kepler,
                                                     Kerr_class, RBH_class, Wormhole_class);

    *metric_det = metric[1][1] * (metric[0][3] * metric[0][3] - metric[0][0] * metric[3][3]);

    double root = sqrt(-metric[0][0] - 2 * metric[0][3] * Kepler - metric[3][3] * Kepler * Kepler);
    double dr_root = (-dr_metric[0][0] - 2 * (dr_metric[0][3] * Kepler + metric[0][3] * dr_Kepler) - dr_metric[3][3] * Kepler * Kepler - 2 * metric[3][3] * Kepler * dr_Kepler);
    *E_disk = -(metric[0][0] + metric[0][3] * Kepler) / root;
    *L_disk = (metric[3][3] * Kepler + metric[0][3]) / root;
    double dr_L = (dr_metric[3][3] * Kepler + metric[3][3] * dr_Kepler + dr_metric[0][3]) / root - *L_disk / root / root / 2 * dr_root;

    return (*E_disk - Kepler * *L_disk) * dr_L;

}

double solve_Flux_integral(Spacetimes e_metric, double lower_bound, double upper_bound, double theta, double tolerance,
                          c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double metric_det, E_disk, L_disk;

    double mid_point = (lower_bound + upper_bound) / 2;
    double left_mid_point = (lower_bound + mid_point) / 2;
    double right_mid_point = (mid_point + upper_bound) / 2;

    double F_lower_bound = Flux_integrand(e_metric, lower_bound, theta, &metric_det, &E_disk, &L_disk,
                                          Kerr_class, RBH_class, Wormhole_class);
    double F_mid_point = Flux_integrand(e_metric, mid_point, theta, &metric_det, &E_disk, &L_disk,
                                        Kerr_class, RBH_class, Wormhole_class);
    double F_upper_bound = Flux_integrand(e_metric, upper_bound, theta, &metric_det, &E_disk, &L_disk,
                                          Kerr_class, RBH_class, Wormhole_class);
 
    double F_left_mid = Flux_integrand(e_metric, left_mid_point, theta, &metric_det, &E_disk, &L_disk,
                                       Kerr_class, RBH_class, Wormhole_class);
    double F_right_mid = Flux_integrand(e_metric, right_mid_point, theta, &metric_det, &E_disk, &L_disk,
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

        double L_value = solve_Flux_integral(e_metric, lower_bound, mid_point, theta, tolerance / 2,
                                             Kerr_class, RBH_class, Wormhole_class);
        double R_value = solve_Flux_integral(e_metric, mid_point, upper_bound, theta, tolerance / 2,
                                             Kerr_class, RBH_class, Wormhole_class);

        integral = L_value + R_value;

    }

    return integral;
}

double get_flux(Spacetimes e_metric, double r, double r_in, double theta,
                c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double metric_det{}, E_disk{}, L_disk{};

    double Kepler = Keplerian_angular_velocity(e_metric, r, theta,
                                               Kerr_class, RBH_class, Wormhole_class);
    double dr_Kepler = dr_Keplerian_angular_velocity(e_metric, r, theta, Kepler,
                                                     Kerr_class, RBH_class, Wormhole_class);

    Flux_integrand(e_metric, r, theta, &metric_det, &E_disk, &L_disk, 
                   Kerr_class, RBH_class, Wormhole_class);

    double Flux_coeff = -dr_Kepler / ((E_disk - Kepler * L_disk) * (E_disk - Kepler * L_disk)) / (4 * M_PI * sqrt(metric_det));

    double r_throat = Wormhole_class.get_r_throat();
        
    if (e_metric == Wormhole) {

        r = sqrt(r * r + r_throat * r_throat);

    }
        
    double Flux_integral = solve_Flux_integral(e_metric, r_in, r, theta, INTEGRAL_ACCURACY,
                                               Kerr_class, RBH_class, Wormhole_class);

    return Flux_coeff * Flux_integral;

}

int get_Radiative_Transfer(double State_Vector[], double Derivatives[], int iteration, double alpha, double Height_Scale , double r_cut, double omega, double redshift) {

    double r = State_Vector[e_r];
    double h = cos(State_Vector[e_theta]);

    double Height_Cutoff = h * h / Height_Scale / Height_Scale / 2;
    double Radial_Cutoff = (r - r_cut) * (r - r_cut) / omega / omega;

    double emission = pow(r, -alpha) * exp(-Height_Cutoff);

    if (State_Vector[e_r] < r_cut) {

        emission *= exp(-Radial_Cutoff);

    }

    double absorbtion = r * redshift * redshift * emission;

    Derivatives[e_Intensity + iteration * e_State_Number]     = -redshift * redshift * emission*exp(-State_Vector[e_Optical_Depth]);
    Derivatives[e_Optical_Depth + iteration * e_State_Number] = -absorbtion / redshift;

    return  OK;

}

int get_EOM(Spacetimes e_metric, double inter_State_vector[7 * 6], double J, double Derivatives[7 * 6], int iteration,
            c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    switch (e_metric) {

    case Kerr:

        Kerr_class.EOM(inter_State_vector, J, Derivatives, iteration);

        break;

    case Wormhole:

        Wormhole_class.EOM(inter_State_vector, J, Derivatives, iteration);

        break;

    case Reg_Black_Hole:

        RBH_class.EOM(inter_State_vector, J, Derivatives, iteration);

        break;

    
    default:

        std::cout << "Wrong metric!" << '\n';

        return ERROR;
    }

    return OK;
}

void print_ASCII_art() {

    std::cout <<

        " ######   ########     ###    ##     ## #### ########    ###    ######## ####  #######  ##    ##    ###    ##          ########     ###    ##    ##    ######## ########     ###     ######  ######## ########  \n"
        "##    ##  ##     ##   ## ##   ##     ##  ##     ##      ## ##      ##     ##  ##     ## ###   ##   ## ##   ##          ##     ##   ## ##    ##  ##        ##    ##     ##   ## ##   ##    ## ##       ##     ## \n"
        "##        ##     ##  ##   ##  ##     ##  ##     ##     ##   ##     ##     ##  ##     ## ####  ##  ##   ##  ##          ##     ##  ##   ##    ####         ##    ##     ##  ##   ##  ##       ##       ##     ## \n"
        "##   #### ########  ##     ## ##     ##  ##     ##    ##     ##    ##     ##  ##     ## ## ## ## ##     ## ##          ########  ##     ##    ##          ##    ########  ##     ## ##       ######   ########  \n"
        "##    ##  ##   ##   #########  ##   ##   ##     ##    #########    ##     ##  ##     ## ##  #### ######### ##          ##   ##   #########    ##          ##    ##   ##   ######### ##       ##       ##   ##   \n"
        "##    ##  ##    ##  ##     ##   ## ##    ##     ##    ##     ##    ##     ##  ##     ## ##   ### ##     ## ##          ##    ##  ##     ##    ##          ##    ##    ##  ##     ## ##    ## ##       ##    ##  \n"
        " ######   ##     ## ##     ##    ###    ####    ##    ##     ##    ##    ####  #######  ##    ## ##     ## ########    ##     ## ##     ##    ##          ##    ##     ## ##     ##  ######  ######## ##     ## \n";

    std::cout << '\n' << '\n';

}

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

    for (int i = 0; i <= max_digits + current_digits; i += 1) {

        std::cout << "\b";

    }

    std::cout << current + 1 << "/" << max + 1;

}

Disk_Intersection Disk_event(Disk_Models e_Disk_Model, double State_Vector[], double Old_State_Vector[], double r_in, double r_out,
                             double disk_alpha, double disk_height_scale, double disk_rad_cutoff, double disk_omega) {

    bool inside_disk{};
    bool is_inside{};
    bool was_inside{};

    switch (e_Disk_Model) {

    case(Novikov_Thorne):

        inside_disk = State_Vector[e_r] * State_Vector[e_r] > r_in * r_in &&
                      State_Vector[e_r] * State_Vector[e_r] < r_out * r_out &&
                      crossed_equatior(State_Vector, Old_State_Vector);

        if (inside_disk) {

            return Inside_Disk;
        }

        return Outside_Disk;

    /*

    case(Optically_Thin_Toroidal):

        is_inside = get_Radiative_Transfer(State_Vector, disk_alpha, disk_height_scale, disk_rad_cutoff, disk_omega, 1) > TOROIDAL_DISK_BOUNDY_DENSITY;
        was_inside = get_Radiative_Transfer(Old_State_Vector, disk_alpha, disk_height_scale, disk_rad_cutoff, disk_omega, 1) > TOROIDAL_DISK_BOUNDY_DENSITY;

        if (is_inside == true) {

            return Inside_Disk;

        }
        else if (is_inside == false && was_inside == true) {

            return Exiting_Disk;

        }

        return Outside_Disk;
        */
    }

    
}