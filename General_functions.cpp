#pragma once

#include "Constants.h"
#include "Enumerations.h"
#include "Spacetimes.h"
#include "Disk_Models.h"

#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>

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

double get_ISCO(e_Spacetimes e_metric, c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    switch (e_metric) {

    case Kerr:

        return Kerr_class.get_ISCO();

    case Wormhole:

        return Wormhole_class.get_ISCO();

    case Reg_Black_Hole:

        return RBH_class.get_ISCO();

    }

}

int get_metric(e_Spacetimes e_metric, double metric[4][4], double* N_metric, double* omega, double r, double theta,
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

double get_metric_det(double metric[4][4]) {

    return metric[1][1] * (metric[0][3] * metric[0][3] - metric[0][0] * metric[3][3]);

}

int get_metric_fist_derivatives(e_Spacetimes e_metric, double dr_metric[4][4], double* dr_N_metric, double* dr_omega,
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

int get_metric_second_derivatives(e_Spacetimes e_metric, double d2r_metric[4][4], double* d2r_N_metric, double* d2r_omega,
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

int get_initial_conditions_from_file(e_Spacetimes e_metric, double* J, double J_data[], double* p_theta, double p_theta_data[],
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

double Redshift(e_Spacetimes e_metric, Disk_Models Disk_Model, double J, double State_Vector[], double r_obs, double theta_obs,
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

    if (Disk_Model == Optically_Thin_Toroidal) {

        /*
        Here u_t, u_r, u_theta, u_phi are covariant components
        */

        double rho = r_source * fabs(sin(State_Vector[e_theta]));
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

    return ERROR;

}


int Loretnz_boost() {






    return OK;
}

int get_Radiative_Transfer(double State_Vector[], double Derivatives[], int iteration, Optically_Thin_Toroidal_Model OTT_Model, double redshift) {

    double r = State_Vector[e_r];
    double h = cos(State_Vector[e_theta]);

    double Height_Scale = OTT_Model.get_disk_height_scale();
    double r_cut        = OTT_Model.get_disk_rad_cutoff();
    double omega        = OTT_Model.get_disk_omega();
    double alpha        = OTT_Model.get_disk_alpha();

    double Height_Cutoff = h * h / Height_Scale / Height_Scale / 2;
    double Radial_Cutoff = (r - r_cut) * (r - r_cut) / omega / omega;

    double electron_density =  pow(r, -alpha) * exp(-Height_Cutoff);

    if (State_Vector[e_r] < r_cut) {

        electron_density *= exp(-Radial_Cutoff);

    }

    double T_electron_cgs = T_ELECTRON_CGS * 2 / r;
    double T_electron_dim = BOLTZMANN_CONST_CGS * T_electron_cgs / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;

    double B_CGS = sqrt(0.01 * C_LIGHT_CGS * C_LIGHT_CGS * electron_density * M_PROTON_CGS * 4 * M_PI);

    double f_cyclo = Q_ELECTRON_CGS * B_CGS / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);
    double f_s = 2. / 9 * f_cyclo * T_electron_dim * T_electron_dim; // Lorentz boost needed to calculate angle between photon and B-field

    if (r < 6) {

        int test = 1;

    }

    double absorbtion = r * redshift * redshift * electron_density; 
    Derivatives[e_Intensity + iteration * e_State_Number] = -redshift * redshift * electron_density * exp(-State_Vector[e_Optical_Depth]);
    Derivatives[e_Optical_Depth + iteration * e_State_Number] = -absorbtion / redshift;

    return  OK;

}

int get_EOM(e_Spacetimes e_metric, double inter_State_vector[7 * 6], double J, double Derivatives[7 * 6], int iteration,
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

Disk_Intersection Disk_event(Disk_Models e_Disk_Model, double State_Vector[], double Old_State_Vector[],
    Novikov_Thorne_Model NT_Model, Optically_Thin_Toroidal_Model OTT_Model) {

    bool inside_disk{};

    double r_in = NT_Model.get_r_in();
    double r_out = NT_Model.get_r_out();

    switch (e_Disk_Model) {

    case(Novikov_Thorne):

        inside_disk = State_Vector[e_r] * State_Vector[e_r] > r_in * r_in &&
            State_Vector[e_r] * State_Vector[e_r] < r_out * r_out &&
            crossed_equatior(State_Vector, Old_State_Vector);

        if (inside_disk) {

            return Inside_Disk;
        }

        return Outside_Disk;

    }

}