#pragma once

#include "Enumerations.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "Disk_Models.h"

#include <iostream>
#include <cmath>
#include <vector>

extern e_Spacetimes e_metric;
extern std::vector<c_Spacetime_Base*> Spacetimes;
extern c_Observer Observer_class;
extern Optically_Thin_Toroidal_Model OTT_Model;
extern Novikov_Thorne_Model NT_Model;

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

double get_metric_det(double metric[4][4]) {

    return metric[1][1] * (metric[0][3] * metric[0][3] - metric[0][0] * metric[3][3]);

}

int get_intitial_conditions_from_angles(Initial_conditions_type* p_Initial_Conditions, double V_angle, double H_angle) {

    double g2, gamma, ksi, L_z, E;

    double (*metric)[4] = p_Initial_Conditions->init_metric;

    g2 = pow(metric[0][3], 2) - metric[0][0] * metric[3][3];
    ksi = sqrt(metric[3][3] / g2);
    gamma = -metric[0][3] / metric[3][3] * ksi;

    L_z = sqrt(metric[3][3]) * sin(H_angle + 2 * M_PI) * cos(V_angle);
    E = (1 + gamma * L_z) / ksi;

    p_Initial_Conditions->init_Three_Momentum[e_phi] = L_z / E;
    p_Initial_Conditions->init_Three_Momentum[e_theta] = sqrt(metric[2][2]) * sin(V_angle) / E;
    p_Initial_Conditions->init_Three_Momentum[e_r] = sqrt(metric[1][1]) * cos(H_angle + 2 * M_PI) * cos(V_angle) / E;

    return OK;
}

double Redshift(double J, double State_Vector[], double U_source[]) {

    double U_obs[4];

    Observer_class.get_obs_velocity(U_obs, Spacetimes);

    return  (-U_obs[0] + U_obs[3] * J) / (-U_source[0] * 1 + 
                                           U_source[1] * State_Vector[e_p_r] + 
                                           U_source[2] * State_Vector[e_p_theta] + 
                                           U_source[3] * J);

}

double get_photon_t_momentum(double State_vector[], double J, double metric[4][4]) {

    /*
    
    State_vector hold covariant momentum components
    
    */

    double p_r = State_vector[e_p_r];
    double p_theta = State_vector[e_p_theta];
    double p_phi = J;
    
    double inv_metric[4][4];

    invert_metric(inv_metric, metric);

    double inv_g2 = inv_metric[e_t_coord][e_phi_coord] * inv_metric[e_t_coord][e_phi_coord] 
                  - inv_metric[e_phi_coord][e_phi_coord] * inv_metric[e_t_coord][e_t_coord];

    double root = p_phi * p_phi * inv_g2 / inv_metric[e_t_coord][e_t_coord] / inv_metric[e_t_coord][e_t_coord]
                + inv_metric[e_theta_coord][e_theta_coord] / inv_metric[e_t_coord][e_t_coord] * p_theta * p_theta
                + inv_metric[e_r_coord][e_r_coord] / inv_metric[e_t_coord][e_t_coord] * p_r * p_r;

    double p_t = -inv_metric[0][3] / inv_metric[0][0] * p_phi - sqrt(root);

    /*
    
    Returns covatiant component
    
    */

    return p_t;

}

int get_ZAMO_tetrad(double e_t[4], double e_r[4], double e_theta[4], double e_phi[4], double metric[4][4]) {

    double ksi = 1.0 / sqrt(-metric[0][0] + metric[0][3] * metric[0][3] / metric[3][3]);
    double sigma = -metric[0][3] / metric[3][3] * ksi;

    e_t[e_t_coord] = ksi;
    e_t[e_r_coord] = 0;
    e_t[e_theta_coord] = 0;
    e_t[e_phi_coord] = sigma;

    e_r[e_t_coord] = 0;
    e_r[e_r_coord] = 1.0 / sqrt(metric[1][1]);
    e_r[e_theta_coord] = 0;
    e_r[e_phi_coord] = 0;

    e_theta[e_t_coord] = 0;
    e_theta[e_r_coord] = 0;
    e_theta[e_theta_coord] = 1.0 / sqrt(metric[2][2]);
    e_theta[e_phi_coord] = 0;

    e_phi[e_t_coord] = 0;
    e_phi[e_r_coord] = 0;
    e_phi[e_theta_coord] = 0;
    e_phi[e_phi_coord] = 1.0 / sqrt(metric[3][3]);

    return OK;

}

int Contravariant_coord_to_ZAMO(double metric[4][4], double Contravariant_Vector[4], double ZAMO_Vector[4]) {

    double alpha = sqrt(-metric[0][0] + metric[0][3] * metric[0][3] / metric[3][3]);
    double beta = metric[0][3] / metric[3][3];

    ZAMO_Vector[0] = alpha * Contravariant_Vector[0];
    ZAMO_Vector[1] = sqrt(metric[1][1]) * Contravariant_Vector[1];
    ZAMO_Vector[2] = sqrt(metric[2][2]) * Contravariant_Vector[2];
    ZAMO_Vector[3] = sqrt(metric[3][3]) * (Contravariant_Vector[3] + beta * Contravariant_Vector[0]);

    return OK;

}

int Lorentz_boost_matrix(double Boost_matrix[4][4], double U_source[4], double metric[4][4]) {

    double V_r_ZAMO     = U_source[1]/U_source[0];
    double V_theta_ZAMO = U_source[2]/U_source[0];
    double V_phi_ZAMO   = U_source[3]/U_source[0];

    double V_squared = V_r_ZAMO * V_r_ZAMO + V_theta_ZAMO * V_theta_ZAMO + V_phi_ZAMO * V_phi_ZAMO;

    double gamma = 1.0 / sqrt(1 - V_squared);

    Boost_matrix[0][0] = gamma;
    Boost_matrix[0][1] = gamma * V_r_ZAMO;
    Boost_matrix[0][2] = gamma * V_theta_ZAMO;
    Boost_matrix[0][3] = gamma * V_phi_ZAMO;

    for (int index = 1; index <= 3; index += 1) {

        Boost_matrix[index][0] = Boost_matrix[0][index];

    }

    Boost_matrix[1][1] = 1 + (gamma - 1) * V_r_ZAMO * V_r_ZAMO / V_squared;
    Boost_matrix[1][2] = (gamma - 1) * V_r_ZAMO * V_theta_ZAMO / V_squared;
    Boost_matrix[1][3] = (gamma - 1) * V_r_ZAMO * V_phi_ZAMO / V_squared;

    Boost_matrix[2][1] = Boost_matrix[1][2];
    Boost_matrix[3][1] = Boost_matrix[1][3];

    Boost_matrix[2][2] = 1 + (gamma - 1) * V_theta_ZAMO * V_theta_ZAMO / V_squared;
    Boost_matrix[2][3] = (gamma - 1) * V_theta_ZAMO * V_phi_ZAMO / V_squared;

    Boost_matrix[3][2] = Boost_matrix[2][3];

    Boost_matrix[3][3] = 1 + (gamma - 1) * V_phi_ZAMO * V_phi_ZAMO / V_squared;


    return OK;
}

int get_Radiative_Transfer(double State_Vector[], double Derivatives[], int iteration, double J) {

    double r;

    switch (e_metric) {

        case Wormhole:

            r = sqrt(State_Vector[e_r] * State_Vector[e_r] + WH_R_THROAT * WH_R_THROAT);

            break;

        default:
         
            r = State_Vector[e_r];

            break;
    }


    
    double electron_density = OTT_Model.get_disk_density(State_Vector);

    /*
    
    Get Disk Cooridinate Velocity
    
    */

    double U_source_coord[4]{};

    OTT_Model.get_disk_velocity(U_source_coord, State_Vector, Spacetimes);

    /*
    
    Transform U_source To The ZAMO Frame
    
    */

    double metric[4][4]{}, N_metric{}, Omega_metric{};

    Spacetimes[e_metric]->get_metric(metric, &N_metric, &Omega_metric, r, State_Vector[e_theta]);

    double U_source_ZAMO[4]{};
    Contravariant_coord_to_ZAMO(metric, U_source_coord, U_source_ZAMO);

    /*

    Boost U_source_ZAMO To The Fluid Frame

    */

    double Boost_matrix[4][4];

    Lorentz_boost_matrix(Boost_matrix, U_source_ZAMO, metric);

    double U_source_Boosted[4]{};

    for (int row = 0; row <= 3; row += 1) {

        for (int column = 0; column <= 3; column += 1) {

            U_source_Boosted[row] += Boost_matrix[row][column] * U_source_ZAMO[column];

        }
    }

    /*
    
    Get Sin Of The Angle Between The Local Fluid Velocity And Magnetic Field In The ZAMO Frame
    
    */

    double U_source_local[3] = { U_source_Boosted[1] / U_source_Boosted[0], U_source_Boosted[2] / U_source_Boosted[0], U_source_Boosted[3] / U_source_Boosted[0] };

    double B_field_local[3];

    double B_CGS = OTT_Model.get_magnetic_field(B_field_local, State_Vector);

    double cos_angle = 1, sin_angle{};

    if (B_CGS != 0 && dot_product(U_source_local, U_source_local) != 0) {

        cos_angle = dot_product(B_field_local, U_source_local) / sqrt(dot_product(B_field_local, B_field_local)) / sqrt(dot_product(U_source_local, U_source_local));
        sin_angle = sqrt(1 - cos_angle * cos_angle);

    }

    /*
    
    Get The Synchotron Frequency And Electron Temperature
    
    */

    double T_electron_cgs = T_ELECTRON_CGS * (1 + sqrt(1 - 0.94 * 0.94)) / r;
    double T_electron_dim = BOLTZMANN_CONST_CGS * T_electron_cgs / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;

    double f_cyclo = Q_ELECTRON_CGS * B_CGS / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);
    double f_s = 2. / 9 * f_cyclo * T_electron_dim * T_electron_dim * sin_angle;

    /*

    Get The Redshift

    */

    double redshift = Redshift(J, State_Vector, U_source_coord);

    /*
   
    Calculate The Emission Function
    
    */

    double X{};

    if (f_s != 0) {

        X = OBS_FREQUENCY_CGS / f_s;

    }

    double X_term = (sqrt(X) + pow(2, 11.0 / 12) * pow(X, 1.0 / 6)) * (sqrt(X) + pow(2, 11.0 / 12) * pow(X, 1.0 / 6));

    double emission_coef = sqrt(2) * M_PI * Q_ELECTRON_CGS * Q_ELECTRON_CGS / 27 / C_LIGHT_CGS;
    double emission = emission_coef * electron_density * f_cyclo * sin_angle * X_term * exp(-pow(X, 1.0 / 3));

    /*

    Calculate the absorbtion

    */

    double absorbtion_coeff = C_LIGHT_CGS * C_LIGHT_CGS / 2 / BOLTZMANN_CONST_CGS / OBS_FREQUENCY_CGS / OBS_FREQUENCY_CGS;
    double absorbtion = absorbtion_coeff * emission / T_electron_cgs;

    /*

    Fill in radiative transfer derivatives

    */

    Derivatives[e_Intensity     + iteration * e_State_Number] = -redshift * redshift * emission * exp(-State_Vector[e_Optical_Depth]) * OBS_FREQUENCY_CGS * 100;
    Derivatives[e_Optical_Depth + iteration * e_State_Number] = -absorbtion / redshift * OBS_FREQUENCY_CGS * 100;

    return  OK;

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

Disk_Intersection Disk_event(double State_Vector[], double Old_State_Vector[]) {

    /*
    
    This function is only relevant for Novikov-Thorne disks
    
    */

    bool inside_NT_disk{};

    double r_in = NT_Model.get_r_in();
    double r_out = NT_Model.get_r_out();


    inside_NT_disk = State_Vector[e_r] * State_Vector[e_r] > r_in * r_in &&
                     State_Vector[e_r] * State_Vector[e_r] < r_out * r_out &&
                     crossed_equatior(State_Vector, Old_State_Vector);

    if (inside_NT_disk) {

        return Inside_Disk;

    }

    return Outside_Disk;

}