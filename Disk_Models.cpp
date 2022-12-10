#pragma once

#include "Enumerations.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "Disk_Models.h"
#include "General_GR_functions.h"
#include "General_math_functions.h"

#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>

/***************************************************
|                                                  |
| Novikov-Thorne Model Class Functions Definitions |
|                                                  |
***************************************************/

extern e_Spacetimes e_metric;

tag_Novikov_Thorne_Model::tag_Novikov_Thorne_Model(double x, double y) {

    r_in = x;
    r_out = y;

};

double tag_Novikov_Thorne_Model::Keplerian_angular_velocity(double r, std::vector<c_Spacetime_Base*> Spacetimes) {

    double dr_metric[4][4], dr_N, dr_omega;

    Spacetimes[e_metric]->get_dr_metric(dr_metric, &dr_N, &dr_omega, r, M_PI_2);

    return (-dr_metric[0][3] + sqrt(dr_metric[0][3] * dr_metric[0][3] - dr_metric[0][0] * dr_metric[3][3])) / dr_metric[3][3];

}

double tag_Novikov_Thorne_Model::dr_Keplerian_angular_velocity(double r, std::vector<c_Spacetime_Base*> Spacetimes) {

    double dr_metric[4][4], dr_N, dr_omega;

    double d2r_metric[4][4], d2r_N, d2r_omega;

    Spacetimes[e_metric]->get_dr_metric(dr_metric, &dr_N, &dr_omega, r, M_PI_2);

    Spacetimes[e_metric]->get_d2r_metric(d2r_metric, &d2r_N, &d2r_omega, r, M_PI_2);

    double root = sqrt(dr_metric[0][3] * dr_metric[0][3] - dr_metric[0][0] * dr_metric[3][3]);

    double Kepler = tag_Novikov_Thorne_Model::Keplerian_angular_velocity(r, Spacetimes);

    return  -Kepler / dr_metric[3][3] * d2r_metric[3][3] + (-d2r_metric[0][3] + 1.0 / root / 2 * (2 * dr_metric[0][3] * d2r_metric[0][3] - dr_metric[0][0] * d2r_metric[3][3] - d2r_metric[0][0] * dr_metric[3][3])) / dr_metric[3][3];

}

double tag_Novikov_Thorne_Model::Redshift(double J, double State_Vector[], double r_obs, double theta_obs,
                                          std::vector<c_Spacetime_Base*> Spacetimes) {

    double r_source = State_Vector[e_r];
    double theta_source = State_Vector[e_theta];

    /*
    Get the observer 4-velocity
    */

    double metric_obs[4][4], N_obs, omega_obs;

    Spacetimes[e_metric]->get_metric(metric_obs, &N_obs, &omega_obs, r_obs, theta_obs);

    double U_obs[4] = { 1.0 / N_obs, 0 ,0 , omega_obs / N_obs };

    /*
    Get the source 4-velocity
    */

    double metric_source[4][4], N_source, omega_source;

    Spacetimes[e_metric]->get_metric(metric_source, &N_source, &omega_source, r_source, M_PI_2);

    double Kepler = Keplerian_angular_velocity(r_source, Spacetimes);

    double gamma = 1 / sqrt(-metric_source[0][0] - 2 * metric_source[0][3] * Kepler - metric_source[3][3] * Kepler * Kepler);

    double U_source[4] = { gamma, 0, 0, gamma * Kepler };

    return  (-U_obs[0] + U_obs[3] * J) / (-U_source[0] + U_source[3] * J);

}

double tag_Novikov_Thorne_Model::disk_Energy(double r, std::vector<c_Spacetime_Base*> Spacetimes) {


    double metric[4][4], N, omega;

    Spacetimes[e_metric]->get_metric(metric, &N, &omega, r, M_PI_2);

    double Kepler = Keplerian_angular_velocity(r, Spacetimes);

    double root = sqrt(-metric[0][0] - 2 * metric[0][3] * Kepler - metric[3][3] * Kepler * Kepler);

    return  -(metric[0][0] + metric[0][3] * Kepler) / root;

}

double tag_Novikov_Thorne_Model::disk_Angular_Momentum(double r, std::vector<c_Spacetime_Base*> Spacetimes) {

    double metric[4][4], N, omega;

    Spacetimes[e_metric]->get_metric(metric, &N, &omega, r, M_PI_2);

    double Kepler = Keplerian_angular_velocity(r, Spacetimes);

    double root = sqrt(-metric[0][0] - 2 * metric[0][3] * Kepler - metric[3][3] * Kepler * Kepler);

    return  (metric[3][3] * Kepler + metric[0][3]) / root;

}

double tag_Novikov_Thorne_Model::Flux_integrand(double r, std::vector<c_Spacetime_Base*> Spacetimes) {

    double metric[4][4], N, omega;

    double dr_metric[4][4], dr_N, dr_omega;

    Spacetimes[e_metric]->get_metric(metric, &N, &omega, r, M_PI_2);

    Spacetimes[e_metric]->get_dr_metric(dr_metric, &dr_N, &dr_omega, r, M_PI_2);

    double Kepler = Keplerian_angular_velocity(r, Spacetimes);
    double dr_Kepler = dr_Keplerian_angular_velocity(r, Spacetimes);

    double metric_det = get_metric_det(metric);

    double root = sqrt(-metric[0][0] - 2 * metric[0][3] * Kepler - metric[3][3] * Kepler * Kepler);
    double dr_root = (-dr_metric[0][0] - 2 * (dr_metric[0][3] * Kepler + metric[0][3] * dr_Kepler) - dr_metric[3][3] * Kepler * Kepler - 2 * metric[3][3] * Kepler * dr_Kepler);

    double E = disk_Energy(r, Spacetimes);
    double L = disk_Angular_Momentum(r, Spacetimes);

    double dr_L = (dr_metric[3][3] * Kepler + metric[3][3] * dr_Kepler + dr_metric[0][3]) / root - L / root / root / 2 * dr_root;

    return (E - Kepler * L) * dr_L;

}

double tag_Novikov_Thorne_Model::solve_Flux_integral(double lower_bound, double upper_bound, double tolerance, std::vector<c_Spacetime_Base*> Spacetimes) {

    double mid_point       = (lower_bound + upper_bound) / 2;
    double left_mid_point  = (lower_bound + mid_point) / 2;
    double right_mid_point = (mid_point + upper_bound) / 2;

    double F_lower_bound = Flux_integrand(lower_bound, Spacetimes);
    double F_mid_point   = Flux_integrand(mid_point, Spacetimes);
    double F_upper_bound = Flux_integrand(upper_bound, Spacetimes);

    double F_left_mid = Flux_integrand(left_mid_point, Spacetimes);
    double F_right_mid = Flux_integrand(right_mid_point, Spacetimes);

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

        double L_value = solve_Flux_integral(lower_bound, mid_point, tolerance / 2, Spacetimes);
        double R_value = solve_Flux_integral(mid_point, upper_bound, tolerance / 2, Spacetimes);

        integral = L_value + R_value;

    }

    return integral;
}

double tag_Novikov_Thorne_Model::get_flux(double r, std::vector<c_Spacetime_Base*> Spacetimes) {

    double metric[4][4], N, omega;

    Spacetimes[e_metric]->get_metric(metric, &N, &omega, r, M_PI_2);

    double metric_det = get_metric_det(metric);
    double E_disk = disk_Energy(r, Spacetimes);
    double L_disk = disk_Angular_Momentum(r, Spacetimes);

    double Kepler    =    Keplerian_angular_velocity(r, Spacetimes);
    double dr_Kepler = dr_Keplerian_angular_velocity(r, Spacetimes);

    double Flux_coeff = -dr_Kepler / ((E_disk - Kepler * L_disk) * (E_disk - Kepler * L_disk)) / (4 * M_PI * sqrt(metric_det));

    if (e_metric == Wormhole) {

        r = sqrt(r * r + WH_R_THROAT * WH_R_THROAT);

    }

    double Flux_integral = solve_Flux_integral(tag_Novikov_Thorne_Model::r_in, r, INTEGRAL_ACCURACY, Spacetimes);

    return Flux_coeff * Flux_integral;

}

double tag_Novikov_Thorne_Model::get_r_in()  { return r_in; };
double tag_Novikov_Thorne_Model::get_r_out() { return r_out; };

/************************************************************
|                                                           |
| Optically Thin Toroidal Model Class Functions Definitions |
|                                                           |
************************************************************/

tag_Optically_Thin_Toroidal_Model::tag_Optically_Thin_Toroidal_Model(Real alpha, Real height_scale, Real rad_cutoff, Real cutoff_scale,
                                                                     Real magnetization, Real mag_field[3]) {

    DISK_ALPHA = alpha;
    DISK_HEIGHT_SCALE = height_scale;
    DISK_RAD_CUTOFF = rad_cutoff;
    DISK_CUTOFF_SCALE = cutoff_scale;
    DISK_MAGNETIZATION = magnetization;
    MAG_FIELD_GEOMETRY[0] = mag_field[0];
    MAG_FIELD_GEOMETRY[1] = mag_field[1];
    MAG_FIELD_GEOMETRY[2] = mag_field[2];

}

double tag_Optically_Thin_Toroidal_Model::get_disk_alpha()         { return DISK_ALPHA; };
double tag_Optically_Thin_Toroidal_Model::get_disk_height_scale()  { return DISK_HEIGHT_SCALE; };
double tag_Optically_Thin_Toroidal_Model::get_disk_rad_cutoff()    { return DISK_RAD_CUTOFF; };;
double tag_Optically_Thin_Toroidal_Model::get_disk_omega()         { return DISK_CUTOFF_SCALE; };
double tag_Optically_Thin_Toroidal_Model::get_disk_magnetization() { return DISK_MAGNETIZATION; };

double tag_Optically_Thin_Toroidal_Model::get_disk_temperature(double State_vector[]) {

    double r = State_vector[e_r];

    double T = T_ELECTRON_CGS * (1 + sqrt(1 - SPIN * SPIN)) / r;

    if (e_metric != Kerr) {

        r = sqrt(r * r + WH_R_THROAT * WH_R_THROAT);

        T = T_ELECTRON_CGS / r;

    }

    return T;

}

int tag_Optically_Thin_Toroidal_Model::get_disk_velocity(double Disk_velocity[], double State_Vector[], std::vector<c_Spacetime_Base*> Spacetimes) {


    double r_source = State_Vector[e_r];

    if (e_metric == Wormhole) {

        r_source = sqrt(State_Vector[e_r] * State_Vector[e_r] + WH_R_THROAT * WH_R_THROAT);

    }

    double& theta_source = State_Vector[e_theta];

    double metric_source[4][4], N_source, omega_source;

    Spacetimes[e_metric]->get_metric(metric_source, &N_source, &omega_source, r_source, theta_source);

    double rho = r_source * sin(theta_source);

    if (rho < 0.0) {

        rho *= -1.0;

    }

    double ell = pow(sqrt(rho), 3) / (1 + rho);

    double u_t{}, u_phi{};

    double inv_metric[4][4]{};

    invert_metric(inv_metric, metric_source);

    u_t = -1.0 / sqrt(-(inv_metric[0][0] - 2 * inv_metric[0][3] * ell + inv_metric[3][3] * ell * ell));
    u_phi = -u_t * ell;

    /*
    Convert U_source to contravariant components
    */

    Disk_velocity[0] = inv_metric[0][0] * u_t + inv_metric[0][3] * u_phi;
    Disk_velocity[1] = 0;
    Disk_velocity[2] = 0;
    Disk_velocity[3] = inv_metric[3][3] * u_phi + inv_metric[3][0] * u_t;

    return OK;

}

double tag_Optically_Thin_Toroidal_Model::get_disk_density(double State_Vector[]) {

    double& r = State_Vector[e_r];
    double rho = sin(State_Vector[e_theta]);
    double h = cos(State_Vector[e_theta]);

    double Height_Cutoff = h * h / (2 * (DISK_ALPHA * rho) * (DISK_ALPHA * rho));
    double Radial_Cutoff = (r - DISK_RAD_CUTOFF) * (r - DISK_RAD_CUTOFF) * DISK_CUTOFF_SCALE;

    double electron_density = N_ELECTRON_CGS * pow(r / (1. + sqrt(1 - SPIN * SPIN)), -2) * exp(-Height_Cutoff);

    //double electron_density = exp(-((r / 10) * (r / 10) + 100. / 3 * 100. / 3 * h * h) / 2);

    //if (State_Vector[e_r] < DISK_RAD_CUTOFF) {

    //    electron_density *= exp(-Radial_Cutoff);

    //}

    return electron_density;

}

double tag_Optically_Thin_Toroidal_Model::get_magnetic_field(double B_field[3], double State_vector[]) {

    /*

    Everything is in GCS!

    */

    double electron_density = get_disk_density(State_vector);

    double B_CGS = sqrt(DISK_MAGNETIZATION * C_LIGHT_CGS * C_LIGHT_CGS * electron_density * M_PROTON_CGS * 4 * M_PI);

    B_field[x] = B_CGS * MAG_FIELD_GEOMETRY[x];
    B_field[y] = B_CGS * MAG_FIELD_GEOMETRY[y];
    B_field[z] = B_CGS * MAG_FIELD_GEOMETRY[z];

    return B_CGS;

}

double tag_Optically_Thin_Toroidal_Model::get_electrron_pitch_angle(double State_vector[], double B_field_local[], std::vector<c_Spacetime_Base*> Spacetimes) {

    double U_source_coord[4];
    get_disk_velocity(U_source_coord, State_vector, Spacetimes);

    /*

    Transform U_source To The ZAMO Frame

    */

    double metric[4][4]{}, N_metric{}, Omega_metric{};

    Spacetimes[e_metric]->get_metric(metric, &N_metric, &Omega_metric, State_vector[e_r], State_vector[e_theta]);

    double U_source_ZAMO[4]{};
    Contravariant_coord_to_ZAMO(metric, U_source_coord, U_source_ZAMO);

    /*

    Boost U_source_ZAMO To The Fluid Frame

    */

    double Boost_matrix[4][4];

    Lorentz_boost_matrix(Boost_matrix, U_source_ZAMO, metric);

    double* U_source_Boosted  = mat_vec_multiply_4D(Boost_matrix, U_source_ZAMO);
    double  U_source_local[3] = { U_source_Boosted[1] / U_source_Boosted[0], U_source_Boosted[2] / U_source_Boosted[0], U_source_Boosted[3] / U_source_Boosted[0] };

    /*

    Get Sin Of The Angle Between The Local Fluid Velocity And Magnetic Field In The ZAMO Frame

    */

    double cos_angle = 1, sin_angle{};

    if (vector_norm(B_field_local, 3) > 1e-10 && dot_product(U_source_local, U_source_local) != 0) {

        cos_angle = dot_product(B_field_local, U_source_local) / sqrt(dot_product(B_field_local, B_field_local) * dot_product(U_source_local, U_source_local));
        sin_angle = sqrt(1 - cos_angle * cos_angle);

    }

    return sin_angle;

}

double tag_Optically_Thin_Toroidal_Model::get_emission_fucntion(double State_vector[], double J, std::vector<c_Spacetime_Base*> Spacetimes) {

    /* Electron Density in CGS */

    double electron_density = get_disk_density(State_vector);

    /* Dimentionless Electron Temperature */

    double T_electron       = get_disk_temperature(State_vector);
    double T_electron_dim   = BOLTZMANN_CONST_CGS * T_electron / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;

    /* Magnetic Field */

    double B_field_local[3];
    double B_CGS = get_magnetic_field(B_field_local, State_vector);

    /* Disk Coordinate Velocity */

    double U_source_coord[4]{};
    get_disk_velocity(U_source_coord, State_vector, Spacetimes);

    /* Synchotron Frequency */

    double f_cyclo         = Q_ELECTRON_CGS * B_CGS / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);
    double sin_pitch_angle = get_electrron_pitch_angle(State_vector, B_field_local, Spacetimes);

    double f_s = 2. / 9 * f_cyclo * T_electron_dim * T_electron_dim * sin_pitch_angle;

    /* Dimentionless Redshifted Synchotron Frequency */

    double redshift = Redshift(J, State_vector, U_source_coord);

    double X = 1e100;

    if (f_s != 0) {

        X = OBS_FREQUENCY_CGS / f_s / redshift;

    }

    /* Emission Function */

    double X_term = (sqrt(X) + pow(2, 11.0 / 12) * pow(X, 1.0 / 6)) * (sqrt(X) + pow(2, 11.0 / 12) * pow(X, 1.0 / 6));
    double constant_coeff = sqrt(2) * M_PI * Q_ELECTRON_CGS * Q_ELECTRON_CGS / 3 / C_LIGHT_CGS;

    return constant_coeff * electron_density * f_s * X_term * exp(-pow(X, 1.0 / 3)) / std::cyl_bessel_k(2.0, 1.0 / T_electron_dim);

}

double tag_Optically_Thin_Toroidal_Model::get_absorbtion_fucntion(double Emission_Function, double Frequency, double Temperature) {

    double Planck_function_CGS = get_planck_function_CGS(Frequency, Temperature);

    return Emission_Function/Planck_function_CGS;

}

