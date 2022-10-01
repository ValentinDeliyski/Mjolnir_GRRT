#pragma once

#include "Enumerations.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "Disk_Models.h"
#include "General_functions.h"

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

double tag_Novikov_Thorne_Model::Flux_integrand(e_Spacetimes e_metric, double r, std::vector<c_Spacetime_Base*> Spacetimes) {

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

    double F_lower_bound = Flux_integrand(e_metric, lower_bound, Spacetimes);
    double F_mid_point   = Flux_integrand(e_metric, mid_point, Spacetimes);
    double F_upper_bound = Flux_integrand(e_metric, upper_bound, Spacetimes);

    double F_left_mid = Flux_integrand(e_metric, left_mid_point, Spacetimes);
    double F_right_mid = Flux_integrand(e_metric, right_mid_point, Spacetimes);

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

tag_Optically_Thin_Toroidal_Model::tag_Optically_Thin_Toroidal_Model(double alpha, double height_scale, double rad_cutoff, double omega,
                                                                double magnetization, double mag_field[3]) {

    DISK_ALPHA = alpha;
    DISK_HEIGHT_SCALE = height_scale;
    DISK_RAD_CUTOFF = rad_cutoff;
    DISK_OMEGA = omega;
    DISK_MAGNETIZATION = magnetization;
    MAG_FIELD_GEOMETRY[0] = mag_field[0];
    MAG_FIELD_GEOMETRY[1] = mag_field[1];
    MAG_FIELD_GEOMETRY[2] = mag_field[2];

}

double tag_Optically_Thin_Toroidal_Model::get_disk_alpha()         { return DISK_ALPHA; };
double tag_Optically_Thin_Toroidal_Model::get_disk_height_scale()  { return DISK_HEIGHT_SCALE; };
double tag_Optically_Thin_Toroidal_Model::get_disk_rad_cutoff()    { return DISK_RAD_CUTOFF; };;
double tag_Optically_Thin_Toroidal_Model::get_disk_omega()         { return DISK_OMEGA; };
double tag_Optically_Thin_Toroidal_Model::get_disk_magnetization() { return DISK_MAGNETIZATION; };

int tag_Optically_Thin_Toroidal_Model::get_disk_velocity(double Disk_velocity[], double State_Vector[], std::vector<c_Spacetime_Base*> Spacetimes) {

    double r_source = State_Vector[e_r];
    double theta_source = State_Vector[e_theta];

    double metric_source[4][4], N_source, omega_source;

    Spacetimes[e_metric]->get_metric(metric_source, &N_source, &omega_source, r_source, theta_source);
  
    double rho = r_source * sin(theta_source);

    if (rho < 0.0) {

        rho *= -1.0;

    }

    double ell = pow(sqrt(rho), 3) / (1 + rho);

    double u_t{}, u_phi{};

    double inv_metric[4][4];

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

    double r   = State_Vector[e_r];
    double rho = sin(State_Vector[e_theta]);
    double h   = cos(State_Vector[e_theta]);

    double Height_Cutoff = h * h / (2 * (1 * rho) * (1 * rho));
    double Radial_Cutoff = (r - DISK_RAD_CUTOFF) * (r - DISK_RAD_CUTOFF) / DISK_OMEGA / DISK_OMEGA;

    double electron_density = N_ELECTRON_CGS * pow(r / (1 + sqrt(1 - 0.94 * 0.94)), -2) * exp(-Height_Cutoff);

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
