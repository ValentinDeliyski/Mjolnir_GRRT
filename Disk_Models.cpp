#pragma once

#include "Constants.h"
#include "Enumerations.h"
#include "Spacetimes.h"
#include "Disk_Models.h"
#include "General_functions.h"

#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>

Novikov_Thorne_Model::tag_Novikov_Thorne_Model(double x, double y) {

    r_in = x;
    r_out = y;

};

double Novikov_Thorne_Model::Keplerian_angular_velocity(e_Spacetimes e_metric, double r,
                                                        c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double dr_metric[4][4], dr_N, dr_omega;

    if (r < get_ISCO(e_metric, Kerr_class, RBH_class, Wormhole_class)) {

        r = get_ISCO(e_metric, Kerr_class, RBH_class, Wormhole_class);

    }

    get_metric_fist_derivatives(e_metric, dr_metric, &dr_N, &dr_omega, r, M_PI_2,
        Kerr_class, RBH_class, Wormhole_class);

    return (-dr_metric[0][3] + sqrt(dr_metric[0][3] * dr_metric[0][3] - dr_metric[0][0] * dr_metric[3][3])) / dr_metric[3][3];

}

double Novikov_Thorne_Model::dr_Keplerian_angular_velocity(e_Spacetimes e_metric, double r, double Kepler,
                                                           c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double dr_metric[4][4], dr_N, dr_omega;

    double d2r_metric[4][4], d2r_N, d2r_omega;

    get_metric_fist_derivatives(e_metric, dr_metric, &dr_N, &dr_omega, r, M_PI_2,
        Kerr_class, RBH_class, Wormhole_class);

    get_metric_second_derivatives(e_metric, d2r_metric, &d2r_N, &d2r_omega, r, M_PI_2,
        Kerr_class, RBH_class, Wormhole_class);

    double root = sqrt(dr_metric[0][3] * dr_metric[0][3] - dr_metric[0][0] * dr_metric[3][3]);

    return  -Kepler / dr_metric[3][3] * d2r_metric[3][3] + (-d2r_metric[0][3] + 1.0 / root / 2 * (2 * dr_metric[0][3] * d2r_metric[0][3] - dr_metric[0][0] * d2r_metric[3][3] - d2r_metric[0][0] * dr_metric[3][3])) / dr_metric[3][3];

}

double Novikov_Thorne_Model::Redshift(e_Spacetimes e_metric, double J, double State_Vector[], double r_obs, double theta_obs,
                                      c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double r_source = State_Vector[e_r];
    double theta_source = State_Vector[e_theta];

    /*
    Get the observer 4-velocity
    */

    double metric_obs[4][4], N_obs, omega_obs;

    get_metric(e_metric, metric_obs, &N_obs, &omega_obs, r_obs, theta_obs,
        Kerr_class, RBH_class, Wormhole_class);

    double U_obs[4] = { 1.0 / N_obs, 0 ,0 , omega_obs / N_obs };

    /*
    Get the source 4-velocity
    */

    double metric_source[4][4], N_source, omega_source;

    get_metric(e_metric, metric_source, &N_source, &omega_source, r_source, theta_source,
        Kerr_class, RBH_class, Wormhole_class);

    double Kepler = Keplerian_angular_velocity(e_metric, r_source, Kerr_class, RBH_class, Wormhole_class);

    double gamma = 1 / sqrt(-metric_source[0][0] - 2 * metric_source[0][3] * Kepler - metric_source[3][3] * Kepler * Kepler);

    double U_source[4] = { gamma, 0, 0, gamma * Kepler };

    /*
    Compute the redshift
    */

    return  (-U_obs[0] + U_obs[3] * J) / (-U_source[0] + U_source[3] * J);

}

double Novikov_Thorne_Model::disk_Energy(e_Spacetimes e_metric, double r,
                                         c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {


    double metric[4][4], N, omega;

    get_metric(e_metric, metric, &N, &omega, r, M_PI_2,
        Kerr_class, RBH_class, Wormhole_class);

    double Kepler = Keplerian_angular_velocity(e_metric, r,
        Kerr_class, RBH_class, Wormhole_class);

    double root = sqrt(-metric[0][0] - 2 * metric[0][3] * Kepler - metric[3][3] * Kepler * Kepler);

    return  -(metric[0][0] + metric[0][3] * Kepler) / root;

}

double Novikov_Thorne_Model::disk_Angular_Momentum(e_Spacetimes e_metric, double r,
                                                   c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double metric[4][4], N, omega;

    get_metric(e_metric, metric, &N, &omega, r, M_PI_2,
        Kerr_class, RBH_class, Wormhole_class);

    double Kepler = Keplerian_angular_velocity(e_metric, r,
        Kerr_class, RBH_class, Wormhole_class);

    double root = sqrt(-metric[0][0] - 2 * metric[0][3] * Kepler - metric[3][3] * Kepler * Kepler);

    return  (metric[3][3] * Kepler + metric[0][3]) / root;

}

double Novikov_Thorne_Model::Flux_integrand(e_Spacetimes e_metric, double r,
                                            c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double metric[4][4], N, omega;

    double dr_metric[4][4], dr_N, dr_omega;

    get_metric(e_metric, metric, &N, &omega, r, M_PI_2,
        Kerr_class, RBH_class, Wormhole_class);

    get_metric_fist_derivatives(e_metric, dr_metric, &dr_N, &dr_omega, r, M_PI_2,
        Kerr_class, RBH_class, Wormhole_class);

    double Kepler = Keplerian_angular_velocity(e_metric, r,
        Kerr_class, RBH_class, Wormhole_class);
    double dr_Kepler = dr_Keplerian_angular_velocity(e_metric, r, Kepler,
        Kerr_class, RBH_class, Wormhole_class);

    double metric_det = get_metric_det(metric);

    double root = sqrt(-metric[0][0] - 2 * metric[0][3] * Kepler - metric[3][3] * Kepler * Kepler);
    double dr_root = (-dr_metric[0][0] - 2 * (dr_metric[0][3] * Kepler + metric[0][3] * dr_Kepler) - dr_metric[3][3] * Kepler * Kepler - 2 * metric[3][3] * Kepler * dr_Kepler);

    double E = disk_Energy(e_metric, r, Kerr_class, RBH_class, Wormhole_class);
    double L = disk_Angular_Momentum(e_metric, r, Kerr_class, RBH_class, Wormhole_class);

    double dr_L = (dr_metric[3][3] * Kepler + metric[3][3] * dr_Kepler + dr_metric[0][3]) / root - L / root / root / 2 * dr_root;

    return (E - Kepler * L) * dr_L;

}

double Novikov_Thorne_Model::solve_Flux_integral(e_Spacetimes e_metric, double lower_bound, double upper_bound, double tolerance,
                                                 c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double mid_point = (lower_bound + upper_bound) / 2;
    double left_mid_point = (lower_bound + mid_point) / 2;
    double right_mid_point = (mid_point + upper_bound) / 2;

    double F_lower_bound = Flux_integrand(e_metric, lower_bound,
        Kerr_class, RBH_class, Wormhole_class);
    double F_mid_point = Flux_integrand(e_metric, mid_point,
        Kerr_class, RBH_class, Wormhole_class);
    double F_upper_bound = Flux_integrand(e_metric, upper_bound,
        Kerr_class, RBH_class, Wormhole_class);

    double F_left_mid = Flux_integrand(e_metric, left_mid_point,
        Kerr_class, RBH_class, Wormhole_class);
    double F_right_mid = Flux_integrand(e_metric, right_mid_point,
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

        double L_value = solve_Flux_integral(e_metric, lower_bound, mid_point, tolerance / 2,
            Kerr_class, RBH_class, Wormhole_class);
        double R_value = solve_Flux_integral(e_metric, mid_point, upper_bound, tolerance / 2,
            Kerr_class, RBH_class, Wormhole_class);

        integral = L_value + R_value;

    }

    return integral;
}

double Novikov_Thorne_Model::get_flux(e_Spacetimes e_metric, double r, double r_in,
                                      c_Kerr Kerr_class, c_RBH RBH_class, c_Wormhole Wormhole_class) {

    double metric[4][4], N, omega;

    get_metric(e_metric, metric, &N, &omega, r, M_PI_2,
        Kerr_class, RBH_class, Wormhole_class);

    double metric_det = get_metric_det(metric);
    double E_disk = disk_Energy(e_metric, r, Kerr_class, RBH_class, Wormhole_class);
    double L_disk = disk_Angular_Momentum(e_metric, r, Kerr_class, RBH_class, Wormhole_class);

    double Kepler = Keplerian_angular_velocity(e_metric, r,
        Kerr_class, RBH_class, Wormhole_class);
    double dr_Kepler = dr_Keplerian_angular_velocity(e_metric, r, Kepler,
        Kerr_class, RBH_class, Wormhole_class);

    double Flux_coeff = -dr_Kepler / ((E_disk - Kepler * L_disk) * (E_disk - Kepler * L_disk)) / (4 * M_PI * sqrt(metric_det));

    double r_throat = Wormhole_class.get_r_throat();

    if (e_metric == Wormhole) {

        r = sqrt(r * r + r_throat * r_throat);

    }

    double Flux_integral = solve_Flux_integral(e_metric, r_in, r, INTEGRAL_ACCURACY,
        Kerr_class, RBH_class, Wormhole_class);

    return Flux_coeff * Flux_integral;

}

double Novikov_Thorne_Model::get_r_in() { return r_in; };

double Novikov_Thorne_Model::get_r_out() { return r_out; };

Optically_Thin_Toroidal_Model::tag_Optically_Thin_Toroidal_Model(double alpha, double height_scale, double rad_cutoff, double omega) {

    DISK_ALPHA = alpha;
    DISK_HEIGHT_SCALE = height_scale;
    DISK_RAD_CUTOFF = rad_cutoff;
    DISK_OMEGA = omega;

}

double Optically_Thin_Toroidal_Model::get_disk_alpha()        { return DISK_ALPHA; };
double Optically_Thin_Toroidal_Model::get_disk_height_scale() { return DISK_HEIGHT_SCALE; };
double Optically_Thin_Toroidal_Model::get_disk_rad_cutoff()   { return DISK_RAD_CUTOFF; };;
double Optically_Thin_Toroidal_Model::get_disk_omega()        { return DISK_OMEGA; };
