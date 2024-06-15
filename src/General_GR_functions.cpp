﻿#pragma once

#include "General_GR_functions.h"
#include "Constants.h"
#include "Spacetimes.h"

#include <cmath>

extern Observer_class Observer;

/* Metric Related Functions */

void invert_metric(double inv_metric[4][4], double metric[4][4]) {

    double g2 = metric[0][3] * metric[0][3] - metric[3][3] * metric[0][0];

    inv_metric[0][0] = -metric[3][3] / g2;
    inv_metric[0][3] = metric[0][3] / g2;
    inv_metric[3][0] = inv_metric[0][3];
    inv_metric[1][1] = 1. / metric[1][1];
    inv_metric[2][2] = 1. / metric[2][2];
    inv_metric[3][3] = -metric[0][0] / g2;

}

double get_metric_det(double metric[4][4]) {

    return metric[0][0] * metric[1][1] * metric[2][2] * metric[3][3] * (1 - metric[0][3] * metric[0][3] / (metric[0][0] * metric[3][3]));

}

double get_eq_induced_metric_det(double metric[4][4]) {

    return metric[0][0] * metric[1][1] * metric[3][3] * (1 - metric[0][3] * metric[0][3] / (metric[0][0] * metric[3][3]));

}

/* Basis Related Functions */

void get_ZAMO_tetrad(double e_t[4], double e_r[4], double e_theta[4], double e_phi[4], double metric[4][4]) {

    double ksi = 1.0 / sqrt(-metric[0][0] + metric[0][3] * metric[0][3] / metric[3][3]);
    double sigma = -metric[0][3] / metric[3][3] * ksi;

    e_t[e_t_coord]	   = ksi;
    e_t[e_r_coord]	   = 0;
    e_t[e_theta_coord] = 0;
    e_t[e_phi_coord]   = sigma;

    e_r[e_t_coord]	   = 0;
    e_r[e_r_coord]     = 1.0 / sqrt(metric[1][1]);
    e_r[e_theta_coord] = 0;
    e_r[e_phi_coord]   = 0;

    e_theta[e_t_coord]     = 0;
    e_theta[e_r_coord]     = 0;
    e_theta[e_theta_coord] = 1.0 / sqrt(metric[2][2]);
    e_theta[e_phi_coord]   = 0;

    e_phi[e_t_coord]     = 0;
    e_phi[e_r_coord]     = 0;
    e_phi[e_theta_coord] = 0;
    e_phi[e_phi_coord]   = 1.0 / sqrt(metric[3][3]);

}

void Contravariant_coord_to_ZAMO(double metric[4][4], double Contravariant_Vector[4], double ZAMO_Vector[4]) {

    // TODO: Make this more obvious

    double alpha = sqrt(-metric[0][0] + metric[0][3] * metric[0][3] / metric[3][3]);
    double beta = metric[0][3] / metric[3][3];

    double p_t = metric[0][0] * Contravariant_Vector[0] + metric[0][3] * Contravariant_Vector[3];
    double p_phi = metric[3][3] * Contravariant_Vector[3] + metric[0][3] * Contravariant_Vector[0];

    ZAMO_Vector[0] = -(1 / alpha * p_t - beta / alpha * p_phi);
    ZAMO_Vector[1] = sqrt(metric[1][1]) * Contravariant_Vector[1];
    ZAMO_Vector[2] = sqrt(metric[2][2]) * Contravariant_Vector[2];
    ZAMO_Vector[3] = 1 / sqrt(metric[3][3]) * p_phi;

}

/* Other Functions */

void get_intitial_conditions_from_angles(Initial_conditions_type* p_Initial_Conditions, double V_angle_cam, double H_angle_cam) {

    /*
    
     n_cam is the direction vector of the light ray in the camera frame ( {r, theta, phi} components)
     n_FIDO is the "rotated" observer (a.e. his camera "y" axis is NOT aligned with the spin axis of the central object)
    
    */

    double n_cam[3] = { cos(V_angle_cam) * cos(H_angle_cam), sin(V_angle_cam), sin(H_angle_cam) * cos(V_angle_cam)}; 
    double n_FIDO[3]{};


    n_FIDO[e_phi]   =  n_cam[e_phi] * cos(obs_cam_rotation_angle) + n_cam[e_theta] * sin(obs_cam_rotation_angle);
    n_FIDO[e_theta] = -n_cam[e_phi] * sin(obs_cam_rotation_angle) + n_cam[e_theta] * cos(obs_cam_rotation_angle);
    n_FIDO[e_r]     =  n_cam[e_r];

    double V_angle = asin(n_FIDO[e_theta]);
    double H_angle = atan2(n_FIDO[e_phi], n_FIDO[e_r]);

    double g2, gamma, ksi, L_z, E;

    double(*metric)[4] = p_Initial_Conditions->init_metric;

    g2 = pow(metric[0][3], 2) - metric[0][0] * metric[3][3];
    ksi = sqrt(metric[3][3] / g2);
    gamma = -metric[0][3] / metric[3][3] * ksi;

    L_z = sqrt(metric[3][3]) * sin(H_angle + 2 * M_PI) * cos(V_angle);
    E = (1 + gamma * L_z) / ksi;

    p_Initial_Conditions->init_Three_Momentum[e_phi]   = L_z / E;
    p_Initial_Conditions->init_Three_Momentum[e_theta] = sqrt(metric[2][2]) * sin(V_angle) / E;
    p_Initial_Conditions->init_Three_Momentum[e_r]     = sqrt(metric[1][1]) * cos(H_angle + 2 * M_PI) * cos(V_angle) / E;

}

void get_impact_parameters(Initial_conditions_type* p_Initial_Conditions, double Image_coords[2]) {

    double(*metric)[4] = p_Initial_Conditions->init_metric;

    double g2, gamma, ksi;

    g2 = pow(metric[0][3], 2) - metric[0][0] * metric[3][3];
    ksi = sqrt(metric[3][3] / g2);
    gamma = -metric[0][3] / metric[3][3] * ksi;

    double& r_0 = p_Initial_Conditions->init_Pos[e_r];
    double& J = p_Initial_Conditions->init_Three_Momentum[e_phi];
    double& p_th = p_Initial_Conditions->init_Three_Momentum[e_theta];

    Image_coords[x] = -r_0 *  J   / (ksi - gamma * J) / sqrt(metric[3][3]);
    Image_coords[y] =  r_0 * p_th / (ksi - gamma * J) / sqrt(metric[2][2]);

}

double Redshift(double State_Vector[], double U_source[]) {

    /******************************************************************************
    |																			  |
    |   @ Description: Computes redshift for for a ray at a point, specified by   |
    |	State_Vector, for an observer, specified by Observer_class				  |
    |																			  |
    |   @ Inputs:                                                                 |
    |     * J: Covariant, azimuthal component of the ray / photon 4 - momentum	  |
    |	  * State_Vector: State vector of the ray / photon						  |
    |	  * U_source: Contravariant 4 - velocity of the emmiting medium			  |
    |																			  |
    |   @ Ouput: The redshift of the ray / photon								  |
    |                                                                             |
    ******************************************************************************/

    double U_obs[4];

    Observer.get_obs_velocity(U_obs);

    double redshift = (U_obs[0] * (-1) + U_obs[3] * State_Vector[e_p_phi]) / (U_source[0] * (-1) +
                                                                              U_source[1] * State_Vector[e_p_r] +
                                                                              U_source[2] * State_Vector[e_p_theta] +
                                                                              U_source[3] * State_Vector[e_p_phi]);

    return redshift;

}

void Lorentz_boost_matrix(double Boost_matrix[4][4], double U_source[4]) {

    /********************************************************************************
    |                                                                               |
    |   @ Description: Computes a Lorentz boost matrix, given a contravariant       | 
    |	4 - velocity (in an observers basis) U_cource, then stores it in Boost_matrix   |
    |																				|
    |   @ Inputs:																	|
    |     * Boost_matrix: Pointer to where we store the boost matrix				|
    |	  * U_source: ZAMO 4-velocity, that generates the boost						|
    |																				|
    |   @ Ouput: None																|
    |																				|
    ********************************************************************************/


    double V_r     = U_source[1] / U_source[0];
    double V_theta = U_source[2] / U_source[0];
    double V_phi   = U_source[3] / U_source[0];

    double V_squared = V_r * V_r + V_theta * V_theta + V_phi * V_phi;

    double gamma = 1.0 / sqrt(1 - V_squared);

    Boost_matrix[0][0] = gamma;
    Boost_matrix[0][1] = gamma * V_r;
    Boost_matrix[0][2] = gamma * V_theta;
    Boost_matrix[0][3] = gamma * V_phi;

    for (int index = 1; index <= 3; index += 1) {

        Boost_matrix[index][0] = Boost_matrix[0][index];

    }

    Boost_matrix[1][1] = 1 + (gamma - 1) * V_r * V_r / V_squared;
    Boost_matrix[1][2] = (gamma - 1) * V_r * V_theta / V_squared;
    Boost_matrix[1][3] = (gamma - 1) * V_r * V_phi / V_squared;

    Boost_matrix[2][1] = Boost_matrix[1][2];
    Boost_matrix[3][1] = Boost_matrix[1][3];

    Boost_matrix[2][2] = 1 + (gamma - 1) * V_theta * V_theta / V_squared;
    Boost_matrix[2][3] = (gamma - 1) * V_theta * V_phi / V_squared;

    Boost_matrix[3][2] = Boost_matrix[2][3];

    Boost_matrix[3][3] = 1 + (gamma - 1) * V_phi * V_phi / V_squared;

}

int Increment_theta_turning_points(double State_Vector[], double Old_State[]) {

    if (State_Vector[e_p_theta] * Old_State[e_p_theta] < 0){

        return 1;
    
    }
    else {

        return 0;
    }

}

int compute_image_order(int N_theta_turning_points, Initial_conditions_type* p_Initial_Conditions) {

    int order = N_theta_turning_points - bool(p_Initial_Conditions->init_Three_Momentum[e_theta] > 0);

    if (order > 3) {

        order = 3;

    }

    return order * bool(order > 0);

}

void get_connection_coefficients(Metric_type s_Metric, Metric_type s_dr_metric, Metric_type s_dtheta_metric, double Connectrion_Coeffs[4][4][4]) {

    double inv_metric[4][4]{};

    invert_metric(inv_metric, s_Metric.Metric);

  /* ==================================================================== Г^t_{..} coefficients ================================================================ */

    Connectrion_Coeffs[e_t_coord][e_t_coord][e_t_coord]         = 0.0;
    Connectrion_Coeffs[e_t_coord][e_r_coord][e_r_coord]         = 0.0;
    Connectrion_Coeffs[e_t_coord][e_theta_coord][e_theta_coord] = 0.0;
    Connectrion_Coeffs[e_t_coord][e_phi_coord][e_phi_coord]     = 0.0;

    /* ------------------------------------------------------------------ Г^t_{t,r} coefficients --------------------------------------------------------------- */

    Connectrion_Coeffs[e_t_coord][e_t_coord][e_r_coord] = inv_metric[e_t_coord][e_t_coord]   * s_dr_metric.Metric[e_t_coord][e_t_coord] / 2 + 
                                                          inv_metric[e_t_coord][e_phi_coord] * s_dr_metric.Metric[e_t_coord][e_phi_coord] / 2;
    Connectrion_Coeffs[e_t_coord][e_r_coord][e_t_coord] = Connectrion_Coeffs[e_t_coord][e_t_coord][e_r_coord];

    /* ------------------------------------------------------------------ Г^t_{t,phi} coefficients ------------------------------------------------------------- */

    Connectrion_Coeffs[e_t_coord][e_t_coord][e_phi_coord] = 0.0;
    Connectrion_Coeffs[e_t_coord][e_phi_coord][e_t_coord] = 0.0;

    /* ------------------------------------------------------------------ Г^t_{t,theta} coefficients ----------------------------------------------------------- */

    Connectrion_Coeffs[e_t_coord][e_t_coord][e_theta_coord] = inv_metric[e_t_coord][e_t_coord]   * s_dtheta_metric.Metric[e_t_coord][e_t_coord] / 2 + 
                                                              inv_metric[e_t_coord][e_phi_coord] * s_dtheta_metric.Metric[e_t_coord][e_phi_coord] / 2;
    Connectrion_Coeffs[e_t_coord][e_theta_coord][e_t_coord] = Connectrion_Coeffs[e_t_coord][e_t_coord][e_theta_coord];

   /* ------------------------------------------------------------------ Г^t_{phi,r} coefficients -------------------------------------------------------------- */

    Connectrion_Coeffs[e_t_coord][e_phi_coord][e_r_coord] = inv_metric[e_t_coord][e_t_coord]   * s_dr_metric.Metric[e_t_coord][e_phi_coord] / 2 + 
                                                            inv_metric[e_t_coord][e_phi_coord] * s_dr_metric.Metric[e_phi_coord][e_phi_coord] / 2;
    Connectrion_Coeffs[e_t_coord][e_r_coord][e_phi_coord] = Connectrion_Coeffs[e_t_coord][e_phi_coord][e_r_coord];

   /* ------------------------------------------------------------------ Г^t_{phi,theta} coefficients ---------------------------------------------------------- */

    Connectrion_Coeffs[e_t_coord][e_phi_coord][e_theta_coord] = inv_metric[e_t_coord][e_t_coord]   * s_dtheta_metric.Metric[e_t_coord][e_phi_coord] / 2 + 
                                                                inv_metric[e_t_coord][e_phi_coord] * s_dtheta_metric.Metric[e_phi_coord][e_phi_coord] / 2;
    Connectrion_Coeffs[e_t_coord][e_theta_coord][e_phi_coord] = Connectrion_Coeffs[e_t_coord][e_phi_coord][e_theta_coord];

   /* ------------------------------------------------------------------ Г^t_{r,theta} coefficients ------------------------------------------------------------- */

    Connectrion_Coeffs[e_t_coord][e_theta_coord][e_r_coord] = 0.0;
    Connectrion_Coeffs[e_t_coord][e_r_coord][e_theta_coord] = 0.0;


  /* ==================================================================== Г^r_{..} coefficients ==================================================================== */

    Connectrion_Coeffs[e_r_coord][e_t_coord][e_t_coord]         = -inv_metric[e_r_coord][e_r_coord] * s_dr_metric.Metric[e_t_coord][e_t_coord] / 2;
    Connectrion_Coeffs[e_r_coord][e_r_coord][e_r_coord]         =  inv_metric[e_r_coord][e_r_coord] * s_dr_metric.Metric[e_r_coord][e_r_coord] / 2;
    Connectrion_Coeffs[e_r_coord][e_theta_coord][e_theta_coord] = -inv_metric[e_r_coord][e_r_coord] * s_dr_metric.Metric[e_theta_coord][e_theta_coord] / 2;
    Connectrion_Coeffs[e_r_coord][e_phi_coord][e_phi_coord]     = -inv_metric[e_r_coord][e_r_coord] * s_dr_metric.Metric[e_phi_coord][e_phi_coord] / 2;

   /* ------------------------------------------------------------------ Г^r_{t,phi} coefficients ---------------------------------------------------------------- */

    Connectrion_Coeffs[e_r_coord][e_t_coord][e_phi_coord] = -inv_metric[e_r_coord][e_r_coord] * s_dr_metric.Metric[e_t_coord][e_phi_coord] / 2;
    Connectrion_Coeffs[e_r_coord][e_phi_coord][e_t_coord] = Connectrion_Coeffs[e_r_coord][e_t_coord][e_phi_coord];

   /* ------------------------------------------------------------------ Г^r_{t,r} coefficients ------------------------------------------------------------------ */

    Connectrion_Coeffs[e_r_coord][e_t_coord][e_r_coord] = 0.0;
    Connectrion_Coeffs[e_r_coord][e_r_coord][e_t_coord] = 0.0;

   /* ------------------------------------------------------------------ Г^r_{theta,r} coefficients -------------------------------------------------------------- */

    Connectrion_Coeffs[e_r_coord][e_theta_coord][e_r_coord] = inv_metric[e_r_coord][e_r_coord] * s_dtheta_metric.Metric[e_r_coord][e_r_coord] / 2;
    Connectrion_Coeffs[e_r_coord][e_r_coord][e_theta_coord] = Connectrion_Coeffs[e_r_coord][e_theta_coord][e_r_coord];

   /* ------------------------------------------------------------------ Г^r_{r,phi} coefficients ---------------------------------------------------------------- */

    Connectrion_Coeffs[e_r_coord][e_r_coord][e_phi_coord] = 0.0;
    Connectrion_Coeffs[e_r_coord][e_phi_coord][e_r_coord] = 0.0;

   /* ------------------------------------------------------------------ Г^r_{t,theta} coefficients -------------------------------------------------------------- */

    Connectrion_Coeffs[e_r_coord][e_t_coord][e_theta_coord] = 0.0;
    Connectrion_Coeffs[e_r_coord][e_theta_coord][e_t_coord] = 0.0;
     
   /* ------------------------------------------------------------------ Г^r_{phi,theta} coefficients ------------------------------------------------------------ */

    Connectrion_Coeffs[e_r_coord][e_phi_coord][e_theta_coord] = 0.0;
    Connectrion_Coeffs[e_r_coord][e_theta_coord][e_phi_coord] = 0.0;

    /* ==================================================================== Г^theta_{..} coefficients ================================================================ */

    Connectrion_Coeffs[e_theta_coord][e_t_coord][e_t_coord] = -inv_metric[e_theta_coord][e_theta_coord] * s_dtheta_metric.Metric[e_t_coord][e_t_coord] / 2;
    Connectrion_Coeffs[e_theta_coord][e_r_coord][e_r_coord] = 0.0;
    Connectrion_Coeffs[e_theta_coord][e_theta_coord][e_theta_coord] = inv_metric[e_theta_coord][e_theta_coord] * s_dtheta_metric.Metric[e_theta_coord][e_theta_coord] / 2;
    Connectrion_Coeffs[e_theta_coord][e_phi_coord][e_phi_coord] = 0.0;

    /* ------------------------------------------------------------------ Г^theta_{t,phi} coefficients ---------------------------------------------------------------- */

    Connectrion_Coeffs[e_theta_coord][e_t_coord][e_phi_coord] = -inv_metric[e_theta_coord][e_theta_coord] * s_dtheta_metric.Metric[e_t_coord][e_phi_coord] / 2;
    Connectrion_Coeffs[e_theta_coord][e_phi_coord][e_t_coord] = Connectrion_Coeffs[e_theta_coord][e_t_coord][e_phi_coord];

    /* ------------------------------------------------------------------ Г^theta_{t,r} coefficients ------------------------------------------------------------------ */

    Connectrion_Coeffs[e_theta_coord][e_t_coord][e_r_coord] = 0.0;
    Connectrion_Coeffs[e_theta_coord][e_r_coord][e_t_coord] = 0.0;

    /* ------------------------------------------------------------------ Г^theta_{theta,r} coefficients -------------------------------------------------------------- */

    Connectrion_Coeffs[e_theta_coord][e_theta_coord][e_r_coord] = inv_metric[e_theta_coord][e_theta_coord] * s_dr_metric.Metric[e_theta_coord][e_theta_coord] / 2;
    Connectrion_Coeffs[e_theta_coord][e_r_coord][e_theta_coord] = Connectrion_Coeffs[e_theta_coord][e_theta_coord][e_r_coord];

    /* ------------------------------------------------------------------ Г^theta_{r,phi} coefficients ---------------------------------------------------------------- */

    Connectrion_Coeffs[e_theta_coord][e_r_coord][e_phi_coord] = 0.0;
    Connectrion_Coeffs[e_theta_coord][e_phi_coord][e_r_coord] = 0.0;

    /* ------------------------------------------------------------------ Г^theta_{t,theta} coefficients -------------------------------------------------------------- */

    Connectrion_Coeffs[e_theta_coord][e_t_coord][e_theta_coord] = 0.0;
    Connectrion_Coeffs[e_theta_coord][e_theta_coord][e_t_coord] = 0.0;

    /* ------------------------------------------------------------------ Г^theta_{phi,theta} coefficients ------------------------------------------------------------ */

    Connectrion_Coeffs[e_theta_coord][e_phi_coord][e_theta_coord] = 0.0;
    Connectrion_Coeffs[e_theta_coord][e_theta_coord][e_phi_coord] = 0.0;

  /* ==================================================================== Г^phi_{..} coefficients ==================================================================== */

    Connectrion_Coeffs[e_phi_coord][e_t_coord][e_t_coord] = 0.0;
    Connectrion_Coeffs[e_phi_coord][e_r_coord][e_r_coord] = 0.0;
    Connectrion_Coeffs[e_phi_coord][e_theta_coord][e_theta_coord] = 0.0;
    Connectrion_Coeffs[e_phi_coord][e_phi_coord][e_phi_coord] = 0.0;

   /* ------------------------------------------------------------------ Г^phi_{t,phi} coefficients ---------------------------------------------------------------- */

    Connectrion_Coeffs[e_phi_coord][e_t_coord][e_phi_coord] = 0.0;
    Connectrion_Coeffs[e_phi_coord][e_phi_coord][e_t_coord] = 0.0;

   /* ------------------------------------------------------------------ Г^phi_{r,phi} coefficients ---------------------------------------------------------------- */

    Connectrion_Coeffs[e_phi_coord][e_r_coord][e_phi_coord] = inv_metric[e_phi_coord][e_phi_coord] * s_dr_metric.Metric[e_phi_coord][e_phi_coord] / 2 + 
                                                              inv_metric[e_phi_coord][e_t_coord] * s_dr_metric.Metric[e_phi_coord][e_t_coord] / 2;
    Connectrion_Coeffs[e_phi_coord][e_phi_coord][e_r_coord] = Connectrion_Coeffs[e_r_coord][e_r_coord][e_phi_coord];

   /* ------------------------------------------------------------------ Г^phi_{theta,phi} coefficients ------------------------------------------------------------ */

    Connectrion_Coeffs[e_phi_coord][e_theta_coord][e_phi_coord] = inv_metric[e_phi_coord][e_phi_coord] * s_dtheta_metric.Metric[e_phi_coord][e_phi_coord] / 2 +
                                                                  inv_metric[e_phi_coord][e_t_coord] * s_dtheta_metric.Metric[e_phi_coord][e_t_coord] / 2;
    Connectrion_Coeffs[e_phi_coord][e_phi_coord][e_theta_coord] = Connectrion_Coeffs[e_phi_coord][e_theta_coord][e_phi_coord];
    
   /* ------------------------------------------------------------------ Г^phi_{t,r} coefficients ------------------------------------------------------------------ */

    Connectrion_Coeffs[e_phi_coord][e_t_coord][e_r_coord] = inv_metric[e_phi_coord][e_phi_coord] * s_dr_metric.Metric[e_phi_coord][e_t_coord] / 2 +
                                                            inv_metric[e_phi_coord][e_t_coord] * s_dr_metric.Metric[e_t_coord][e_t_coord] / 2;
    Connectrion_Coeffs[e_phi_coord][e_r_coord][e_t_coord] = Connectrion_Coeffs[e_phi_coord][e_t_coord][e_r_coord];

   /* ------------------------------------------------------------------ Г^phi_{theta,r} coefficients -------------------------------------------------------------- */

    Connectrion_Coeffs[e_phi_coord][e_theta_coord][e_r_coord] = 0.0;
    Connectrion_Coeffs[e_phi_coord][e_r_coord][e_theta_coord] = 0.0;
    
   /* ------------------------------------------------------------------ Г^phi_{t,r} coefficients ------------------------------------------------------------------ */

    Connectrion_Coeffs[e_phi_coord][e_t_coord][e_theta_coord] = inv_metric[e_phi_coord][e_phi_coord] * s_dtheta_metric.Metric[e_phi_coord][e_t_coord] / 2 +
                                                                inv_metric[e_phi_coord][e_t_coord] * s_dtheta_metric.Metric[e_t_coord][e_t_coord] / 2;
    Connectrion_Coeffs[e_phi_coord][e_theta_coord][e_t_coord] = Connectrion_Coeffs[e_phi_coord][e_t_coord][e_theta_coord];

}