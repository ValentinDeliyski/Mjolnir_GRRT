#pragma once

#include "General_GR_functions.h"
#include "Constants.h"
#include <cmath>

extern std::vector<Spacetime_Base_Class*> Spacetimes;
extern c_Observer Observer_class;

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

    return metric[1][1] * (metric[0][3] * metric[0][3] - metric[0][0] * metric[3][3]);

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

void get_intitial_conditions_from_angles(Initial_conditions_type* p_Initial_Conditions, double V_angle, double H_angle) {

    double g2, gamma, ksi, L_z, E;

    double(*metric)[4] = p_Initial_Conditions->init_metric;

    g2 = pow(metric[0][3], 2) - metric[0][0] * metric[3][3];
    ksi = sqrt(metric[3][3] / g2);
    gamma = -metric[0][3] / metric[3][3] * ksi;

    L_z = sqrt(metric[3][3]) * sin(H_angle + 2 * M_PI) * cos(V_angle);
    E = (1 + gamma * L_z) / ksi;

    p_Initial_Conditions->init_Three_Momentum[e_phi] = L_z / E;
    p_Initial_Conditions->init_Three_Momentum[e_theta] = sqrt(metric[2][2]) * sin(V_angle) / E;
    p_Initial_Conditions->init_Three_Momentum[e_r] = sqrt(metric[1][1]) * cos(H_angle + 2 * M_PI) * cos(V_angle) / E;

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

double Redshift(double J, double State_Vector[], double U_source[]) {

    /******************************************************************************
    |																			  |
    |   @ Description: Computes redshift for for a ray at a point, specified by   |
    |	State_Vector, for an observer, specified bt Observer_class				  |
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

    Observer_class.get_obs_velocity(U_obs, Spacetimes);

    return  (-U_obs[0] + U_obs[3] * J) / (-U_source[0] * 1 +
                                           U_source[1] * State_Vector[e_p_r] +
                                           U_source[2] * State_Vector[e_p_theta] +
                                           U_source[3] * J);

}

double get_photon_t_momentum(double State_vector[], double J, double metric[4][4]) {

    /*

    State_vector holds covariant momentum components

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

void Lorentz_boost_matrix(double Boost_matrix[4][4], double U_source[4]) {

    /********************************************************************************
    |                                                                               |
    |   @ Description: Computes a Lorentz boost matrix, given a contravariant       | 
    |	4 - velocity (in the ZAMO basis) U_cource, then stores it in Boost_matrix   |
    |																				|
    |   @ Inputs:																	|
    |     * Boost_matrix: Pointer to where we store the boost matrix				|
    |	  * U_source: ZAMO 4-velocity, that generates the boost						|
    |																				|
    |   @ Ouput: None																|
    |																				|
    ********************************************************************************/


    double V_r_ZAMO = U_source[1] / U_source[0];
    double V_theta_ZAMO = U_source[2] / U_source[0];
    double V_phi_ZAMO = U_source[3] / U_source[0];

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

}

int compute_order_no_spin(double State_Vector[]) {

    /***************************************************************************
    |                                                                          |
    |   @ Description: Computes the order of the image, generated by the ray   |
    |                                                                          |
    |   @ Inputs:                                                              |
    |     * State_Vector: State vector of the ray / photon					   |
    |																		   |
    |   @ Ouput: The order of the image								           |
    |                                                                          |
    ***************************************************************************/

    double abs_phi   = sqrt(State_Vector[e_phi] * State_Vector[e_phi]);

        if (abs_phi <= M_PI) {

            return 0;

        }
        else if (abs_phi <= 2 * M_PI) {

            return 1;

        }
        else if (abs_phi <= 3 * M_PI) {

            return 2;

        }
        else {

            return 3;

        }

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

