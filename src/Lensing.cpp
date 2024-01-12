#pragma once

#define _USE_MATH_DEFINES

#include "Enumerations.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "Disk_Models.h"
#include "IO_files.h"
#include <iostream>
#include "General_GR_functions.h"
#include "General_math_functions.h"
#include "Rendering_Engine.h"

#include "Lensing.h"
#include <iostream>

void log_ray_path(double State_Vector[], Results_type* s_Ray_Results, Step_controller Controller){

    int &log_offset = s_Ray_Results->Ray_log_struct.Log_offset;

    // The loop continues up untill e_path_log_number - 2 in order to exclude the log step from the loop,
    // because its not parat of the photon state vector - I take care of it after the loop "by hand".

    for (int index = e_r; index <= e_path_log_number - 2; index++) {

        s_Ray_Results->Ray_log_struct.Ray_path_log[index + log_offset * e_path_log_number] = State_Vector[index];

    }

    s_Ray_Results->Ray_log_struct.Ray_path_log[e_path_log_step + log_offset * e_path_log_number] = Controller.previous_step;

    // The wormhole metric works with a "global" radial coordinate, that goes negative on the other side of the throat.
    // The emission model can't work with this coordinate, so I log the normal spherical radial coordinate instead.

    if (e_metric == Wormhole) {

        double& WH_radial_coord = s_Ray_Results->Ray_log_struct.Ray_path_log[e_r + log_offset * e_path_log_number];

        WH_radial_coord = sqrt(WH_R_THROAT * WH_R_THROAT + WH_radial_coord * WH_radial_coord);

    }

}

void log_ray_emission(double Intensity[STOKES_PARAM_NUM], double Optical_depth, Results_type* s_Ray_Results, int log_index) {
    
    for (int stokes_idx = 0; stokes_idx <= STOKES_PARAM_NUM - 1; stokes_idx++) {

        s_Ray_Results->Ray_log_struct.Ray_emission_log[0 + 2 * log_index][stokes_idx] = Intensity[stokes_idx];
        s_Ray_Results->Ray_log_struct.Ray_emission_log[1 + 2 * log_index][stokes_idx] = Optical_depth;
    }

}

void Evaluate_Novikov_Thorne_disk(Initial_conditions_type* const s_Initial_Conditions,
                                  Results_type* const s_Ray_results,
                                  double State_vector[],
                                  double Old_state[], 
                                  int N_theta_turning_points) {

    double crossing_coords[3]{}, crossing_momenta[3]{};

    if (interpolate_crossing(State_vector, Old_state, crossing_coords, crossing_momenta, s_Initial_Conditions->NT_model)) {

        int Image_Order = compute_image_order(N_theta_turning_points, s_Initial_Conditions);

        double r_crossing = std::sqrt(crossing_coords[x] * crossing_coords[x] + crossing_coords[y] * crossing_coords[y]);
        double state_crossing[2] = { r_crossing, M_PI_2 };

        s_Ray_results->Redshift_NT[Image_Order] = s_Initial_Conditions->NT_model->Redshift(State_vector[e_p_phi], 
                                                                                           state_crossing, 
                                                                                           r_obs, 
                                                                                           theta_obs, 
                                                                                           s_Initial_Conditions->Spacetimes);

        if (s_Ray_results->Redshift_NT[Image_Order] > std::numeric_limits<double>::min()) {

            s_Ray_results->Flux_NT[Image_Order] = s_Initial_Conditions->NT_model->get_flux(r_crossing, s_Initial_Conditions->Spacetimes);

        }

        s_Ray_results->Source_Coords[e_r][Image_Order] = r_crossing;
        s_Ray_results->Source_Coords[e_phi][Image_Order] = State_vector[e_phi];

        s_Ray_results->Three_Momentum[e_r][Image_Order] = crossing_momenta[e_r];
        s_Ray_results->Three_Momentum[e_theta][Image_Order] = crossing_momenta[e_theta];
    }

}

void Seperate_Image_into_orders(int const Max_theta_turning_points, 
                                int const N_theta_turning_points,
                                Initial_conditions_type* const s_Initial_Conditions, 
                                Results_type* const p_Ray_results,
                                double const Intensity[STOKES_PARAM_NUM]) {

        int Image_Order = compute_image_order(Max_theta_turning_points - N_theta_turning_points, s_Initial_Conditions);

        for (int stokes_idx = 0; stokes_idx <= STOKES_PARAM_NUM - 1; stokes_idx++) {

            p_Ray_results->Intensity[Image_Order][stokes_idx] = Intensity[stokes_idx];

            for (int order_scan = Image_Order + 1; order_scan <= compute_image_order(Max_theta_turning_points, s_Initial_Conditions); order_scan++) {

                p_Ray_results->Intensity[Image_Order][stokes_idx] -= p_Ray_results->Intensity[order_scan][stokes_idx];
                p_Ray_results->Intensity[Image_Order][stokes_idx] *= bool(p_Ray_results->Intensity[Image_Order][stokes_idx] > 0);

            }

        }

}

void Get_radiative_transfer_matrix(double const redshift, 
                                   double absorbtion_functions[STOKES_PARAM_NUM],
                                   double faradey_functions[STOKES_PARAM_NUM], 
                                   double step, 
                                   double Transfer_Operator[STOKES_PARAM_NUM][STOKES_PARAM_NUM],
                                   double Integrated_Transfer_Operator[STOKES_PARAM_NUM][STOKES_PARAM_NUM]) {

    /* The reference for this implementation is from appendix D in https://arxiv.org/pdf/1602.03184.pdf, originally derived in https://doi.org/10.1007/BF00165988 */

    for (int row_idx = 0; row_idx <= STOKES_PARAM_NUM - 1; row_idx++) {

        for (int colum_idx = 0; colum_idx <= STOKES_PARAM_NUM - 1; colum_idx++) {

            Transfer_Operator[row_idx][colum_idx] = 0;
            Integrated_Transfer_Operator[row_idx][colum_idx] = 0;

        }

    }


    // Here I define a bunch of references, because its going to get hairy if I don't...
    double* alpha = absorbtion_functions;
    double* rho   = faradey_functions;

    if (redshift > std::numeric_limits<double>::min()) {

        for (int index = 0; index <= STOKES_PARAM_NUM - 1; index++) {

            alpha[index] /= redshift;
            rho[index] /= redshift;

        }

    }
    else {

        return;

    }

    /* These are the variables defined in D8 - D13, used in calculating the M matricies */
    
    double const alpha_squared = alpha[Q] * alpha[Q] + 
                                 alpha[U] * alpha[U] + 
                                 alpha[V] * alpha[V];
    
    double const rho_squared = rho[Q] * rho[Q] +
                               rho[U] * rho[U] +
                               rho[V] * rho[V];

    double const alpha_rho = alpha[Q] * rho[Q] +
                             alpha[U] * rho[U] +
                             alpha[V] * rho[V];

    // sigma is the sign of the variable alpha_rho
    int const sigma = (alpha_rho >= 0) ? 1 : -1;

    // These quantities can go ever so sligtly negative, which physically should not happen, but nmerically it does.
    // This breaks the sqrt() functions, and so guards have to be put in place
    double Theta = (alpha_squared - rho_squared) * (alpha_squared - rho_squared) / 4 + alpha_rho * alpha_rho;

    if (Theta > 0) {

        Theta = sqrt(Theta);

    }
    else {

        Theta = 0;

    }
    
    double Lambda[2] = { (Theta / 2 + (alpha_squared - rho_squared) / 2),
                         (Theta / 2 - (alpha_squared - rho_squared) / 2) };

    if (Lambda[0] > 0 && Lambda[1] > 0) {

        Lambda[0] = sqrt(Lambda[0]);
        Lambda[1] = sqrt(Lambda[1]);

    }
    else {

        Lambda[0] = 0;
        Lambda[1] = 0;

    }


    /* Thesse are used in the "scaling factors" infront of the M matricies */

    double const exp_I = exp(-alpha[I] * step);

    double const cosh_term = cosh(Lambda[0] * step);
    double const cos_term  = cos(Lambda[1] * step);

    double const sinh_term = sinh(Lambda[0] * step);
    double const sin_term  = sin(Lambda[1] * step);

    /* ========================== M_1 Matrix calculation ========================== */

    double const M_1_scale_factor = exp_I * (cosh_term + cos_term) / 2;

    double const M_1[4][4] = { {1, 0, 0, 0},
                               {0, 1, 0, 0},
                               {0, 0, 1, 0},
                               {0, 0, 0, 1} };

    /* ========================== M_2 Matrix calculation ========================== */

    double M_2_scale_factor{};

    if (Theta > std::numeric_limits<double>::min()) {

        M_2_scale_factor = -exp_I * sin_term / Theta;
    }

    double const M_2[4][4] = { {                         0,                          ( Lambda[1] * alpha[Q] - sigma * Lambda[0] * rho[Q]), ( Lambda[1] * alpha[U] - sigma * Lambda[0] * rho[U]), ( Lambda[1] * alpha[V] - sigma * Lambda[0] * rho[V])},
                               {(Lambda[1] * alpha[Q] - sigma * Lambda[0] * rho[Q]),                           0,                          ( sigma * Lambda[0] * alpha[V] + Lambda[1] * rho[V]), (-sigma * Lambda[0] * alpha[U] - Lambda[1] * rho[U])},
                               {(Lambda[1] * alpha[U] - sigma * Lambda[0] * rho[U]), (-sigma * Lambda[0] * alpha[V] - Lambda[1] * rho[V]),                           0,                          ( sigma * Lambda[0] * alpha[Q] + Lambda[1] * rho[Q])},
                               {(Lambda[1] * alpha[V] - sigma * Lambda[0] * rho[V]), ( sigma * Lambda[0] * alpha[U] + Lambda[1] * rho[U]), (-sigma * Lambda[0] * alpha[Q] - Lambda[1] * rho[Q]),                           0                         } };

    /* ========================== M_3 Matrix calculation ========================== */

    double M_3_scale_factor{};

    if (Theta > std::numeric_limits<double>::min()) {

        M_3_scale_factor = -exp_I * sinh_term / Theta;
    }

    double const M_3[4][4] = { {						 0,							 ( Lambda[0] * alpha[Q] + sigma * Lambda[1] * rho[Q]), ( Lambda[0] * alpha[U] + sigma * Lambda[1] * rho[Q]), ( Lambda[0] * alpha[V] + sigma * Lambda[1] * rho[V])},
                               {(Lambda[0] * alpha[Q] + sigma * Lambda[1] * rho[Q]),	 		               0,                          (-sigma * Lambda[1] * alpha[V] + Lambda[0] * rho[V]), ( sigma * Lambda[1] * alpha[U] - Lambda[0] * rho[U])},
                               {(Lambda[0] * alpha[U] + sigma * Lambda[1] * rho[U]), ( sigma * Lambda[1] * alpha[V] - Lambda[0] * rho[V]),	                         0,	                         (-sigma * Lambda[1] * alpha[Q] + Lambda[0] * rho[Q])},
                               {(Lambda[0] * alpha[V] + sigma * Lambda[1] * rho[V]), (-sigma * Lambda[1] * alpha[U] + Lambda[0] * rho[U]), ( sigma * Lambda[1] * alpha[Q] - Lambda[0] * rho[Q]),                           0                         } };

    /* ========================== M_4 Matrix calculation ========================== */

    double M_4_scale_factor{};

    if (Theta > std::numeric_limits<double>::min()) {

        M_4_scale_factor = exp_I * (cosh_term - cos_term) / Theta;
    }

    double const M_4[4][4] = { {   (alpha_squared + rho_squared) / 2,                      (alpha[V] * rho[U] - alpha[U] * rho[V]),                                     (alpha[Q] * rho[V] - alpha[V] * rho[Q]),                                     (alpha[U] * rho[Q] - alpha[Q] * rho[U])},
                               {(alpha[U] * rho[V] - alpha[V] * rho[U]), (alpha[Q] * alpha[Q] + rho[Q] * rho[Q] - (alpha_squared + rho_squared) / 2),                   (alpha[Q] * alpha[U] + rho[Q] * rho[U]),                                     (alpha[V] * alpha[Q] + rho[V] * rho[Q])},
                               {(alpha[V] * rho[Q] - alpha[Q] * rho[V]),                   (alpha[Q] * alpha[U] + rho[Q] * rho[U]),                   (alpha[U] * alpha[U] + rho[U] * rho[U] - (alpha_squared + rho_squared) / 2),                   (alpha[U] * alpha[V] + rho[U] * rho[V])},
                               {(alpha[Q] * rho[U] - alpha[U] * rho[Q]),                   (alpha[V] * alpha[Q] + rho[V] * rho[Q]),                                     (alpha[U] * alpha[V] + rho[U] * rho[V]),                   (alpha[V] * alpha[V] + rho[V] * rho[V] - (alpha_squared + rho_squared) / 2)}};

    /* ========================== This is the formal operator O(s,s') - the solution to D1 ========================== */

    for (int row_idx = 0; row_idx <= STOKES_PARAM_NUM - 1; row_idx++) {

        for (int colum_idx = 0; colum_idx <= STOKES_PARAM_NUM - 1; colum_idx++) {

            Transfer_Operator[row_idx][colum_idx] = M_1_scale_factor * M_1[row_idx][colum_idx] +
                                                    M_2_scale_factor * M_2[row_idx][colum_idx] + 
                                                    M_3_scale_factor * M_3[row_idx][colum_idx] + 
                                                    M_4_scale_factor * M_4[row_idx][colum_idx];

        }

    }

    /* ========================== The intergral of O(s,s') for constant M matricies ========================== */

    /* This part of the implementation is adapted from equation (24) of https://academic.oup.com/mnras/article/475/1/43/4712230 */

    bool math_guard_1 = fabs(alpha[I] * alpha[I] - Lambda[0] * Lambda[0]) > std::numeric_limits<double>::min();
    bool math_guard_2 = fabs(alpha[I] * alpha[I] + Lambda[1] * Lambda[1]) > std::numeric_limits<double>::min();

    if (math_guard_1 && math_guard_2){

        double f_1 = 1.0 / (alpha[I] * alpha[I] - Lambda[0] * Lambda[0]);
        double f_2 = 1.0 / (alpha[I] * alpha[I] + Lambda[1] * Lambda[1]);

        for (int row_idx = 0; row_idx <= STOKES_PARAM_NUM - 1; row_idx++) {

            for (int colum_idx = 0; colum_idx <= STOKES_PARAM_NUM - 1; colum_idx++) {

            
                Integrated_Transfer_Operator[row_idx][colum_idx] = -Lambda[0] * f_1 * M_3[row_idx][colum_idx] + alpha[I] * f_1 / 2 * (M_1[row_idx][colum_idx] + M_4[row_idx][colum_idx]) +
                                                                    Lambda[1] * f_2 * M_2[row_idx][colum_idx] + alpha[I] * f_2 / 2 * (M_1[row_idx][colum_idx] - M_4[row_idx][colum_idx]) -
                                                                    exp_I * ( ( -Lambda[0] * f_1 * M_3[row_idx][colum_idx] + alpha[I] * f_1 / 2 * (M_1[row_idx][colum_idx] + M_4[row_idx][colum_idx]) ) * cosh_term + 
                                                                              ( -Lambda[1] * f_2 * M_2[row_idx][colum_idx] + alpha[I] * f_2 / 2 * (M_1[row_idx][colum_idx] - M_4[row_idx][colum_idx]) ) * cos_term + 
                                                                              ( -alpha[I] * f_2 * M_2[row_idx][colum_idx] - Lambda[1] * f_2 / 2 * (M_1[row_idx][colum_idx] - M_4[row_idx][colum_idx]) ) * sin_term -
                                                                              (  alpha[I] * f_1 * M_3[row_idx][colum_idx] - Lambda[0] * f_1 / 2 * (M_1[row_idx][colum_idx] + M_4[row_idx][colum_idx]) ) * sinh_term);
                                                                        
            }

        }

    }

}

void Propagate_Stokes_vector(double redshift, 
                             double emission_functions[STOKES_PARAM_NUM], 
                             double absorbtion_functions[STOKES_PARAM_NUM],
                             double faradey_functions[STOKES_PARAM_NUM], 
                             double step, 
                             double Intensity[STOKES_PARAM_NUM]) {

    double transfer_operator[4][4]{};
    double integrated_transfer_operator[4][4];

    Get_radiative_transfer_matrix(redshift, absorbtion_functions, faradey_functions, step, transfer_operator, integrated_transfer_operator);

    double Transfered_emission_vector[STOKES_PARAM_NUM]{};
    double temp_Intensity[STOKES_PARAM_NUM]{};

    // Placeholder vector for use in the mat_vec_multiply_4D() function
    for (int index = 0; index <= STOKES_PARAM_NUM - 1; index++) {

        temp_Intensity[index] = Intensity[index];

    }

    mat_vec_multiply_4D(integrated_transfer_operator, emission_functions, Transfered_emission_vector);
    mat_vec_multiply_4D(transfer_operator, temp_Intensity, Intensity);

    for (int index = 0; index <= STOKES_PARAM_NUM - 1; index++) {

        Transfered_emission_vector[index] *= redshift * redshift;

        Intensity[index] += Transfered_emission_vector[index];

    }
}

void Propagate_forward_emission(Initial_conditions_type* const s_Initial_Conditions, Results_type* const s_Ray_results, int* const N_theta_turning_points, int const log_length) {

    int const Max_theta_turning_points = *N_theta_turning_points;
    *N_theta_turning_points = 0;

    int Image_Order = compute_image_order(Max_theta_turning_points, s_Initial_Conditions);

    double Intensity[STOKES_PARAM_NUM]{};

    // TODO: Propagate this aswell
    double Optical_Depth{};

    double* Logged_ray_path[INTERPOLATION_NUM];
    double* U_source_coord[INTERPOLATION_NUM];

    double redshift[INTERPOLATION_NUM]{};
    double emission_functions[INTERPOLATION_NUM][STOKES_PARAM_NUM]{};
    double absorbtion_functions[INTERPOLATION_NUM][STOKES_PARAM_NUM]{};

    /* Initialize the arrays that hold the "current" emission function, 
    because in the loop only the "next" ones are evaluated, and then reused as the "current" on the next iteration */

    Logged_ray_path[Current] = &(s_Ray_results->Ray_log_struct.Ray_path_log[log_length * e_path_log_number]);

    U_source_coord[Current] = s_Initial_Conditions->OTT_model->get_disk_velocity(Logged_ray_path[Current], s_Initial_Conditions);
    redshift[Current]       = Redshift(Logged_ray_path[Current], U_source_coord[Current]);

    s_Initial_Conditions->OTT_model->get_emission_function_synchotron_exact(Logged_ray_path[Current], s_Initial_Conditions, emission_functions[Current]);
    s_Initial_Conditions->OTT_model->get_absorbtion_function(emission_functions[Current], Logged_ray_path[Current], redshift[Current], OBS_FREQUENCY_CGS / redshift[Current], absorbtion_functions[Current]);
    
    for (int index = log_length; index > 0; index--) {

        /* Pick out the ray position / momenta from the Log, at the given log index */

        Logged_ray_path[Current]  = &(s_Ray_results->Ray_log_struct.Ray_path_log[index * e_path_log_number]);
        Logged_ray_path[Next]     = &(s_Ray_results->Ray_log_struct.Ray_path_log[(index - 1) * e_path_log_number]);

        *N_theta_turning_points += Increment_theta_turning_points(Logged_ray_path[Current], Logged_ray_path[Next]);

        /* Reset the bitmask that checks weather the metric was calculated this step */

        s_Initial_Conditions->Spacetimes[e_metric]->reset_eval_bitmask();

        /* Get the "next" disk cooridinate velocity and interpolate with the "current" point along the ray */

        U_source_coord[Next] = s_Initial_Conditions->OTT_model->get_disk_velocity(Logged_ray_path[Next], s_Initial_Conditions);

        double const interpolated_U_source_coord[4] = { 0.5 * (U_source_coord[Current][e_t_coord]     + U_source_coord[Next][e_t_coord]),
                                                        0.5 * (U_source_coord[Current][e_r_coord]     + U_source_coord[Next][e_r_coord]),
                                                        0.5 * (U_source_coord[Current][e_theta_coord] + U_source_coord[Next][e_theta_coord]),
                                                        0.5 * (U_source_coord[Current][e_phi_coord]   + U_source_coord[Next][e_phi_coord]) };

        /* Get the "next" redshift and interpolate with the "current" point along the ray */

        redshift[Next] = Redshift(Logged_ray_path[Next], U_source_coord[Next]);
        double const interpolated_redshit = 0.5 * (redshift[Current] + redshift[Next]);

        /* Get the "next" emission function and interpolate with the current point along the ray */

        s_Initial_Conditions->OTT_model->get_emission_function_synchotron_exact(Logged_ray_path[Next], s_Initial_Conditions, emission_functions[Next]);

        double interpolated_emission_function[4] = { 0.5 * (emission_functions[Current][I] + emission_functions[Next][I]),
                                                     0.5 * (emission_functions[Current][Q] + emission_functions[Next][Q]),
                                                     0.5 * (emission_functions[Current][V] + emission_functions[Next][V]),
                                                     0.5 * (emission_functions[Current][U] + emission_functions[Next][U]) };

        /* Get the "next" absorbtion and interpolate with the current point along the ray */

        s_Initial_Conditions->OTT_model->get_absorbtion_function(emission_functions[Next], Logged_ray_path[Next], redshift[Next], OBS_FREQUENCY_CGS / redshift[Next], absorbtion_functions[Next]);

        double  interpolated_absorbtion_function[4] = { 0.5 * (absorbtion_functions[Current][I] + absorbtion_functions[Next][I]),
                                                        0.5 * (absorbtion_functions[Current][Q] + absorbtion_functions[Next][Q]),
                                                        0.5 * (absorbtion_functions[Current][V] + absorbtion_functions[Next][V]),
                                                        0.5 * (absorbtion_functions[Current][U] + absorbtion_functions[Next][U]) };

        double interpolation_step = Logged_ray_path[Current][e_path_log_step] * MASS_TO_CM / 2;

        // TODO: Add these to the simulation!
        double faradey_functions[STOKES_PARAM_NUM]{};

        // First half-step
        Propagate_Stokes_vector(redshift[Current],    emission_functions[Current],    absorbtion_functions[Current],    faradey_functions, interpolation_step, Intensity);

        // Second half-step
        Propagate_Stokes_vector(interpolated_redshit, interpolated_emission_function, interpolated_absorbtion_function, faradey_functions, interpolation_step, Intensity);

        /* Re-use the "next" emission and absorbtion functions as the "current" for the following iteration of the emission propagation */
        for (int stokes_index = 0; stokes_index <= STOKES_PARAM_NUM - 1; stokes_index++) {

            emission_functions[Current][stokes_index]   = emission_functions[Next][stokes_index];
            absorbtion_functions[Current][stokes_index] = absorbtion_functions[Next][stokes_index];

        }

        /* Re-use the "next" redshift functions and source velocity as the "current" ones for the following iteration of the emission propagation */

        redshift[Current]        = redshift[Next];
        U_source_coord[Current]  = U_source_coord[Next];

        log_ray_emission(Intensity, Optical_Depth, s_Ray_results, index);

        Seperate_Image_into_orders(Max_theta_turning_points, *N_theta_turning_points, s_Initial_Conditions, s_Ray_results, Intensity);

    }

}

Results_type* Propagate_ray(Initial_conditions_type* s_Initial_Conditions) {

    /*************************************************************************************************
    |                                                                                                |
    |   @ Description: Propagates one light ray, specified by the struct "s_Initial_Conditions",     |
    |     and stores the results in the "s_Ray_results" struct                                       |
    |                                                                                                |
    |   @ Inputs:                                                                                    |
    |     * s_Initial_Conditions: Struct that holds the initial position and momenta of the photon,  |
    |       as well as instances of the spacetime, and physical medium classes                       |
    |                                                                                                |
    |   @ Ouput: None                                                                                |
    |                                                                                                |
    *************************************************************************************************/

    double& r_obs     = s_Initial_Conditions->init_Pos[e_r];
    double& theta_obs = s_Initial_Conditions->init_Pos[e_theta];
    double& phi_obs   = s_Initial_Conditions->init_Pos[e_phi];
    double& J         = s_Initial_Conditions->init_Three_Momentum[e_phi];
    double& p_theta_0 = s_Initial_Conditions->init_Three_Momentum[e_theta];
    double& p_r_0     = s_Initial_Conditions->init_Three_Momentum[e_r];

    // Initialize the struct that holds the ray results (as static in order to not blow up the stack -> this must always be passed around as a pointer outside of this function!)
    static Results_type s_Ray_results{};

    // The final results must be manually set to 0s because the s_Ray_results struct is STATIC, and therefore not automatically re-initialized to 0s!
    memset(&s_Ray_results.Intensity,      0, sizeof(s_Ray_results.Intensity));
    memset(&s_Ray_results.Flux_NT,        0, sizeof(s_Ray_results.Flux_NT));
    memset(&s_Ray_results.Redshift_NT,    0, sizeof(s_Ray_results.Redshift_NT));
    memset(&s_Ray_results.Source_Coords,  0, sizeof(s_Ray_results.Source_Coords));
    memset(&s_Ray_results.Three_Momentum, 0, sizeof(s_Ray_results.Three_Momentum));

    for (int Image_order = direct; Image_order <= ORDER_NUM - 1; Image_order += 1) {

        s_Ray_results.Three_Momentum[e_phi][Image_order] = J;

    }

    s_Ray_results.Parameters[Kerr]              = SPIN;
    s_Ray_results.Parameters[Wormhole]          = WH_REDSHIFT;
    s_Ray_results.Parameters[Reg_Black_Hole]    = RBH_PARAM;
    s_Ray_results.Parameters[Naked_Singularity] = JNW_GAMMA;
    s_Ray_results.Parameters[Gauss_Bonnet]      = GAUSS_BONNET_GAMMA;
    s_Ray_results.Parameters[BH_w_Dark_Matter]  = M_HALO / A_0;

    // Initialize the State Vector
    double State_vector[] = {r_obs, theta_obs, phi_obs, J, p_theta_0, p_r_0, 0., 0.};

    // Initialize a vector that stores the old state
    double Old_state[e_State_Number]{};

    // Set the old State Vector and the Test State Vector to the Initial State Vector
    for (int vector_indexer = e_r; vector_indexer <= e_p_r; vector_indexer += 1) {

        Old_state[vector_indexer] = State_vector[vector_indexer];

    }

    // Initialize counters for the Number Of Integration Steps, the Image Order and the Number Of Equator Crossings
    int integration_count{}, Image_Order{}, N_theta_turning_points{};

    // Calculate the image coordinates from the initial conditions
    get_impact_parameters(s_Initial_Conditions, s_Ray_results.Image_Coords);

    Step_controller controller(INIT_STEPSIZE);

    s_Ray_results.Ray_log_struct.Log_offset = 0;
    log_ray_path(State_vector, &s_Ray_results, controller);

    while (true) {

        RK45(State_vector, &controller, s_Initial_Conditions);

        // If error estimate, returned from RK45 < RK45_ACCURACY
        if (controller.continue_integration) {

            integration_count += 1;
            s_Ray_results.Ray_log_struct.Log_offset = integration_count;

            log_ray_path(State_vector, &s_Ray_results, controller);

            N_theta_turning_points += Increment_theta_turning_points(State_vector, Old_state);

            if (Evaluate_NT_disk){

                Evaluate_Novikov_Thorne_disk(s_Initial_Conditions, &s_Ray_results, State_vector, Old_state, N_theta_turning_points);
                
            }

            // Evaluate logical flags for terminating the integration

            if (controller.integration_complete || integration_count >= MAX_INTEGRATION_COUNT) { 

                if (integration_count >= MAX_INTEGRATION_COUNT) {

                    std::cout << "Max iterations reached!" << '\n';

                }

                s_Ray_results.Ray_log_struct.Log_length = integration_count;

                Propagate_forward_emission(s_Initial_Conditions, &s_Ray_results, &N_theta_turning_points, s_Ray_results.Ray_log_struct.Log_length);

                integration_count = 0;

                break;

            }

            for (int vector_indexer = 0; vector_indexer <= e_State_Number - 1; vector_indexer += 1) {

                Old_state[vector_indexer] = State_vector[vector_indexer];

            }

        }

    }

    return &s_Ray_results;

}