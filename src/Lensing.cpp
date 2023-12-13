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
    // because its not parat of the photon state vector - I take care of it after the loop "by hand"

    for (int index = e_r; index <= e_path_log_number - 2; index++) {

        s_Ray_Results->Ray_log_struct.Ray_path_log[index + log_offset * e_path_log_number] = State_Vector[index];

    }

    s_Ray_Results->Ray_log_struct.Ray_path_log[e_path_log_step + log_offset * e_path_log_number] = Controller.previous_step;

}

void log_ray_emission(double Intensity, double Optical_depth, Results_type* s_Ray_Results, int log_index) {

    s_Ray_Results->Ray_log_struct.Ray_emission_log[0 + 2 * log_index] = Intensity;
    s_Ray_Results->Ray_log_struct.Ray_emission_log[1 + 2 * log_index] = Optical_depth;

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

        if (s_Ray_results->Redshift_NT[Image_Order] != 0) {

            s_Ray_results->Flux_NT[Image_Order] = 0;

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
                                double const Intensity) {

        int Image_Order = compute_image_order(Max_theta_turning_points - N_theta_turning_points, s_Initial_Conditions);

        p_Ray_results->Intensity[Image_Order] = Intensity;

        for (int order_scan = Image_Order + 1; order_scan <= compute_image_order(Max_theta_turning_points, s_Initial_Conditions); order_scan++) {

            p_Ray_results->Intensity[Image_Order] -= p_Ray_results->Intensity[order_scan];

        }

}

void Propagate_forward_emission(Initial_conditions_type* const s_Initial_Conditions, Results_type* const s_Ray_results, int* const N_theta_turning_points, int integration_count) {

    int const Max_theta_turning_points = *N_theta_turning_points;
    *N_theta_turning_points = 0;

    int Image_Order = compute_image_order(Max_theta_turning_points, s_Initial_Conditions);

    double Intensity{};
    double Optical_Depth{};

    double* Logged_ray_path[INTERPOLATION_NUM];
    double* U_source_coord[INTERPOLATION_NUM];

    double redshift[INTERPOLATION_NUM]{};
    double emission_function[INTERPOLATION_NUM]{};
    double absorbtion_function[INTERPOLATION_NUM]{};

    /* Initialize the arrays that hold the "current" emission function, 
    because in the loop only the "next" ones are evaluated, and then reused as the "current" on the next iteration */

    Logged_ray_path[Current] = &(s_Ray_results->Ray_log_struct.Ray_path_log[integration_count * e_path_log_number]);

    U_source_coord[Current]      = s_Initial_Conditions->OTT_model->get_disk_velocity(Logged_ray_path[Current], s_Initial_Conditions);
    redshift[Current]            = Redshift(Logged_ray_path[Current], U_source_coord[Current]);

    emission_function[Current]   = s_Initial_Conditions->OTT_model->get_emission_function_synchotron_exact(Logged_ray_path[Current], s_Initial_Conditions);
    absorbtion_function[Current] = s_Initial_Conditions->OTT_model->get_absorbtion_function(emission_function[Current], Logged_ray_path[Current], redshift[Current], OBS_FREQUENCY_CGS / redshift[Current]);
    
    for (int index = integration_count; index > 0; index--) {

        /* Pick out the ray position / momenta from the Log, at the given log index */

        Logged_ray_path[Current]  = &(s_Ray_results->Ray_log_struct.Ray_path_log[index * e_path_log_number]);
        Logged_ray_path[Next]     = &(s_Ray_results->Ray_log_struct.Ray_path_log[(index - 1) * e_path_log_number]);

        *N_theta_turning_points += Increment_theta_turning_points(Logged_ray_path[Current], Logged_ray_path[Next]);

        /* Reset the ebitmask that checks weather the metric was calculated this step */

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

        /* Get the next emission function and interpolate with the current point along the ray */

        emission_function[Next] = s_Initial_Conditions->OTT_model->get_emission_function_synchotron_exact(Logged_ray_path[Next], s_Initial_Conditions);
        double const interpolated_emission_function = 0.5 * (emission_function[Current] + emission_function[Next]);

        /* Get the next absorbtion and interpolate with the current point along the ray */

        absorbtion_function[Next] = s_Initial_Conditions->OTT_model->get_absorbtion_function(emission_function[Next], Logged_ray_path[Next], redshift[Next], OBS_FREQUENCY_CGS / redshift[Next]);
        double const interpolated_absorbtion_function = 0.5 * (absorbtion_function[Current] + absorbtion_function[Next]);

        double interpolation_step = Logged_ray_path[Current][e_path_log_step] / 2;

        if (absorbtion_function[Current] != 0) {

            double absorbtion_exponent = exp(-absorbtion_function[Current] / redshift[Current] * interpolation_step * MASS_TO_CM);

            Intensity *= absorbtion_exponent;
            Intensity += redshift[Current] * redshift[Current] * redshift[Current] * emission_function[Current] / absorbtion_function[Current] * (1 - absorbtion_exponent);
            Optical_Depth = emission_function[Current] * Logged_ray_path[Current][e_path_log_step];

        }

        if (interpolated_absorbtion_function != 0) {

            double absorbtion_exponent = exp(-interpolated_absorbtion_function / interpolated_redshit * interpolation_step * MASS_TO_CM);

            Intensity *= absorbtion_exponent;
            Intensity += interpolated_redshit * interpolated_redshit * interpolated_redshit * interpolated_emission_function / interpolated_absorbtion_function * (1 - absorbtion_exponent);
            Optical_Depth = interpolated_absorbtion_function * interpolation_step;


        }

        /* Re-use the "next" emission functions as the "current" for the next iteration of the emission propagation */

        redshift[Current]            = redshift[Next];
        U_source_coord[Current]      = U_source_coord[Next];
        emission_function[Current]   = emission_function[Next];
        absorbtion_function[Current] = absorbtion_function[Next];

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

                Propagate_forward_emission(s_Initial_Conditions, &s_Ray_results, &N_theta_turning_points, integration_count);

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