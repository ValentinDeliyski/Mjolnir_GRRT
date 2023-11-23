#pragma once

#define _USE_MATH_DEFINES

#include "Enumerations.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "Disk_Models.h"
#include "IO_files.h"
#include "General_GR_functions.h"
#include "General_math_functions.h"
#include "Rendering_Engine.h"

#include "Lensing.h"
#include <iostream>


void log_ray(double State_Vector[], std::vector<double> *State_log, Results_type *s_Ray_Results) {

    for (int index = e_r; index <= e_State_Number - 1; index++) {

        State_log->at(index) = State_Vector[index];

    }

    s_Ray_Results->Ray_log.push_back(*State_log);

}

Results_type Propagate_ray(Initial_conditions_type* s_Initial_Conditions) {

    /*************************************************************************************************
    |                                                                                                |
    |   @ Description: Propagates one light ray, specified by the struct "s_Initial_Conditions",     |
    |     and stores the results in the files "s_Ray_results" struct                                 |
    |                                                                                                |
    |   @ Inputs:                                                                                    |
    |     * s_Initial_Conditions: Struct that holds the initial position and momenta of the photon   |
    |     * data: Pointer to a series of output files for the images                                 |
    |     * momentum_data: Pointer to a series of ouput files for the momenta                        |
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

    // Initialize the struct that holds the ray results
    Results_type s_Ray_results{};

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

    std::vector<double> State_log(e_State_Number, 0);

    log_ray(State_vector, &State_log, &s_Ray_results);

    // Initialize a vector that stores the old state
    double Old_state[e_State_Number]{};

    // Initialize arrays that store the states and derivatives of the intermidiate integration steps
    double Derivatives[RK45_size * e_State_Number]{};

    // Set the old State Vector and the Test State Vector to the Initial State Vector
    for (int vector_indexer = e_r; vector_indexer <= e_p_r; vector_indexer += 1) {

        Old_state[vector_indexer] = State_vector[vector_indexer];

    }

    // Initialize counters for the Number Of Integration Steps, the Image Order and the Number Of Equator Crossings
    int integration_count{}, Image_Order[DISK_MODEL_NUM]{}, N_theta_turning_points{};

    // Calculate the image coordinates from the initial conditions
    get_impact_parameters(s_Initial_Conditions, s_Ray_results.Image_Coords);

    Step_controller controller(INIT_STEPSIZE);

    while (true) {

        RK45(State_vector, Derivatives, &controller, s_Initial_Conditions);

        // If error estimate, returned from RK45_EOM < RK45_ACCURACY
        if (controller.continue_integration) {

            // Novikov-Thorne Model Evaluation

            double crossing_coords[3]{}, crossing_momenta[3]{};

            N_theta_turning_points += Increment_theta_turning_points(State_vector, Old_state);

            if (Evaluate_NT_disk && interpolate_crossing(State_vector, Old_state, crossing_coords, crossing_momenta, s_Initial_Conditions->NT_model)) {

                Image_Order[Novikov_Thorne] = compute_image_order(N_theta_turning_points, s_Initial_Conditions);

                double r_crossing = std::sqrt(crossing_coords[x] * crossing_coords[x] + crossing_coords[y] * crossing_coords[y]);
                double state_crossing[2] = { r_crossing, M_PI_2 };

                s_Ray_results.Redshift_NT[Image_Order[Novikov_Thorne]] = s_Initial_Conditions->NT_model->Redshift(J, state_crossing, r_obs, theta_obs, s_Initial_Conditions->Spacetimes);
                
                if (s_Ray_results.Redshift_NT[Image_Order[Novikov_Thorne]] != 0) {

                    s_Ray_results.Flux_NT[Image_Order[Novikov_Thorne]] = 0;

                }

                s_Ray_results.Source_Coords[e_r][Image_Order[Novikov_Thorne]]      = r_crossing;
                s_Ray_results.Source_Coords[e_phi][Image_Order[Novikov_Thorne]]    = State_vector[e_phi];

                s_Ray_results.Three_Momentum[e_r][Image_Order[Novikov_Thorne]]     = crossing_momenta[e_r];
                s_Ray_results.Three_Momentum[e_theta][Image_Order[Novikov_Thorne]] = crossing_momenta[e_theta];

            }

            if (compute_image_order(N_theta_turning_points, s_Initial_Conditions) != Image_Order[Optically_Thin_Toroidal]) {

                s_Ray_results.Intensity[Image_Order[Optically_Thin_Toroidal]] = Old_state[e_Intensity]; 
                
                for (int order_scan = Image_Order[Optically_Thin_Toroidal] - 1; order_scan >= 0; order_scan -= 1) {

                    s_Ray_results.Intensity[Image_Order[Optically_Thin_Toroidal]] -= s_Ray_results.Intensity[order_scan];

                }      

                Image_Order[Optically_Thin_Toroidal] = compute_image_order(N_theta_turning_points, s_Initial_Conditions);

            }

            if (Active_Sim_Mode == 4) {

                log_ray(State_vector, &State_log, &s_Ray_results);

            }

            // Evaluate logical flags for terminating the integration

            if (s_Initial_Conditions->Spacetimes[e_metric]->terminate_integration(State_vector, Derivatives) ||
                integration_count >= MAX_INTEGRATION_COUNT) {    

                Image_Order[Optically_Thin_Toroidal] = compute_image_order(N_theta_turning_points, s_Initial_Conditions);

                s_Ray_results.Intensity[Image_Order[Optically_Thin_Toroidal]] = State_vector[e_Intensity];

                for (int order_scan = Image_Order[Optically_Thin_Toroidal] - 1; order_scan >= 0; order_scan -= 1) {

                    s_Ray_results.Intensity[Image_Order[Optically_Thin_Toroidal]] -= s_Ray_results.Intensity[order_scan];

                }

                s_Ray_results.Optical_Depth = State_vector[e_Optical_Depth];

                if (integration_count >= MAX_INTEGRATION_COUNT) {

                    std::cout << "Max iterations reached!" << '\n';

                }

                integration_count = 0;

                break;

            }

            for (int vector_indexer = 0; vector_indexer <= e_State_Number - 1; vector_indexer += 1) {

                Old_state[vector_indexer] = State_vector[vector_indexer];

            }

            integration_count += 1;

        }

    }

    return s_Ray_results;

}