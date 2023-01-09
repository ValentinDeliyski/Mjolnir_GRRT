#pragma once

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

extern e_Spacetimes e_metric;
extern std::vector<c_Spacetime_Base*> Spacetimes;
extern c_Observer Observer_class;
extern Optically_Thin_Toroidal_Model OTT_Model;
extern Novikov_Thorne_Model NT_Model;
extern int texture_indexer;
extern float Max_Intensity;
extern bool Normalizing_colormap;

void Lens(Initial_conditions_type* s_Initial_Conditions, std::ofstream data[], std::ofstream momentum_data[]) {

    double& r_obs     = s_Initial_Conditions->init_Pos[e_r];
    double& theta_obs = s_Initial_Conditions->init_Pos[e_theta];
    double& phi_obs   = s_Initial_Conditions->init_Pos[e_phi];
    double& J         = s_Initial_Conditions->init_Three_Momentum[e_phi];
    double& p_theta_0 = s_Initial_Conditions->init_Three_Momentum[e_theta];
    double& p_r_0     = s_Initial_Conditions->init_Three_Momentum[e_r];

    // Initialize the struct that holds the ray results
    s_Results Ray_results{};

    for (int Image_order = direct; Image_order <= ORDER_NUM - 1; Image_order += 1) {

        Ray_results.Three_Momentum[e_phi][Image_order] = J;

    }

    Ray_results.Parameters[Kerr] = SPIN;
    Ray_results.Parameters[Wormhole] = WH_REDSHIFT;
    Ray_results.Parameters[Reg_Black_Hole] = RBH_PARAM;
    Ray_results.Parameters[Naked_Singularity] = JNW_GAMMA;

    // Initialize initial State Vector
    double State_vector[e_State_Number] = { r_obs, theta_obs, phi_obs, 0., p_theta_0, p_r_0, 0., 0. };

    // Initialize a vector that stores the old state
    double Old_state[e_State_Number]{};

    // Initialize arrays that store the states and derivatives of the intermidiate integration steps
    double Derivatives[RK45_size * e_State_Number]{};

    // Set the old State Vector and the Test State Vector to the Initial State Vector
    for (int vector_indexer = e_r; vector_indexer <= e_p_r; vector_indexer += 1) {

        Old_state[vector_indexer] = State_vector[vector_indexer];

    }

    // Initialize two 3D position vectors, used for computing the photon tangent at the observer
    double r2[3]{}, r1[3]{
                          State_vector[e_r] * cos(State_vector[e_phi]) * sin(State_vector[e_theta]),
                          State_vector[e_r] * sin(State_vector[e_phi]) * sin(State_vector[e_theta]),
                          State_vector[e_r] * cos(State_vector[e_theta])
    };

    // Initialize counters for the Number Of Integration Steps, the Image Order and the Number Of Equator Crossings
    int integration_count{}, n_equator_crossings{}, Image_Order[DISK_MODEL_NUM]{};

    // Initialize the logical flags and error enums
    bool continue_integration = false;

    // Calculate the image coordinates from the initial conditions
    get_impact_parameters(s_Initial_Conditions, Ray_results.Image_Coords);

    Step_controller controller(INIT_STEPSIZE);

    while (integration_count < MAX_INTEGRATION_COUNT) {

        RK45(State_vector, Derivatives, J, &controller);

        // If error estimate, returned from RK45_EOM < RK45_ACCURACY
        if (controller.continue_integration) {

            // Novikov-Thorne Model Evaluation

            double crossing_coords[3];

            bool inside_NT_disk = State_vector[e_r] * State_vector[e_r] > NT_Model.get_r_in()  * NT_Model.get_r_in() &&
                                  State_vector[e_r] * State_vector[e_r] < NT_Model.get_r_out() * NT_Model.get_r_out() &&
                                  crossed_equatior(State_vector, Old_state);

            if (interpolate_crossing(State_vector, Old_state, crossing_coords)) {

                Image_Order[Novikov_Thorne] = n_equator_crossings;

                double r_crossing = sqrt(crossing_coords[0] * crossing_coords[0] + crossing_coords[1] * crossing_coords[1]);
                double state_crossing[2] = { r_crossing, 3.141 / 2 };

                Ray_results.Redshift_NT[Image_Order[Novikov_Thorne]]             = NT_Model.Redshift(J, state_crossing, r_obs, theta_obs, Spacetimes);
                Ray_results.Flux_NT[Image_Order[Novikov_Thorne]]                 = NT_Model.get_flux(r_crossing, Spacetimes);
                Ray_results.Three_Momentum[e_r][Image_Order[Novikov_Thorne]]     = State_vector[e_p_r];
                Ray_results.Three_Momentum[e_theta][Image_Order[Novikov_Thorne]] = State_vector[e_p_theta];

            }

            if (crossed_equatior(State_vector, Old_state)) {

                if (n_equator_crossings < ORDER_NUM - 1) {

                    n_equator_crossings += 1;

                }
            }

            // Evaluate logical flags for terminating the integration

            if (Spacetimes[e_metric]->terminate_integration(State_vector, Derivatives)) {

                Ray_results.Intensity[direct] = State_vector[e_Intensity];
                Ray_results.Optical_Depth     = State_vector[e_Optical_Depth];

                if (!Normalizing_colormap) {

                    set_pixel_color(State_vector[e_Intensity], texture_indexer);
                    Ray_results.Source_Coords[e_theta][direct] = State_vector[e_theta];
                    Ray_results.Source_Coords[e_phi][direct]   = State_vector[e_phi];

                    write_to_file(Ray_results, data, momentum_data);
                    
                    texture_indexer += 3;

                }
                else if(State_vector[e_Intensity] > Max_Intensity) {

                    Max_Intensity = State_vector[e_Intensity];

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

}

