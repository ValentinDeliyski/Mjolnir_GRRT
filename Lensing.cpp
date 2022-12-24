#pragma once

#include "Enumerations.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "Disk_Models.h"
#include "IO_files.h"
#include "General_math_functions.h"
#include "Rendering_Engine.h"

#include <iostream>

extern e_Spacetimes e_metric;
extern std::vector<c_Spacetime_Base*> Spacetimes;
extern c_Observer Observer_class;
extern Optically_Thin_Toroidal_Model OTT_Model;
extern Novikov_Thorne_Model NT_Model;
extern int texture_indexer;

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
    double inter_State_vector[RK45_size * e_State_Number]{}, Derivatives[RK45_size * e_State_Number]{};

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

    // Initialize the Initial Step Size
    double step = INIT_STEPSIZE;

    // Initialize the logical flags and error enums
    bool continue_integration = false;

    while (integration_count < MAX_INTEGRATION_COUNT) {

        RK45(State_vector, Derivatives, &step, J, &continue_integration);

        // If error estimate, returned from RK45_EOM < RK45_ACCURACY
        if (continue_integration) {

            // Initialize the light ray
            if (integration_count == 1) {

                r2[x] = State_vector[e_r] * cos(State_vector[e_phi]) * sin(State_vector[e_theta]);
                r2[y] = State_vector[e_r] * sin(State_vector[e_phi]) * sin(State_vector[e_theta]);
                r2[z] = State_vector[e_r] * cos(State_vector[e_theta]);

                double photon_tangent[3] = { r1[x] - r2[x], r1[y] - r2[y], r1[z] - r2[z] };
                double photon_LOS_parameter = -dot_product(r1, r1) / dot_product(r1, photon_tangent);
                double obs_plane_intersection[3] = { r1[x] + photon_LOS_parameter * photon_tangent[x],
                                                     r1[y] + photon_LOS_parameter * photon_tangent[y],
                                                     r1[z] + photon_LOS_parameter * photon_tangent[z] };

                double Image_coordiantes[3]{}; // Temporary array to store the rotated obs_plane_intersection vector
                Rorate_to_obs_plane(theta_obs, phi_obs, obs_plane_intersection, Image_coordiantes);

                Ray_results.Image_Coords[x] = -Image_coordiantes[0];
                Ray_results.Image_Coords[y] =  Image_coordiantes[2];

            }

            // Novikov-Thorne Model Evaluation

            bool inside_NT_disk = State_vector[e_r] * State_vector[e_r] > NT_Model.get_r_in()  * NT_Model.get_r_in() &&
                                  State_vector[e_r] * State_vector[e_r] < NT_Model.get_r_out() * NT_Model.get_r_out() &&
                                  crossed_equatior(State_vector, Old_state);

            if (inside_NT_disk) {

                Image_Order[Novikov_Thorne] = n_equator_crossings;

                Ray_results.Redshift_NT[Image_Order[Novikov_Thorne]]             = NT_Model.Redshift(J, State_vector, r_obs, theta_obs, Spacetimes);
                Ray_results.Flux_NT[Image_Order[Novikov_Thorne]]                 = NT_Model.get_flux(State_vector[e_r], Spacetimes);
                Ray_results.Source_Coords[e_r][Image_Order[Novikov_Thorne]]      = State_vector[e_r];
                Ray_results.Source_Coords[e_phi][Image_Order[Novikov_Thorne]]    = State_vector[e_phi];
                Ray_results.Three_Momentum[e_r][Image_Order[Novikov_Thorne]]     = State_vector[e_p_r];
                Ray_results.Three_Momentum[e_theta][Image_Order[Novikov_Thorne]] = State_vector[e_p_theta];

            }

            if (crossed_equatior(State_vector, Old_state)) {

                if (n_equator_crossings < ORDER_NUM - 1) {

                    n_equator_crossings += 1;

                }
            }

            for (int vector_indexer = 0; vector_indexer <= e_State_Number - 1; vector_indexer += 1) {

                Old_state[vector_indexer] = State_vector[vector_indexer];

            }

            // Evaluate logical flags for terminating the integration

            if (Spacetimes[e_metric]->terminate_integration(State_vector, Derivatives)) {

                Ray_results.Intensity[direct] = State_vector[e_Intensity];
                Ray_results.Optical_Depth     = State_vector[e_Optical_Depth];

                write_to_file(Ray_results, data, momentum_data);

                set_pixel_color(State_vector[e_Intensity], texture_indexer);
                texture_indexer += 3;

                integration_count = 0;

                break;

            }

            integration_count += 1;

        }

    }

}
