#pragma once
#define _USE_MATH_DEFINES
#include "IO_files.h"
#include "Enumerations.h"
#include "Constants.h"
#include "Spacetimes.h"

#include "Page_Thorne_model.h"
#include "Emission_Models.h"
#include "Disk_Models.h"

#include "General_GR_functions.h"
#include "General_math_functions.h"

#include "Rendering_Engine.h"
#include "Radiative_Transfer.h"

#include "Lensing.h"

#include <iostream>
#include <complex>

void static log_ray_path(double State_Vector[], Results_type* s_Ray_Results, Step_controller Controller, Initial_conditions_type* p_Init_Conditions){

    int& log_offset = s_Ray_Results->Ray_log_struct.Log_offset;
    double& R_throat = p_Init_Conditions->Metric_parameters.R_throat;

    for (int index = 0; index <= e_State_Number - 1; index++) {

        s_Ray_Results->Ray_log_struct.Ray_path_log[index + log_offset * e_State_Number] = State_Vector[index];

        // The wormhole metric works with a "global" radial coordinate, that goes negative on the other side of the throat.
        // The emission model can't work with this coordinate, so I log the normal spherical radial coordinate instead. 

        if (Wormhole == p_Init_Conditions->Metric_parameters.e_Spacetime && e_r == index) {

            s_Ray_Results->Ray_log_struct.Ray_path_log[e_r + log_offset * e_State_Number] = sqrt(State_Vector[e_r] * State_Vector[e_r] + R_throat * R_throat);

        }

    }

    s_Ray_Results->Ray_log_struct.Ray_path_log[e_step + log_offset * e_State_Number] = Controller.step;

}

void static log_ray_emission(double Stokes_Vector[e_Stokes_param_num], double Optical_depth, Results_type* p_Ray_Results, int log_index) {
    
    for (int stokes_idx = 0; stokes_idx <= e_Stokes_param_num - 1; stokes_idx++) {

        p_Ray_Results->Ray_log_struct.Ray_emission_log[stokes_idx][0 + 2 * log_index] = Stokes_Vector[stokes_idx] * p_Ray_Results->Intensity_scale;
        p_Ray_Results->Ray_log_struct.Ray_emission_log[stokes_idx][1 + 2 * log_index] = Optical_depth;
    }

}

void static Evaluate_Equatorial_Disk(const Simulation_Context_type* const p_Sim_Context,
                                     Results_type* const s_Ray_results,
                                     const double* const State_vector,
                                     const double* const Old_state, 
                                     int N_theta_turning_points) {

    /* ------------ The number of components is e_State_Number - 1 because we do not include the integration step. */
    double Crossing_State[e_State_Number - 1]{};
    double& R_throat = p_Sim_Context->p_Init_Conditions->Metric_parameters.R_throat;

    if (interpolate_crossing(State_vector, Old_state, Crossing_State)) {

        int Image_Order = compute_image_order(N_theta_turning_points, p_Sim_Context->p_Init_Conditions);

        if (Wormhole == p_Sim_Context->p_Init_Conditions->Metric_parameters.e_Spacetime) {

            // The wormhole metric uses the global coordinate ell^2 = r^2 - r_throat^2
            // Here I convert back to the r coordinate for the NT model evaluation

            Crossing_State[e_r] = sqrt(Crossing_State[e_r] * Crossing_State[e_r] - R_throat * R_throat);

        }

        double& r_in = p_Sim_Context->p_Init_Conditions->Disk_params.Page_Thorne_params.r_in;
        double& r_out = p_Sim_Context->p_Init_Conditions->Disk_params.Page_Thorne_params.r_out;

        if (Crossing_State[e_r] < r_out && Crossing_State[e_r] > r_in){

            s_Ray_results->Redshift_PT[Image_Order] = p_Sim_Context->p_PT_model->Redshift(Crossing_State, p_Sim_Context->p_Observer);
            s_Ray_results->Flux_PT[Image_Order]     = p_Sim_Context->p_PT_model->get_flux(Crossing_State);

        }

        s_Ray_results->Source_Coords[e_r][Image_Order]   = Crossing_State[e_r];
        s_Ray_results->Source_Coords[e_phi][Image_Order] = Crossing_State[e_phi];

        s_Ray_results->Photon_Momentum[e_r][Image_Order]     = Crossing_State[e_p_r];
        s_Ray_results->Photon_Momentum[e_theta][Image_Order] = Crossing_State[e_p_theta];
    }

}

void static Seperate_Image_into_orders(int const Max_theta_turning_points, 
                                       int const N_theta_turning_points,
                                       Initial_conditions_type* const s_Initial_Conditions, 
                                       Results_type* const p_Ray_results,
                                       double const Intensity[e_Stokes_param_num]) {

    int Image_Order = compute_image_order(Max_theta_turning_points - N_theta_turning_points, s_Initial_Conditions);
    
    for (int stokes_idx = 0; stokes_idx <= e_Stokes_param_num - 1; stokes_idx++) {
    
        p_Ray_results->Intensity[Image_Order][stokes_idx] = Intensity[stokes_idx] * p_Ray_results->Intensity_scale;
    
        for (int order_scan = Image_Order; order_scan < compute_image_order(Max_theta_turning_points, s_Initial_Conditions); order_scan++) {
    
            p_Ray_results->Intensity[Image_Order][stokes_idx] -= p_Ray_results->Intensity[order_scan + 1][stokes_idx];
    
        }
    
    }

}

void static Propagate_Stokes_vector(Radiative_Transfer_Integrator Integrator,
                                    const Transfer_functions_type Transfer_functions,
                                    double const step, 
                                    double* const Intensity) {

    switch (Integrator) {

    case Analytic:

        Analytic_Radiative_Transfer(const_cast<double*>(Transfer_functions.Emission_functions), 
                                    const_cast<double*>(Transfer_functions.Absorbtion_functions), 
                                    const_cast<double*>(Transfer_functions.Faradey_functions), 
                                    step, 
                                    Intensity);

        break;

    case Implicit_Trapezoid:

        Implicit_Trapezoid_Radiative_Transfer(const_cast<double*>(Transfer_functions.Emission_functions), 
                                              const_cast<double*>(Transfer_functions.Absorbtion_functions), 
                                              const_cast<double*>(Transfer_functions.Faradey_functions), 
                                              step, 
                                              Intensity);

        break;

    default:

        std::cout << "Integration method not supported for the radiative transfer equations!" << '\n';

        exit(ERROR);

    }

}

Return_Values static Construct_Stokes_Tetrad(double Tetrad[4][4],
                                             double inv_Tetrad[4][4],
                                             const Simulation_Context_type* const p_Sim_Context, 
                                             const double* const State_vector) {

    /* The reference of this implementation is the second RAPTOR paper: https://arxiv.org/pdf/2007.03045.pdf 
       Expressions 10 - 11. Note that their definition of 9d is wrong... the "g" in the denominator should 
       be omega, as defined in 9c. */


    /* ============================= Get the three 4-vectors from which we will construct the tetrad ============================= */

    // ---------------- The velocity vector -> this is chosen to be the local emitter's velocity when we are inside the emission medium, 
    // otherwise it is chosen to be the observer's velocity

    double Plasma_velocity_contravariant[4]{};
    double Plasma_velocity_covariant[4]{};

    // ---------------- The tetrad requires a spacelike vector -> this is chosen to be the local magnetic field 4-vector when we are inside the emission medium,
    // otherwise it is fixed to a constant vector.
    // NOTE: Check the conventions for the tetrad -> the choice of this constant vector can effectively rotate the observer's basis.

    double Spacelike_vector_contravariant[4] = {0, 0, 1, 0};
    double Spacelike_vector_covariant[4]{};

    // ---------------- We will need the state of all emission media
    Emission_medium_state_type s_Disk_state{};
    Emission_medium_state_type s_Hotspot_state{};

    const Metric_type s_Metric_photon = p_Sim_Context->p_Spacetime->get_metric(State_vector);

    // ---------------- Certain hotspot quantities are evaluated at the spot center (regardless of where the photon is), this metric variable accounts for 
    Metric_type s_Metric_emission = s_Metric_photon;

    // ---------------- Determine which part of the emission medium we are in -> this determines the local magnetic field. The hotspot and jet models are allowed to have their own 
    // local magneic field, while outside them the field is considered due to the accretion disk.

    p_Sim_Context->p_Emission_Model->get_plasma_velocity(p_Sim_Context->p_Init_Conditions->Hotspot_params.Position,
                                                         p_Sim_Context,
                                                         p_Sim_Context->p_Init_Conditions->Hotspot_params.Velocity_profile_type, 
                                                         p_Sim_Context->p_Init_Conditions->Hotspot_params.Radial_velocity_fraction,
                                                         s_Hotspot_state.Plasma_Velocity);

     p_Sim_Context->p_Emission_Model->get_plasma_velocity(State_vector,
                                                          p_Sim_Context,
                                                          p_Sim_Context->p_Init_Conditions->Disk_params.Velocity_profile_type,
                                                          p_Sim_Context->p_Init_Conditions->Disk_params.Radial_velocity_fraction,
                                                          s_Disk_state.Plasma_Velocity);

    const bool In_hotspot = p_Sim_Context->p_Emission_Model->p_Hotspot_Model->is_inside_hotspot(State_vector, 
                                                                                                s_Hotspot_state.Plasma_Velocity,
                                                                                               &s_Hotspot_state);

    const bool In_disk = p_Sim_Context->p_Emission_Model->p_Disk_Model->is_inside_disk(State_vector, 
                                                                                       p_Sim_Context->p_Emission_Model->p_Disk_Model->s_Disk_params.e_Disk_model, 
                                                                                      &s_Disk_state);

    if (In_hotspot && NULL != s_Hotspot_state.Plasma_Velocity) {

        /* We are inside the hotspot - we assume the dominant magnetic field here is whatever the local field of the spot is - 
           a.e. inisde the hotspot, the disk magnetic field is "screened" by the spot. */

        p_Sim_Context->p_Emission_Model->p_Hotspot_Model->get_density_and_temperature(State_vector, s_Hotspot_state.Plasma_Velocity, &s_Hotspot_state);

        s_Hotspot_state.Magnetization = p_Sim_Context->p_Init_Conditions->Hotspot_params.Magnetization;
        memcpy(s_Hotspot_state.Magnetic_fields.Mag_field_geometry_vector, p_Sim_Context->p_Init_Conditions->Hotspot_params.Mag_field_geometry, 3 * sizeof(double));

        s_Hotspot_state.Magnetic_fields.e_Mag_field_geometry          = p_Sim_Context->p_Init_Conditions->Hotspot_params.e_Mag_field_geometry;
        s_Hotspot_state.Magnetic_fields.e_Mag_field_magnitude_profile = p_Sim_Context->p_Init_Conditions->Hotspot_params.e_Mag_field_magnitude_profile;

        p_Sim_Context->p_Emission_Model->get_plasma_velocity(State_vector,
                                                             p_Sim_Context,
                                                             p_Sim_Context->p_Init_Conditions->Hotspot_params.Velocity_profile_type, 
                                                             p_Sim_Context->p_Init_Conditions->Hotspot_params.Radial_velocity_fraction,
                                                             s_Hotspot_state.Plasma_Velocity);

        //s_Metric_emission = p_Sim_Context->p_Spacetime->get_metric(p_Sim_Context->p_Init_Conditions->Hotspot_params.Position);
        p_Sim_Context->p_Emission_Model->get_magnetic_field(State_vector, &s_Metric_emission, &s_Hotspot_state);
        
        memcpy(Spacelike_vector_contravariant, s_Hotspot_state.Magnetic_fields.B_field_plasma_frame, 4 * sizeof(double));
        memcpy(Plasma_velocity_contravariant, s_Hotspot_state.Plasma_Velocity, 4 * sizeof(double));

    }
    else if ((!In_hotspot && In_disk) && NULL != s_Disk_state.Plasma_Velocity){

        /* We are outside the hotspot - we assume the dominant magnetic field here is due to the background accretion disk. */

        p_Sim_Context->p_Emission_Model->p_Disk_Model->get_density_and_temperature(State_vector, 
                                                                                   p_Sim_Context->p_Emission_Model->p_Disk_Model->s_Disk_params.e_Disk_model, 
                                                                                  &s_Disk_state);

        s_Disk_state.Magnetization = p_Sim_Context->p_Init_Conditions->Disk_params.Magnetization;
        memcpy(s_Disk_state.Magnetic_fields.Mag_field_geometry_vector, p_Sim_Context->p_Init_Conditions->Disk_params.Mag_field_geometry, 3 * sizeof(double));
        
        s_Disk_state.Magnetic_fields.e_Mag_field_geometry = p_Sim_Context->p_Init_Conditions->Disk_params.e_Mag_field_geometry;
        s_Disk_state.Magnetic_fields.e_Mag_field_magnitude_profile = p_Sim_Context->p_Init_Conditions->Disk_params.e_Mag_field_magnitude_profile;
        p_Sim_Context->p_Emission_Model->get_magnetic_field(State_vector, &s_Metric_emission, &s_Disk_state);
        
        memcpy(Spacelike_vector_contravariant, s_Disk_state.Magnetic_fields.B_field_plasma_frame, 4 * sizeof(double));
        memcpy(Plasma_velocity_contravariant, s_Disk_state.Plasma_Velocity, 4 * sizeof(double));
    }
    else if (NULL != s_Disk_state.Plasma_Velocity) {

        //memcpy(Spacelike_vector_contravariant, s_Hotspot_state.Magnetic_fields.B_field_plasma_frame, 4 * sizeof(double));
        memcpy(Plasma_velocity_contravariant, p_Sim_Context->p_Observer->get_fiducial_obs_velocity(&s_Metric_emission), 4 * sizeof(double));

    }

    const double Wave_Vector_covariant[4] = { State_vector[e_p_t], State_vector[e_p_r], State_vector[e_p_theta], State_vector[e_p_phi] };

    /* --------------------- Evaluate the inner products (expressions 9a - 9d), and compute the contravariant Wave-Vector --------------------- */

    double Wave_vec_dot_Plasma_vel{}, Plasma_vel_dot_Spacelike_vec{}, Spacelike_vec_norm_squared{}, Wave_vec_dot_Spacelike_vec{};

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        Wave_vec_dot_Plasma_vel    += Wave_Vector_covariant[left_idx] * Plasma_velocity_contravariant[left_idx];
        Wave_vec_dot_Spacelike_vec += Wave_Vector_covariant[left_idx] * Spacelike_vector_contravariant[left_idx];

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            Plasma_vel_dot_Spacelike_vec += s_Metric_photon.Metric[left_idx][right_idx] * Plasma_velocity_contravariant[left_idx] * Spacelike_vector_contravariant[right_idx];
            Spacelike_vec_norm_squared   += s_Metric_photon.Metric[left_idx][right_idx] * Spacelike_vector_contravariant[left_idx] * Spacelike_vector_contravariant[right_idx];

        }
        
    }

    const double metric_determinant = get_metric_det(s_Metric_photon.Metric);
    const double sqrt_determinant = sqrt(-metric_determinant);

    const double C_coeff = -Wave_vec_dot_Spacelike_vec / Wave_vec_dot_Plasma_vel - Plasma_vel_dot_Spacelike_vec;
    const double N_coeff = sqrt(Spacelike_vec_norm_squared + Plasma_vel_dot_Spacelike_vec * Plasma_vel_dot_Spacelike_vec - C_coeff * C_coeff);

    /* --- Perform checks on these coefficients, because for very low plasma densities they blow up --- */

    if (isnan(N_coeff) || isinf(1.0 / N_coeff) || isnan(C_coeff) || isinf(1.0 / C_coeff)) {

        return ERROR;

    }

    /* --- Raise / Lower indicies on the three main 4-vectors - this is needed for computing the final tetrad vector --- */

    double inv_Metric[4][4]{};
    invert_metric(inv_Metric, s_Metric_photon.Metric);

    double Wave_vector_contravariant[4]{};

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            Wave_vector_contravariant[left_idx] += inv_Metric[left_idx][right_idx] * Wave_Vector_covariant[right_idx];

            // These indecies are lowered with the s_Metric_emission metric, because they were calculated with it
            Spacelike_vector_covariant[left_idx] += s_Metric_photon.Metric[left_idx][right_idx] * Spacelike_vector_contravariant[right_idx];
            Plasma_velocity_covariant[left_idx]  += s_Metric_photon.Metric[left_idx][right_idx] * Plasma_velocity_contravariant[right_idx];

        }

    }

    /* --------------------- Compute the first two tetrad basis vectors --------------------- */

    // Wave_vec_dot_Plasma_vel is safe to divide by as it can never go to zero - Plasma_vel is never null.

    for (int index = 0; index <= 3; index++) {

        Tetrad[e_t][index]   = Plasma_velocity_contravariant[index];
        Tetrad[e_phi][index] = -Wave_vector_contravariant[index] / Wave_vec_dot_Plasma_vel - Plasma_velocity_contravariant[index];

    }

    /* -------- Compute the contravariant 4D permutation symbol (this I lifted straight from the RAPTOR code...) -------- */

    double Levi_Cevita_tensor[4][4][4][4]{};

    for (int i = 0; i <= 3; i++) {

        for (int j = 0; j <= 3; j++) {

            for (int k = 0; k <= 3; k++) {

                for (int l = 0; l <= 3; l++) {

                    Levi_Cevita_tensor[i][j][k][l] = -((i - j) * (i - k) * (i - l) * (j - k) * (j - l) * (k - l) / 12.) / sqrt_determinant;
                                                   
                }
            }
        }
    }

    /* --------------------- Compute the last two tetrad basis vectors --------------------- */
 
    for (int index = 0; index <= 3; index++) {

        Tetrad[e_theta][index] = (Spacelike_vector_contravariant[index] + Plasma_vel_dot_Spacelike_vec * Plasma_velocity_contravariant[index] - C_coeff * Tetrad[e_phi][index]) / N_coeff;

        for (int i = 0; i <= 3; i++) {

            for (int j = 0; j <= 3; j++) {

                for (int k = 0; k <= 3; k++) {

                    Tetrad[e_r][index] += -Levi_Cevita_tensor[index][i][j][k] *
                                                     Plasma_velocity_covariant[i] *
                                                     Wave_Vector_covariant[j] *
                                                     Spacelike_vector_covariant[k] / Wave_vec_dot_Plasma_vel / N_coeff;

                }

            }

        }

    }

    /* --------------------- Compute the inverse tetrad --------------------- */

    // Technically this is the inverse Minkowski metric (with upper indicies)
    const double Minkowski_Metric[4][4] = { {-1., 0., 0., 0.},
                                            { 0., 1., 0., 0.},
                                            { 0., 0., 1., 0.},
                                            { 0., 0., 0., 1.} };

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            for (int m = 0; m <= 3; m++) {

                for (int g = 0; g <= 3; g++) {

                    inv_Tetrad[left_idx][right_idx] += Minkowski_Metric[left_idx][m] * s_Metric_photon.Metric[right_idx][g] * Tetrad[m][g];

                    if (isnan(inv_Tetrad[left_idx][right_idx]) || isinf(inv_Tetrad[left_idx][right_idx]) ||
                        isnan(Tetrad[left_idx][right_idx])     || isinf(Tetrad[left_idx][right_idx])) {

                        return ERROR;

                    }
                        
                }
            }
        }
    }

    //double test[4][4]{};

    //for (int a = 0; a <= 3; a++) {

    //    for (int b = 0; b <= 3; b++) {

    //        for (int left_idx = 0; left_idx <= 3; left_idx++) {

    //            for (int right_idx = 0; right_idx <= 3; right_idx++) {

    //                test[a][b] += s_Metric_photon.Metric[right_idx][left_idx] * Tetrad[a][left_idx] * Tetrad[b][right_idx];


    //            }
    //        }
    //    }
    //}

    return OK;

}

void static Parallel_Transport_Polarization_Vector(double State_Vector[], 
                                                   Spacetime_Base_Class* const Spacetime, 
                                                   std::complex<double> Polarization_Vector[]) {

    Metric_type s_Metric        = Spacetime->get_metric(State_Vector);
    Metric_type s_dr_Metric     = Spacetime->get_dr_metric(State_Vector);
    Metric_type s_dtheta_Metric = Spacetime->get_dtheta_metric(State_Vector);

    double inv_Metric[4][4]{};
    invert_metric(inv_Metric, s_Metric.Metric);
    
    double Connection_Coefficients[4][4][4]{};

    get_connection_coefficients(s_Metric, s_dr_Metric, s_dtheta_Metric, Connection_Coefficients);

    std::complex<double> Polarization_Vector_Derivative[4]{};

  /* ========================== Construct the full CONTRVARIANT photon wave vector ========================== */

    double photon_wave_vector[4]{};

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            photon_wave_vector[left_idx] += inv_Metric[left_idx][right_idx] * State_Vector[right_idx + e_p_t];

        }

    }

  /* ========================== Compute the derivative of the polarization vector from the parallel transport ========================== */

    for (int derivative_index = 0; derivative_index <= e_Stokes_param_num - 1; derivative_index++) {

        for (int polarization_index = 0; polarization_index <= e_Stokes_param_num - 1; polarization_index++) {

            for (int wave_vector_index = 0; wave_vector_index <= e_Stokes_param_num - 1; wave_vector_index++) {

                Polarization_Vector_Derivative[derivative_index] = -Connection_Coefficients[derivative_index][wave_vector_index][polarization_index] * photon_wave_vector[wave_vector_index] * Polarization_Vector[polarization_index];

            }

        }

    }

    for (int index = 0; index <= 3; index++) {

        Polarization_Vector[index] += Polarization_Vector_Derivative[index] * State_Vector[e_step];

    }

    double Norm{};

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            Norm += (s_Metric.Metric[left_idx][right_idx] * Polarization_Vector[left_idx] * std::conj(Polarization_Vector[right_idx])).real();

        }

    }

    for (int index = 0; index <= 3; index++) {

        Polarization_Vector[index] /= Norm;

    }

}

void static Map_Polarization_Vector_to_Stokes(const double inv_Stokes_Tetrad[4][4],
                                              std::complex<double> Coord_Basis_Pol_vec[e_Stokes_param_num],
                                              double Stokes_Vector[e_Stokes_param_num]) {

            std::complex<double> Stokes_Basis_Pol_vec[4]{};

            for (int stokes_idx = 0; stokes_idx <= 3; stokes_idx++) {

                Stokes_Basis_Pol_vec[stokes_idx] = (0.0, 0.0);

                for (int coord_idx = 0; coord_idx <= 3; coord_idx++) {

                    Stokes_Basis_Pol_vec[stokes_idx] += inv_Stokes_Tetrad[stokes_idx][coord_idx] * Coord_Basis_Pol_vec[coord_idx];


                }

            }

            double Polarized_Intensity_before = sqrt(Stokes_Vector[Q] * Stokes_Vector[Q] +
                                                     Stokes_Vector[U] * Stokes_Vector[U] +
                                                     Stokes_Vector[V] * Stokes_Vector[V]);

           
            Stokes_Vector[Q] =  Polarized_Intensity_before * (Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[1]) -
                                                              Stokes_Basis_Pol_vec[2] * std::conj(Stokes_Basis_Pol_vec[2])).real();

            Stokes_Vector[U] =  Polarized_Intensity_before * (Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[2]) +
                                                              Stokes_Basis_Pol_vec[2] * std::conj(Stokes_Basis_Pol_vec[1])).real();

            Stokes_Vector[V] = -Polarized_Intensity_before * (complex_i * (Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[2]) -
                                                                           Stokes_Basis_Pol_vec[2] * std::conj(Stokes_Basis_Pol_vec[1]))).real();
            
            /* The numerics seem to introduce a surprisingly large error in the norm of this vector, which results in a polarization fraction > 1.
               Normalizing here the Polarization vector in the Stokes basis seems to resolve the issue. 

               TODO: Check if this is really numerics or an analytical error! */

            double Polarized_Intensity_after = sqrt(Stokes_Vector[Q] * Stokes_Vector[Q] +
                                                    Stokes_Vector[U] * Stokes_Vector[U] +
                                                    Stokes_Vector[V] * Stokes_Vector[V]);


            if (!isinf(1.0 / Polarized_Intensity_after)) {

                Stokes_Vector[Q] *= Polarized_Intensity_before / Polarized_Intensity_after;
                Stokes_Vector[U] *= Polarized_Intensity_before / Polarized_Intensity_after;
                Stokes_Vector[V] *= Polarized_Intensity_before / Polarized_Intensity_after;

            }

}

void static Map_Stokes_to_Polarization_Vector(const double Stokes_Vector[e_Stokes_param_num],
                                              const double Stokes_Tetrad[4][4],
                                              std::complex<double> Coord_Basis_Pol_vec[e_Stokes_param_num]) {

    for (int index = 0; index <= 3; index++) {

        Coord_Basis_Pol_vec[index] = (0, 0);

    }

    std::complex<double> Stokes_Basis_Pol_vec[4];

    double Polarized_Intensity = sqrt(Stokes_Vector[Q] * Stokes_Vector[Q] +
                                      Stokes_Vector[U] * Stokes_Vector[U] +
                                      Stokes_Vector[V] * Stokes_Vector[V]);

    if (!isinf(1.0 / Polarized_Intensity) && fabs(Stokes_Vector[Q] / Polarized_Intensity) < 1) {

        Stokes_Basis_Pol_vec[1] = sqrt((1 + Stokes_Vector[Q] / Polarized_Intensity) / 2);

    }

    Stokes_Basis_Pol_vec[2] = 1.0;

    if (!isinf(1.0 / std::norm(Stokes_Basis_Pol_vec[1] * Polarized_Intensity))) {

        Stokes_Basis_Pol_vec[2] = (Stokes_Vector[U] - complex_i * Stokes_Vector[V]) / (2.0 * Stokes_Basis_Pol_vec[1] * Polarized_Intensity);

    }

    for (int coord_idx = 0; coord_idx <= 3; coord_idx++) {

        for (int stokes_idx = 0; stokes_idx <= 3; stokes_idx++) {

            Coord_Basis_Pol_vec[coord_idx] += Stokes_Tetrad[stokes_idx][coord_idx] * Stokes_Basis_Pol_vec[stokes_idx];

        }

    }

}

void static Propagate_forward_emission(const Simulation_Context_type* const p_Sim_Context, 
                                       Results_type* const p_Ray_results,
                                       int* const N_theta_turning_points) {

    int const Max_theta_turning_points = *N_theta_turning_points;
    *N_theta_turning_points = 0;

    double Stokes_Vector[e_Stokes_param_num]{};

    std::complex<double> Coord_Basis_Pol_vec[4]{};

    // TODO: Propagate this aswell
    double Optical_Depth{};

    for (int log_index = p_Ray_results->Ray_log_struct.Log_length; log_index > 0; log_index--) {

        /* =============== Pick out the ray position / momenta from the Log, at the given log index =============== */

        double* Logged_ray_path = &(p_Ray_results->Ray_log_struct.Ray_path_log[log_index * e_State_Number]);

        /* ====================================== Parallel transport the polarization vector ====================================== */

        double Tetrad[4][4]{};
        double inv_Tetrad[4][4]{};

        if (p_Sim_Context->p_Init_Conditions->Observer_params.include_polarization) {

            /* If a Stokes basis cannot be contstructed, skip the parallel transport.
               NOTE: This should never happen. */

            if (OK == Construct_Stokes_Tetrad(Tetrad, inv_Tetrad, p_Sim_Context, Logged_ray_path)) {

                Map_Stokes_to_Polarization_Vector(std::as_const(Stokes_Vector), std::as_const(Tetrad), Coord_Basis_Pol_vec);

                Parallel_Transport_Polarization_Vector(Logged_ray_path, p_Sim_Context->p_Spacetime, Coord_Basis_Pol_vec);

                Map_Polarization_Vector_to_Stokes(std::as_const(inv_Tetrad), Coord_Basis_Pol_vec, Stokes_Vector);

            }
        }

        /* ======================================================================================================================== */

        /* ====================================== Propagate the radiative transfer equations ====================================== */

        Transfer_functions_type total_Transfer_functions{};

        /* Loop trough each emission medium (Disk, Hotspot, Jet and so on) and sum their respective transfer functions */
        for (int emission_medium = Disk; emission_medium <= Hotspot; emission_medium++) {

            Transfer_functions_type temp_Transfer_functions{};

            p_Sim_Context->p_Emission_Model->get_radiative_transfer_functions(Logged_ray_path,
                                                                              p_Sim_Context,
                                                                              static_cast<Emission_medium_enums>(emission_medium),
                                                                             &temp_Transfer_functions);

            add_vectors(temp_Transfer_functions.Emission_functions, total_Transfer_functions.Emission_functions, e_Stokes_param_num, total_Transfer_functions.Emission_functions);
            add_vectors(temp_Transfer_functions.Faradey_functions, total_Transfer_functions.Faradey_functions, e_Stokes_param_num, total_Transfer_functions.Faradey_functions);
            add_vectors(temp_Transfer_functions.Absorbtion_functions, total_Transfer_functions.Absorbtion_functions, e_Stokes_param_num, total_Transfer_functions.Absorbtion_functions);

        }

        Propagate_Stokes_vector(Implicit_Trapezoid, total_Transfer_functions, Logged_ray_path[e_step], Stokes_Vector);

        /* ======================================================================================================================== */

        log_ray_emission(Stokes_Vector, Optical_Depth, p_Ray_results, log_index);

        *N_theta_turning_points += Check_for_theta_turning_point(Logged_ray_path, Logged_ray_path - e_State_Number);
        Seperate_Image_into_orders(Max_theta_turning_points, *N_theta_turning_points, p_Sim_Context->p_Init_Conditions, p_Ray_results, Stokes_Vector);

    }

}

void Propagate_ray(const Simulation_Context_type* const p_Sim_Context, Results_type* const p_Ray_results) {

    // Initialize the State Vectors
    double State_Vector[e_State_Number]{};
    double Old_State_Vector[e_State_Number]{};

    State_Vector[e_t]       = p_Sim_Context->p_Init_Conditions->Observer_params.init_time;
    State_Vector[e_r]       = p_Sim_Context->p_Init_Conditions->Observer_params.distance;
    State_Vector[e_theta]   = p_Sim_Context->p_Init_Conditions->Observer_params.inclination;
    State_Vector[e_phi]     = p_Sim_Context->p_Init_Conditions->Observer_params.azimuth;
    State_Vector[e_p_phi]   = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_phi];
    State_Vector[e_p_theta] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_theta];
    State_Vector[e_p_r]     = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_r];
    State_Vector[e_p_t]     = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_t];

    // Set the Old State Vector to the Initial State Vector
    memcpy(Old_State_Vector, State_Vector, e_State_Number * sizeof(double));

    for (int Image_order = e_direct; Image_order < e_order_number; Image_order += 1) {

        p_Ray_results->Photon_Momentum[e_phi][Image_order] = State_Vector[e_p_phi];
        p_Ray_results->Photon_Momentum[e_t][Image_order]   = State_Vector[e_p_t];

    }

    p_Ray_results->Metric_parameters = p_Sim_Context->p_Init_Conditions->Metric_parameters;

    // Initialize counters for the Number Of Integration Steps and the Number Of Turning points of the Polar Coordinate
    int integration_count{}, N_theta_turning_points{};

    // Calculate the image coordinates from the initial conditions
    get_image_coordinates(p_Sim_Context->p_Init_Conditions, p_Ray_results->Image_Coords);

    Step_controller controller(p_Sim_Context->p_Init_Conditions->Integrator_params);

    p_Ray_results->Ray_log_struct.Log_offset = 0;
    log_ray_path(State_Vector, p_Ray_results, controller, p_Sim_Context->p_Init_Conditions);

    while (!controller.integration_complete && integration_count <= controller.Parameters.Max_integration_count) {

        RK45(State_Vector, &controller, p_Sim_Context);

        if (controller.continue_integration) {

            integration_count += 1;
            p_Ray_results->Ray_log_struct.Log_offset = integration_count;

            log_ray_path(State_Vector, p_Ray_results, controller, p_Sim_Context->p_Init_Conditions);

            N_theta_turning_points += Check_for_theta_turning_point(State_Vector, Old_State_Vector);

            /* ======================================== Evaluate the thin disk models ======================================== */

            if (e_Page_Thorne == p_Sim_Context->p_Init_Conditions->Disk_params.e_Disk_model) {

                Evaluate_Equatorial_Disk(p_Sim_Context, p_Ray_results, State_Vector, Old_State_Vector, N_theta_turning_points);

            }

            /* ============================================================================================================== */

            memcpy(Old_State_Vector, State_Vector, e_State_Number * sizeof(double));

        }

    }

    if (integration_count >= controller.Parameters.Max_integration_count) { std::cout << "Max iterations reached!" << '\n'; }

    p_Ray_results->Ray_log_struct.Log_length = p_Ray_results->Ray_log_struct.Log_offset;

    /* =========== Integrate the radiative transfer equations forward along the ray for the RIAF models =========== */

    if (e_Page_Thorne != p_Sim_Context->p_Init_Conditions->Disk_params.e_Disk_model) {

        Propagate_forward_emission(p_Sim_Context, p_Ray_results, &N_theta_turning_points);

    }

    /* ============================================================================================================ */

}