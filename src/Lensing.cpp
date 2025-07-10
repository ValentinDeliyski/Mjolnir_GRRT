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

void static log_ray_path(double* State_Vector, Results_type* s_Ray_Results, Initial_conditions_type* p_Init_Conditions){

    const int& log_offset = s_Ray_Results->Ray_log_struct.Log_offset;
    const double& R_throat = p_Init_Conditions->Metric_parameters.R_throat;

    memcpy(&s_Ray_Results->Ray_log_struct.Ray_path_log[log_offset * e_Full_state_size], State_Vector, e_Full_state_size * sizeof(double));

    // The wormhole metric works with a "global" radial coordinate, that goes negative on the other side of the throat.
    // The emission model can't work with this coordinate, so I log the normal spherical radial coordinate instead. 
    if (Wormhole == p_Init_Conditions->Metric_parameters.e_Spacetime) {

        s_Ray_Results->Ray_log_struct.Ray_path_log[e_r + log_offset * e_Full_state_size] = sqrt(State_Vector[e_r] * State_Vector[e_r] + R_throat * R_throat);

    }

}

void static log_ray_emission(double Stokes_Vector[e_Stokes_param_num], double Optical_depth, Results_type* p_Ray_Results, int log_index, int Sim_mode) {
    
    for (int stokes_idx = 0; stokes_idx <= e_Stokes_param_num - 1; stokes_idx++) {

        p_Ray_Results->Ray_log_struct.Ray_emission_log[stokes_idx][0 + 2 * log_index] = Stokes_Vector[stokes_idx] ;
        p_Ray_Results->Ray_log_struct.Ray_emission_log[stokes_idx][1 + 2 * log_index] = Optical_depth;

        if (3 != Sim_mode) {

            p_Ray_Results->Ray_log_struct.Ray_emission_log[stokes_idx][0 + 2 * log_index] *= p_Ray_Results->Intensity_scale;

        }

    }

}

void static Evaluate_Equatorial_Disk(const Simulation_Context_type* const p_Sim_Context,
                                     Results_type* const s_Ray_results,
                                     const double* const State_vector,
                                     const double* const Old_state, 
                                     const int N_theta_turning_points) {

    /* ------------ The number of components is e_State_Number - 1 because we do not include the integration step. */
    double Crossing_State[e_Dynamic_state_size - 1]{};
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

void static Seperate_Image_into_orders(const int Max_order,
                                       const double Stokes_Vector_offset[e_Stokes_param_num][e_order_number],
                                       Results_type* const p_Ray_results) {

    for (int stokes_idx = 0; stokes_idx <= e_Stokes_param_num - 1; stokes_idx++) {
    
        for (int order_scan = 0; order_scan <= Max_order; order_scan++) {

            p_Ray_results->Intensity[order_scan][stokes_idx] = Stokes_Vector_offset[stokes_idx][order_scan];

            if (order_scan + 1 <= Max_order) {

                p_Ray_results->Intensity[order_scan][stokes_idx] -= Stokes_Vector_offset[stokes_idx][order_scan + 1];

            }

            p_Ray_results->Intensity[order_scan][stokes_idx] *= p_Ray_results->Intensity_scale;

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
                                             const double* const State_vector,
                                             double* const Trial_vector,
                                             const bool Overwride) {

    memset(Tetrad, 0, 16 * sizeof(double));
    memset(inv_Tetrad, 0, 16 * sizeof(double));

    /* The reference of this implementation is the second RAPTOR paper: https://arxiv.org/pdf/2007.03045.pdf 
       Expressions 10 - 11. Note that their definition of 9d is wrong... the "g" in the denominator should 
       be omega, as defined in 9c. */

    Metric_type s_Metric_photon = p_Sim_Context->p_Spacetime->get_metric(State_vector);

    Return_Values Plasma_velocity_OK{};

    /* ============================= Get the three 4-vectors from which we will construct the tetrad ============================= */

    // ---------------- The velocity vector -> this is chosen to be the local emitter's velocity when we are inside the emission medium, 
    // otherwise it is chosen to be the observer's velocity. I use the theta dependant profile for the observer velocity, because it tends to be well defined below ISCO.

    double Plasma_velocity_contravariant[4]{};
    double Plasma_velocity_covariant[4]{};

    // I define these as variables for the sake of code readability.
    const Velocity_enums Obs_velocity_profile = e_Theta_dependant;
    const double Obs_radial_velocity_fraction = 0.0;

    Plasma_velocity_OK = p_Sim_Context->p_Emission_Model->get_plasma_velocity(State_vector,
                                                                              p_Sim_Context,
                                                                              Obs_velocity_profile,
                                                                              Obs_radial_velocity_fraction,
                                                                              Plasma_velocity_contravariant);

    if (ERROR == Plasma_velocity_OK) { return ERROR; }

    // ---------------- The tetrad requires a spacelike vector -> this is chosen to be the local magnetic field 4-vector when we are inside the emission medium,
    // otherwise it is fixed to a constant vector.
    // NOTE: Check the conventions for the tetrad -> the choice of this constant vector can effectively rotate the observer's basis.

    double Spacelike_vector_contravariant[4]{};
    memcpy(Spacelike_vector_contravariant, Trial_vector, 4 * sizeof(double));

    double Spacelike_vector_covariant[4]{};

    const double Wave_Vector_covariant[4] = { State_vector[e_p_t], State_vector[e_p_r], State_vector[e_p_theta], State_vector[e_p_phi] };

    // ---------------- We will need the state of all emission media
    Emission_medium_state_type s_Disk_state{};
    Emission_medium_state_type s_Hotspot_state{};

    // ---------------- Certain hotspot quantities are evaluated at the spot center (regardless of where the photon is), this metric variable accounts for that
    Metric_type s_Metric_emission = s_Metric_photon;

    // ---------------- Determine which part of the emission medium we are in -> this determines the local magnetic field. The hotspot and jet models are allowed to have their own 
    // local magneic field, while outside them the field is considered due to the accretion disk.

    p_Sim_Context->p_Emission_Model->get_plasma_velocity(p_Sim_Context->p_Init_Conditions->Hotspot_params.Position,
                                                         p_Sim_Context,
                                                         p_Sim_Context->p_Init_Conditions->Hotspot_params.Velocity_profile_type,
                                                         p_Sim_Context->p_Init_Conditions->Hotspot_params.Radial_velocity_fraction,
                                                         s_Hotspot_state.Plasma_Velocity);

    const bool In_hotspot = p_Sim_Context->p_Emission_Model->p_Hotspot_Model->is_inside_hotspot(State_vector,
                                                                                                &s_Hotspot_state);

    const bool In_disk = p_Sim_Context->p_Emission_Model->p_Disk_Model->is_inside_disk(State_vector,
                                                                                       p_Sim_Context->p_Emission_Model->p_Disk_Model->s_Disk_params.e_Disk_model,
                                                                                       &s_Disk_state);
    if (!Overwride) {

        if (In_hotspot) {

            /* We are inside the hotspot - we assume the dominant magnetic field here is whatever the local field of the spot is -
               a.e. inisde the hotspot, the disk magnetic field is "screened" by the spot. */

            p_Sim_Context->p_Emission_Model->p_Hotspot_Model->get_density_and_temperature(State_vector, &s_Hotspot_state);

            s_Hotspot_state.Magnetization = p_Sim_Context->p_Init_Conditions->Hotspot_params.Magnetization;
            memcpy(s_Hotspot_state.Magnetic_fields.Mag_field_geometry_vector, p_Sim_Context->p_Init_Conditions->Hotspot_params.Mag_field_geometry, 3 * sizeof(double));

            s_Hotspot_state.Magnetic_fields.e_Mag_field_geometry = p_Sim_Context->p_Init_Conditions->Hotspot_params.e_Mag_field_geometry;
            s_Hotspot_state.Magnetic_fields.e_Mag_field_magnitude_profile = p_Sim_Context->p_Init_Conditions->Hotspot_params.e_Mag_field_magnitude_profile;

            Plasma_velocity_OK = p_Sim_Context->p_Emission_Model->get_plasma_velocity(State_vector,
                                                                                      p_Sim_Context,
                                                                                      p_Sim_Context->p_Init_Conditions->Hotspot_params.Velocity_profile_type,
                                                                                      p_Sim_Context->p_Init_Conditions->Hotspot_params.Radial_velocity_fraction,
                                                                                      s_Hotspot_state.Plasma_Velocity);

            if (OK == Plasma_velocity_OK) {

                //s_Metric_emission = p_Sim_Context->p_Spacetime->get_metric(p_Sim_Context->p_Init_Conditions->Hotspot_params.Position);
                p_Sim_Context->p_Emission_Model->get_magnetic_field(State_vector, &s_Metric_emission, &s_Hotspot_state);

                memcpy(Spacelike_vector_contravariant, s_Hotspot_state.Magnetic_fields.B_field_plasma_frame, 4 * sizeof(double));
                memcpy(Plasma_velocity_contravariant, s_Hotspot_state.Plasma_Velocity, 4 * sizeof(double));

            }

        }
        else if (!In_hotspot && In_disk) {

            /* We are outside the hotspot - we assume the dominant magnetic field here is due to the background accretion disk. */

            p_Sim_Context->p_Emission_Model->p_Disk_Model->get_density_and_temperature(State_vector,
                                                                                       p_Sim_Context->p_Emission_Model->p_Disk_Model->s_Disk_params.e_Disk_model,
                                                                                       &s_Disk_state);

            s_Disk_state.Magnetization = p_Sim_Context->p_Init_Conditions->Disk_params.Magnetization;
            memcpy(s_Disk_state.Magnetic_fields.Mag_field_geometry_vector, p_Sim_Context->p_Init_Conditions->Disk_params.Mag_field_geometry, 3 * sizeof(double));

            s_Disk_state.Magnetic_fields.e_Mag_field_geometry = p_Sim_Context->p_Init_Conditions->Disk_params.e_Mag_field_geometry;
            s_Disk_state.Magnetic_fields.e_Mag_field_magnitude_profile = p_Sim_Context->p_Init_Conditions->Disk_params.e_Mag_field_magnitude_profile;
            p_Sim_Context->p_Emission_Model->get_magnetic_field(State_vector, &s_Metric_emission, &s_Disk_state);

            Plasma_velocity_OK = p_Sim_Context->p_Emission_Model->get_plasma_velocity(State_vector,
                                                                                      p_Sim_Context,
                                                                                      p_Sim_Context->p_Init_Conditions->Disk_params.Velocity_profile_type,
                                                                                      p_Sim_Context->p_Init_Conditions->Disk_params.Radial_velocity_fraction,
                                                                                      s_Disk_state.Plasma_Velocity);

            if (OK != Plasma_velocity_OK) { return ERROR; }

            memcpy(Spacelike_vector_contravariant, s_Disk_state.Magnetic_fields.B_field_plasma_frame, 4 * sizeof(double));
            memcpy(Plasma_velocity_contravariant, s_Disk_state.Plasma_Velocity, 4 * sizeof(double));

        }

    }

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
    const double N_coeff = sqrt((Spacelike_vec_norm_squared + Plasma_vel_dot_Spacelike_vec * Plasma_vel_dot_Spacelike_vec - C_coeff * C_coeff));

    /* --- Perform checks on these coefficients, because for very low plasma densities they blow up --- */

    if (isnan(N_coeff) || isinf(1.0 / N_coeff) || isnan(C_coeff)) { return ERROR; }

    /* --- Raise / Lower indicies on the three main 4-vectors - this is needed for computing the final tetrad vector --- */

    double inv_Metric[4][4]{};
    invert_metric(inv_Metric, s_Metric_photon.Metric);

    double Wave_vector_contravariant[4]{};

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            Wave_vector_contravariant[left_idx] += inv_Metric[left_idx][right_idx] * Wave_Vector_covariant[right_idx];

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

                        
                }
            }

            if (isnan(inv_Tetrad[left_idx][right_idx]) || isinf(inv_Tetrad[left_idx][right_idx]) ||
                isnan(Tetrad[left_idx][right_idx]) || isinf(Tetrad[left_idx][right_idx])) {

                return ERROR;

            }

        }
    }

    double test[4][4]{};

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            for (int m = 0; m <= 3; m++) {

                    test[left_idx][right_idx] += inv_Tetrad[left_idx][m] * Tetrad[right_idx][m];

            }

        }
    }

    return OK;

}

void static Parallel_Transport_RHS(const double* const State_Vector, 
                                   std::complex<double>* Polarization_Vector,
                                   const Spacetime_Base_Class* const p_Spacetime, 
                                   std::complex<double>* Polarization_Vector_Derivative) {

    Metric_type s_Metric        = p_Spacetime->get_metric(State_Vector);
    Metric_type s_dr_Metric     = p_Spacetime->get_dr_metric(State_Vector);
    Metric_type s_dtheta_Metric = p_Spacetime->get_dtheta_metric(State_Vector);

    double inv_Metric[4][4]{};
    invert_metric(inv_Metric, s_Metric.Metric);

    double Connection_Coefficients[4][4][4]{};

    get_connection_coefficients(s_Metric, s_dr_Metric, s_dtheta_Metric, Connection_Coefficients);

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

                Polarization_Vector_Derivative[derivative_index] += -Connection_Coefficients[derivative_index][wave_vector_index][polarization_index] * photon_wave_vector[wave_vector_index] * Polarization_Vector[polarization_index];

            }

        }

    }

}

void static Parallel_Transport_Polarization_Vector(const double* const State_Vector, 
                                                   Spacetime_Base_Class* const p_Spacetime, 
                                                   std::complex<double>* Polarization_Vector) {

    std::complex<double> RHS1[4]{};
    std::complex<double> RHS2[4]{};
    std::complex<double> RHS3[4]{};
    std::complex<double> RHS4[4]{};

    double EOM1[e_Dynamic_state_size]{};
    double EOM2[e_Dynamic_state_size]{};
    double EOM3[e_Dynamic_state_size]{};
    double EOM4[e_Dynamic_state_size]{};

    std::complex<double> Temp_pol_vec[4]{};
    double Temp_State_Vector[e_Dynamic_state_size]{};

    /* ================================== 1 ========================================== */

    Parallel_Transport_RHS(State_Vector, Polarization_Vector, p_Spacetime, RHS1);
    p_Spacetime->get_EOM(State_Vector, EOM1);

    for (int idx = 0; idx <= 3; idx++) {

        Temp_pol_vec[idx] = Polarization_Vector[idx] + RHS1[idx] * State_Vector[e_step] * 0.5;

    }

    for (int idx = 0; idx <= e_Dynamic_state_size - 1; idx++) {

        Temp_State_Vector[idx] = State_Vector[idx] + EOM1[idx] * State_Vector[e_step] * 0.5;

    }

    /* ================================== 2 ========================================== */

    Parallel_Transport_RHS(Temp_State_Vector, Temp_pol_vec, p_Spacetime, RHS2);
    p_Spacetime->get_EOM(Temp_State_Vector, EOM2);

    for (int idx = 0; idx <= 3; idx++) {

        Temp_pol_vec[idx] = Polarization_Vector[idx] + RHS2[idx] * State_Vector[e_step] * 0.5;

    }

    for (int idx = 0; idx <= e_Dynamic_state_size - 1; idx++) {

        Temp_State_Vector[idx] = State_Vector[idx] + EOM2[idx] * State_Vector[e_step] * 0.5;

    }

    /* ================================== 3 ========================================== */

    Parallel_Transport_RHS(Temp_State_Vector, Temp_pol_vec, p_Spacetime, RHS3);
    p_Spacetime->get_EOM(Temp_State_Vector, EOM3);

    for (int idx = 0; idx <= 3; idx++) {

        Temp_pol_vec[idx] = Polarization_Vector[idx] + RHS3[idx] * State_Vector[e_step];

    }

    for (int idx = 0; idx <= e_Dynamic_state_size - 1; idx++) {

        Temp_State_Vector[idx] = State_Vector[idx] + EOM3[idx] * State_Vector[e_step];

    }

    /* ================================== 4 ========================================== */

    Parallel_Transport_RHS(Temp_State_Vector, Temp_pol_vec, p_Spacetime, RHS4);

    for (int idx = 0; idx <= 3; idx++) {

        Polarization_Vector[idx] += (RHS1[idx] + 2.0 * RHS2[idx] + 2.0 * RHS3[idx] + RHS4[idx]) * State_Vector[e_step] * (1.0 / 6);

    }

    Metric_type s_Metric = p_Spacetime->get_metric(State_Vector);

    std::complex<double> Norm_Squared{};

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            Norm_Squared += s_Metric.Metric[left_idx][right_idx] * Polarization_Vector[left_idx] * std::conj(Polarization_Vector[right_idx]);

        }

    }

    double Norm = sqrt(std::abs(Norm_Squared));

    for (int idx = 0; idx <= 3; idx++) {

        Polarization_Vector[idx] /= Norm;

    }

}

void static Map_Polarization_Vector_to_Stokes(const double inv_Stokes_Tetrad[4][4],
                                              const std::complex<double>* Coord_Basis_Pol_vec,
                                              double* const Stokes_Vector) {

    std::complex<double> Stokes_Basis_Pol_vec[4]{};
    
    for (int stokes_idx = 0; stokes_idx <= 3; stokes_idx++) {
    
        for (int coord_idx = 0; coord_idx <= 3; coord_idx++) {
    
            Stokes_Basis_Pol_vec[stokes_idx] += inv_Stokes_Tetrad[stokes_idx][coord_idx] * Coord_Basis_Pol_vec[coord_idx];

        }
    
    }

    std::complex<double> test{};

    for (int stokes_idx = 0; stokes_idx <= 3; stokes_idx++) {

        test += Stokes_Basis_Pol_vec[stokes_idx] * std::conj(Stokes_Basis_Pol_vec[stokes_idx]);

    }

    double Norm = sqrt(std::abs(test));

    if (!isnan(1.0 / Norm) && !isinf(1.0 / Norm)) {

        for (int stokes_idx = 0; stokes_idx <= 3; stokes_idx++) {

            Stokes_Basis_Pol_vec[stokes_idx] /= Norm;

        }
    }

    double Polarized_Intensity = sqrt(Stokes_Vector[Q] * Stokes_Vector[Q] +
                                      Stokes_Vector[U] * Stokes_Vector[U] +
                                      Stokes_Vector[V] * Stokes_Vector[V]);

    Stokes_Vector[Q] = (Polarized_Intensity * (Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[1]) -
                                               Stokes_Basis_Pol_vec[2] * std::conj(Stokes_Basis_Pol_vec[2]))).real();
    
    Stokes_Vector[U] = (Polarized_Intensity * (Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[2]) +
                                               Stokes_Basis_Pol_vec[2] * std::conj(Stokes_Basis_Pol_vec[1]))).real();
    
    Stokes_Vector[V] = (Polarized_Intensity * (std::conj(Stokes_Basis_Pol_vec[1]) * Stokes_Basis_Pol_vec[2] -
                                               Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[2]))).imag();

    double Polarized_Intensity_after = sqrt(Stokes_Vector[Q] * Stokes_Vector[Q] +
                                            Stokes_Vector[U] * Stokes_Vector[U] +
                                            Stokes_Vector[V] * Stokes_Vector[V]);

    if (fabs(Polarized_Intensity - Polarized_Intensity_after) > 1e-6) {

        int test{};

    }



    if (!isinf(Polarized_Intensity / Polarized_Intensity_after) && !isnan(Polarized_Intensity / Polarized_Intensity_after)) {

        for (int idx = 1; idx <= 3; idx++) {

            Stokes_Vector[idx] *= Polarized_Intensity / Polarized_Intensity_after;

        }

    }

}

void static Map_Stokes_to_Polarization_Vector(const double* const Stokes_Vector,
                                              const double Stokes_Tetrad[4][4],
                                              std::complex<double>* const Coord_Basis_Pol_vec) {

    for (int index = 0; index <= 3; index++) {

        Coord_Basis_Pol_vec[index] = 0;

    }

    std::complex<double> Stokes_Basis_Pol_vec[4];

    double Polarized_Intensity = sqrt(Stokes_Vector[Q] * Stokes_Vector[Q] +
                                      Stokes_Vector[U] * Stokes_Vector[U] +
                                      Stokes_Vector[V] * Stokes_Vector[V]);

    Stokes_Basis_Pol_vec[1] = M_SQRT1_2;

    if (!isinf(Stokes_Vector[Q] / Polarized_Intensity) && !isnan(Stokes_Vector[Q] / Polarized_Intensity)) {

        Stokes_Basis_Pol_vec[1] = sqrt((1 + Stokes_Vector[Q] / Polarized_Intensity) / 2);

    }

    Stokes_Basis_Pol_vec[2] = 1.;

    if (!isinf(1. / std::norm(Stokes_Basis_Pol_vec[1] * Polarized_Intensity)) && !isnan(1. / std::norm(Stokes_Basis_Pol_vec[1] * Polarized_Intensity))) {

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
                                       int const N_theta_turning_points) {

    int Current_theta_turning_points = N_theta_turning_points;
    const int Max_order = compute_image_order(N_theta_turning_points, p_Sim_Context->p_Init_Conditions);
    int Current_order = Max_order;

    double Stokes_Vector[e_Stokes_param_num]{};
    double Stokes_Vector_offset[e_Stokes_param_num][e_order_number]{};

    std::complex<double> Coord_Basis_Pol_vec[4]{};
    std::complex<double> PW_const;
    // TODO: Propagate this aswell
    double Optical_Depth{};

    double LP_fraction{};

    for (int log_index = p_Ray_results->Ray_log_struct.Log_length - 1; log_index > 0; log_index--) {

        /* =============== Pick out the ray position / momenta from the Log, at the given log index =============== */

        double* Logged_ray_path = &(p_Ray_results->Ray_log_struct.Ray_path_log[log_index * e_Full_state_size]);

        Current_theta_turning_points -= Check_for_theta_turning_point(Logged_ray_path, Logged_ray_path - e_Full_state_size);

        if (Current_order != compute_image_order(Current_theta_turning_points, p_Sim_Context->p_Init_Conditions)) {

            for (int Stokes_idx = 0; Stokes_idx <= e_Stokes_param_num - 1; Stokes_idx++) {

                Stokes_Vector_offset[Stokes_idx][Current_order] = Stokes_Vector[Stokes_idx];

            }

            Current_order = compute_image_order(Current_theta_turning_points, p_Sim_Context->p_Init_Conditions);

        }

        /* ====================================== Propagate the radiative transfer equations ====================================== */

        bool Inside_emission_medium = false;

        if (Current_order <= p_Sim_Context->p_Init_Conditions->Max_order) {

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

            if (vector_norm(total_Transfer_functions.Emission_functions, 4) > 0 || 
                vector_norm(total_Transfer_functions.Absorbtion_functions, 4) > 0 || 
                vector_norm(total_Transfer_functions.Faradey_functions, 4) > 0) {

                Inside_emission_medium = true;

                Propagate_Stokes_vector(Analytic, total_Transfer_functions, Logged_ray_path[e_step], Stokes_Vector);

            }

        }

        /* ======================================================================================================================== */

        /* ====================================== Parallel transport the polarization vector ====================================== */

        double Tetrad[4][4]{};
        double inv_Tetrad[4][4]{};

        LP_fraction = sqrt(Stokes_Vector[Q] * Stokes_Vector[Q] + Stokes_Vector[U] * Stokes_Vector[U]) / Stokes_Vector[I];

        if (LP_fraction < 0.6 && !isnan(LP_fraction)) {

            int stop{};

        }

        if (p_Sim_Context->p_Init_Conditions->Observer_params.include_polarization && Stokes_Vector[I] > 0) {

            if (Inside_emission_medium) {

                /* Literally arbitrary... (I picked something that does not align with the plasma velocity, because Ive noticed that such cases can cause problems).
                   NOTE: This vector CAN be an arbitrary spacelike vector. Only the final mapping at the observer needs to follow conventions. */
                double Trial_spacelike_vector[4] = { 0, 0, -1, 0 };

                /* We are still inside the emission medium and need to parallel transport along with evaluating the emission.
                   This nessecitates a mapping back and forth between the local Stokes basis. */
                if (OK != Construct_Stokes_Tetrad(Tetrad, inv_Tetrad, p_Sim_Context, Logged_ray_path, Trial_spacelike_vector, false)) {

                    /* Sometimes the magnetic field lines up with the photon wave vector in a way that makes creating the tetrad awkward. In such cases, attempt to 
                       create the tetrad with the Trial_spacelike_vector. */
                    if (OK != Construct_Stokes_Tetrad(Tetrad, inv_Tetrad, p_Sim_Context, Logged_ray_path, Trial_spacelike_vector, true)) {

                        std::cout << "Could not construct the Stokes basis at the next ray point! \n";

                        exit(ERROR);

                    }

                }

                Map_Stokes_to_Polarization_Vector(std::as_const(Stokes_Vector), std::as_const(Tetrad), Coord_Basis_Pol_vec);

                /* -------------------------------- The actual parallel transport step ------------------------------- */

                Parallel_Transport_Polarization_Vector(Logged_ray_path, p_Sim_Context->p_Spacetime, Coord_Basis_Pol_vec);

                /* --------------------------------------------------------------------------------------------------- */

                /* The inverse mapping needs to be done at the point where the polarization vector is defined. After the parallel 
                   transport step, this point is further along the ray (which in this loop is further backwards along the ray log). */
                if (OK != Construct_Stokes_Tetrad(Tetrad, inv_Tetrad, p_Sim_Context, Logged_ray_path - e_Full_state_size, Trial_spacelike_vector, false)) {

                    /* Sometimes the magnetic field lines up with the photon wave vector in a way that makes creating the tetrad awkward. In such cases, attempt to 
                       create the tetrad with the Trial_spacelike_vector. */
                    if (OK != Construct_Stokes_Tetrad(Tetrad, inv_Tetrad, p_Sim_Context, Logged_ray_path - e_Full_state_size, Trial_spacelike_vector, true)) {

                        std::cout << "Could not construct the Stokes basis at the next ray point! \n";

                        exit(ERROR);

                    }

                }

                Map_Polarization_Vector_to_Stokes(std::as_const(inv_Tetrad), Coord_Basis_Pol_vec, Stokes_Vector);

                LP_fraction = sqrt(Stokes_Vector[Q] * Stokes_Vector[Q] + Stokes_Vector[U] * Stokes_Vector[U]) / Stokes_Vector[I];

            }
            else {

                /* In vacuum the Stokes vector is constant, so we don't need to perform mappings between it and the polarization vector.
                   We thus just run the parallel transport routine. */

                /* -------------------------------- The actual parallel transport step -------------------------------- */

                Parallel_Transport_Polarization_Vector(Logged_ray_path, p_Sim_Context->p_Spacetime, Coord_Basis_Pol_vec);

                /* ---------------------------------------------------------------------------------------------------- */

            }

        }

        /* ======================================================================================================================== */

        log_ray_emission(Stokes_Vector, Optical_Depth, p_Ray_results, log_index, p_Sim_Context->p_Init_Conditions->Simulation_mode);

    }

    if (p_Sim_Context->p_Init_Conditions->Observer_params.include_polarization) {

        /* =============== The final mapping of the polarizatio vector to Stokes parameters at the observer ===================== */

        double Observer_Tetrad[4][4]{};
        double Observer_inv_Tetrad[4][4]{};

        double Local_north_vector[4] = { 0, 0, -1, 0 };

        if (OK != Construct_Stokes_Tetrad(Observer_Tetrad, Observer_inv_Tetrad, p_Sim_Context, p_Ray_results->Ray_log_struct.Ray_path_log, Local_north_vector, true)) {

            /* There is no second attempt to create the observer tetrad, because it should never fail. */

            std::cout << "Could not construct the Stokes basis at the observer! \n";

            exit(ERROR);

        }

        Map_Polarization_Vector_to_Stokes(std::as_const(Observer_inv_Tetrad), Coord_Basis_Pol_vec, Stokes_Vector);

    }
    /* ====================================================================================================================== */

    for (int Stokes_idx = 0; Stokes_idx <= e_Stokes_param_num - 1; Stokes_idx++) {

        Stokes_Vector_offset[Stokes_idx][e_direct] = Stokes_Vector[Stokes_idx];

    }

    Seperate_Image_into_orders(Max_order, Stokes_Vector_offset, p_Ray_results);

}

void Propagate_ray(const Simulation_Context_type* const p_Sim_Context, Results_type* const p_Ray_results) {

    // Initialize the State Vectors
    double State_Vector[e_Full_state_size]{};
    double Old_State_Vector[e_Full_state_size]{};

    State_Vector[e_t]       = p_Sim_Context->p_Init_Conditions->Observer_params.init_time;
    State_Vector[e_r]       = p_Sim_Context->p_Init_Conditions->Observer_params.distance;
    State_Vector[e_theta]   = p_Sim_Context->p_Init_Conditions->Observer_params.inclination;
    State_Vector[e_phi]     = p_Sim_Context->p_Init_Conditions->Observer_params.azimuth;
    State_Vector[e_p_phi]   = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_phi];
    State_Vector[e_p_theta] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_theta];
    State_Vector[e_p_r]     = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_r];
    State_Vector[e_p_t]     = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_t];
    State_Vector[e_step]    = p_Sim_Context->p_Init_Conditions->Integrator_params.Init_stepzie; 
    State_Vector[e_affine_param] = 0;

    // Set the Old State Vector to the Initial State Vector
    memcpy(Old_State_Vector, State_Vector, e_Full_state_size * sizeof(double));

    for (int Image_order = e_direct; Image_order < e_order_number; Image_order++) {

        p_Ray_results->Photon_Momentum[e_phi][Image_order] = State_Vector[e_p_phi];
        p_Ray_results->Photon_Momentum[e_t][Image_order]   = State_Vector[e_p_t];

    }

    p_Ray_results->Metric_parameters = p_Sim_Context->p_Init_Conditions->Metric_parameters;

    // Initialize counters for the Number Of Integration Steps and the Number Of Turning points of the Polar Coordinate
    int integration_count{}, N_theta_turning_points{};

    // Calculate the image coordinates from the initial conditions
    get_image_coordinates(p_Sim_Context->p_Init_Conditions, p_Ray_results->Image_Coords);

    Step_controller controller(p_Sim_Context->p_Init_Conditions->Integrator_params);

    log_ray_path(State_Vector, p_Ray_results, p_Sim_Context->p_Init_Conditions);

    while (!controller.integration_complete && integration_count <= controller.Parameters.Max_integration_count && std::abs(State_Vector[e_affine_param]) <= controller.Parameters.Max_affine_param) {

        RK45(State_Vector, &controller, p_Sim_Context);

        if (controller.continue_integration) {

            integration_count += 1;
            p_Ray_results->Ray_log_struct.Log_offset = integration_count;

            log_ray_path(State_Vector, p_Ray_results, p_Sim_Context->p_Init_Conditions);

            N_theta_turning_points += Check_for_theta_turning_point(State_Vector, Old_State_Vector);

            /* ======================================== Evaluate the thin disk models ======================================== */

            if (e_Page_Thorne == p_Sim_Context->p_Init_Conditions->Disk_params.e_Disk_model) {

                Evaluate_Equatorial_Disk(p_Sim_Context, p_Ray_results, State_Vector, Old_State_Vector, N_theta_turning_points);

            }

            /* ============================================================================================================== */

            memcpy(Old_State_Vector, State_Vector, e_Full_state_size * sizeof(double));

        }

    }

    if (integration_count >= controller.Parameters.Max_integration_count) { std::cout << "Max iterations reached! \n"; }

    if (std::abs(State_Vector[e_affine_param]) >= controller.Parameters.Max_affine_param) { std::cout << "Max affine parameter value reached! \n"; };

    p_Ray_results->Ray_log_struct.Log_length = p_Ray_results->Ray_log_struct.Log_offset + 1;

    /* =========== Integrate the radiative transfer equations forward along the ray for the RIAF models =========== */

    if (e_Page_Thorne != p_Sim_Context->p_Init_Conditions->Disk_params.e_Disk_model) {

        Propagate_forward_emission(p_Sim_Context, p_Ray_results, N_theta_turning_points);

    }

    /* ============================================================================================================ */

}