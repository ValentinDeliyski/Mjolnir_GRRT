#define _USE_MATH_DEFINES
#include "IO_files.h"
#include "Enumerations.h"
#include "Constants.h"
#include "Spacetimes.h"

#include "Novikov_Thorne_model.h"
#include "Emission_Models.h"
#include "Disk_Models.h"

#include "General_GR_functions.h"
#include "General_math_functions.h"

#include "Rendering_Engine.h"
#include "Radiative_Transfer.h"

#include "Lensing.h"
#include "Integrators.h"

#include "Parallel_transport.h"

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
    
    for (int stokes_idx = 0; stokes_idx < e_Stokes_param_num; stokes_idx++) {

        p_Ray_Results->Ray_log_struct.Ray_emission_log[stokes_idx][0 + 2 * log_index] = Stokes_Vector[stokes_idx] ;
        p_Ray_Results->Ray_log_struct.Ray_emission_log[stokes_idx][1 + 2 * log_index] = Optical_depth;

    }

}

static void Propagate_Stokes_vector(Radiative_Transfer_Integrator e_Integrator,
                                    const Simulation_Context_type* p_Sim_Context,
                                    double* const State_Vector,
                                    double* const Stokes_Vector) {

    Transfer_functions_type Total_Transfer_Functions{};

    /* Loop trough each emission medium (Disk, Hotspot, Jet and so on) and sum their respective transfer functions */
    for (int emission_medium = Disk; emission_medium <= Hotspot; emission_medium++) {

        Transfer_functions_type Temp_Transfer_functions{};

        p_Sim_Context->p_Emission_Model->get_radiative_transfer_functions(State_Vector,
                                                                          p_Sim_Context,
                                                                          static_cast<Emission_medium_enums>(emission_medium),
                                                                          &Temp_Transfer_functions);

        add_vectors(Temp_Transfer_functions.Emission_functions, Total_Transfer_Functions.Emission_functions, e_Stokes_param_num, Total_Transfer_Functions.Emission_functions);
        add_vectors(Temp_Transfer_functions.Faradey_functions, Total_Transfer_Functions.Faradey_functions, e_Stokes_param_num, Total_Transfer_Functions.Faradey_functions);
        add_vectors(Temp_Transfer_functions.Absorbtion_functions, Total_Transfer_Functions.Absorbtion_functions, e_Stokes_param_num, Total_Transfer_Functions.Absorbtion_functions);

    }

    if ((vector_norm(Total_Transfer_Functions.Emission_functions, 4) +
        vector_norm(Total_Transfer_Functions.Faradey_functions, 4) +
        vector_norm(Total_Transfer_Functions.Absorbtion_functions, 4)) < std::numeric_limits<double>::min()) {

        return;

    }

    switch (e_Integrator) {

    case Analytic:

        Analytic_Radiative_Transfer(const_cast<double*>(Total_Transfer_Functions.Emission_functions),
                                    const_cast<double*>(Total_Transfer_Functions.Absorbtion_functions),
                                    const_cast<double*>(Total_Transfer_Functions.Faradey_functions),
                                    State_Vector[e_step] * MASS_TO_CM * p_Sim_Context->p_Init_Conditions->central_object_mass,
                                    Stokes_Vector);

        break;

    case Implicit_Trapezoid:

        Implicit_Trapezoid_Radiative_Transfer(const_cast<double*>(Total_Transfer_Functions.Emission_functions),
                                              const_cast<double*>(Total_Transfer_Functions.Absorbtion_functions),
                                              const_cast<double*>(Total_Transfer_Functions.Faradey_functions),
                                              State_Vector[e_step] * MASS_TO_CM * p_Sim_Context->p_Init_Conditions->central_object_mass,
                                              Stokes_Vector);
        break;

    default:

        RK5_radiative_transfer(const_cast<double*>(Total_Transfer_Functions.Emission_functions),
                               const_cast<double*>(Total_Transfer_Functions.Absorbtion_functions),
                               const_cast<double*>(Total_Transfer_Functions.Faradey_functions),
                               State_Vector,
                               p_Sim_Context,
                               Stokes_Vector);

        break;

    }

}

Return_Values static Construct_Stokes_Tetrad(double Tetrad[4][4],
                                             double inv_Tetrad[4][4],
                                             const Simulation_Context_type* const p_Sim_Context,
                                             const bool Overwride,
                                             const double* const State_vector) {

    memset(Tetrad, 0, 16 * sizeof(double));
    memset(inv_Tetrad, 0, 16 * sizeof(double));

    /* The reference of this implementation is the second RAPTOR paper: https://arxiv.org/pdf/2007.03045.pdf 
       Expressions 10 - 11. Note that their definition of 9d is wrong... the "g" in the denominator should 
       be omega, as defined in 9c. */

    Metric_type s_Metric = p_Sim_Context->p_Spacetime->get_metric(State_vector);

    Return_Values Plasma_velocity_OK{};

    // ---------------- The velocity vector -> this is chosen to be the local emitter's velocity when we are inside the emission medium, 
    // otherwise it is chosen to be the observer's velocity. I use the theta dependant profile for the observer velocity, because it tends to be well defined below ISCO.
    double Plasma_velocity_contravariant[4]{};
    double Plasma_velocity_covariant[4]{};

    // ---------------- The tetrad requires a spacelike vector -> this is chosen to be the local magnetic field 4-vector when we are inside the emission medium
    double Spacelike_vector_contravariant[4]{};
    double Spacelike_vector_covariant[4]{};

    const double Wave_Vector_covariant[4] = { State_vector[e_p_t], State_vector[e_p_r], State_vector[e_p_theta], State_vector[e_p_phi] };

    // ---------------- We will need the state of all emission media
    Emission_medium_state_type s_Disk_state{};
    Emission_medium_state_type s_Hotspot_state{};

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

        if (OK != Plasma_velocity_OK) { return ERROR; }

        p_Sim_Context->p_Emission_Model->get_magnetic_field(State_vector, &s_Metric, &s_Hotspot_state);
        memcpy(Spacelike_vector_contravariant, s_Hotspot_state.Magnetic_fields.B_field_plasma_frame, 4 * sizeof(double));
        memcpy(Plasma_velocity_contravariant, s_Hotspot_state.Plasma_Velocity, 4 * sizeof(double));

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
        p_Sim_Context->p_Emission_Model->get_magnetic_field(State_vector, &s_Metric, &s_Disk_state);

        Plasma_velocity_OK = p_Sim_Context->p_Emission_Model->get_plasma_velocity(State_vector,
                                                                                  p_Sim_Context,
                                                                                  p_Sim_Context->p_Init_Conditions->Disk_params.Velocity_profile_type,
                                                                                  p_Sim_Context->p_Init_Conditions->Disk_params.Radial_velocity_fraction,
                                                                                  s_Disk_state.Plasma_Velocity);

        if (OK != Plasma_velocity_OK) { return ERROR; }

        memcpy(Spacelike_vector_contravariant, s_Disk_state.Magnetic_fields.B_field_plasma_frame, 4 * sizeof(double));
        memcpy(Plasma_velocity_contravariant, s_Disk_state.Plasma_Velocity, 4 * sizeof(double));

    }
    else if (Overwride) {

        double Local_north_vector[4] = { 0, 0, -1, 0 };

        memcpy(Spacelike_vector_contravariant, Local_north_vector, 4 * sizeof(double));
        memcpy(Plasma_velocity_contravariant, p_Sim_Context->p_Observer->get_obs_velocity(), 4 * sizeof(double));

    }
    else { return NOT_IN_EMISSION_MEDIUM; }

    /* --------------------- Evaluate the inner products (expressions 9a - 9d), and compute the contravariant Wave-Vector --------------------- */

    double Wave_vec_dot_Plasma_vel{}, Plasma_vel_dot_Spacelike_vec{}, Spacelike_vec_norm_squared{}, Wave_vec_dot_Spacelike_vec{};

    for (int left_idx = 0; left_idx < 4; left_idx++) {

        Wave_vec_dot_Plasma_vel    += Wave_Vector_covariant[left_idx] * Plasma_velocity_contravariant[left_idx];
        Wave_vec_dot_Spacelike_vec += Wave_Vector_covariant[left_idx] * Spacelike_vector_contravariant[left_idx];

        for (int right_idx = 0; right_idx < 4; right_idx++) {

            Plasma_vel_dot_Spacelike_vec += s_Metric.Metric[left_idx][right_idx] * Plasma_velocity_contravariant[left_idx] * Spacelike_vector_contravariant[right_idx];
            Spacelike_vec_norm_squared   += s_Metric.Metric[left_idx][right_idx] * Spacelike_vector_contravariant[left_idx] * Spacelike_vector_contravariant[right_idx];

        }
        
    }

    const double metric_determinant = get_metric_det(s_Metric.Metric);
    const double sqrt_determinant = sqrt(-metric_determinant);

    const double C_coeff = -Wave_vec_dot_Spacelike_vec / Wave_vec_dot_Plasma_vel - Plasma_vel_dot_Spacelike_vec;
    const double N_coeff = sqrt((Spacelike_vec_norm_squared + Plasma_vel_dot_Spacelike_vec * Plasma_vel_dot_Spacelike_vec - C_coeff * C_coeff));

    /* --- Perform checks on these coefficients, because for very low plasma densities they blow up --- */

    if (isnan(N_coeff) || isinf(1.0 / N_coeff) || isnan(C_coeff)) { 
        
        return ERROR; }

    /* --- Raise / Lower indicies on the three main 4-vectors - this is needed for computing the final tetrad vector --- */

    double inv_Metric[4][4]{};
    invert_metric(inv_Metric, s_Metric.Metric);

    double Wave_vector_contravariant[4]{};

    for (int left_idx = 0; left_idx < 4; left_idx++) {

        for (int right_idx = 0; right_idx < 4; right_idx++) {

            Wave_vector_contravariant[left_idx] += inv_Metric[left_idx][right_idx] * Wave_Vector_covariant[right_idx];

            Spacelike_vector_covariant[left_idx] += s_Metric.Metric[left_idx][right_idx] * Spacelike_vector_contravariant[right_idx];
            Plasma_velocity_covariant[left_idx]  += s_Metric.Metric[left_idx][right_idx] * Plasma_velocity_contravariant[right_idx];

        }

    }

    /* --------------------- Compute the first two tetrad basis vectors --------------------- */

    // Wave_vec_dot_Plasma_vel is safe to divide by as it can never go to zero - Plasma_vel is never null.

    for (int index = 0; index < 4; index++) {

        Tetrad[e_t][index]   = Plasma_velocity_contravariant[index];
        Tetrad[e_phi][index] = -Wave_vector_contravariant[index] / Wave_vec_dot_Plasma_vel - Plasma_velocity_contravariant[index];

    }

    /* -------- Compute the contravariant 4D permutation symbol (this I lifted straight from the RAPTOR code...) -------- */

    double Levi_Cevita_tensor[4][4][4][4]{};

    for (int i = 0; i < 4; i++) {

        for (int j = 0; j < 4; j++) {

            for (int k = 0; k < 4; k++) {

                for (int l = 0; l < 4; l++) {

                    Levi_Cevita_tensor[i][j][k][l] = -((i - j) * (i - k) * (i - l) * (j - k) * (j - l) * (k - l) / 12.) / sqrt_determinant;
                                                   
                }
            }
        }
    }

    /* --------------------- Compute the last two tetrad basis vectors --------------------- */
 
    for (int index = 0; index < 4; index++) {

        Tetrad[e_r][index] = (Spacelike_vector_contravariant[index] + Plasma_vel_dot_Spacelike_vec * Plasma_velocity_contravariant[index] - C_coeff * Tetrad[e_phi][index]) / N_coeff;

        for (int i = 0; i < 4; i++) {

            for (int j = 0; j < 4; j++) {

                for (int k = 0; k < 4; k++) {

                    Tetrad[e_theta][index] += Levi_Cevita_tensor[index][i][j][k] *
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

    for (int left_idx = 0; left_idx < 4; left_idx++) {

        for (int right_idx = 0; right_idx < 4; right_idx++) {

            for (int m = 0; m < 4; m++) {

                for (int g = 0; g < 4; g++) {

                    inv_Tetrad[left_idx][right_idx] += Minkowski_Metric[left_idx][m] * s_Metric.Metric[right_idx][g] * Tetrad[m][g];

                        
                }
            }

            if (isnan(inv_Tetrad[left_idx][right_idx]) || isinf(inv_Tetrad[left_idx][right_idx]) ||
                isnan(Tetrad[left_idx][right_idx]) || isinf(Tetrad[left_idx][right_idx])) {

                return ERROR;

            }

        }
    }

    return OK;

}

void static Map_Polarization_Vector_to_Stokes(const double inv_Stokes_Tetrad[4][4],
                                              const std::complex<double>* Coord_Basis_Pol_vec,
                                              double* const Stokes_Vector) {

    std::complex<double> Stokes_Basis_Pol_vec[4]{};
    
    for (int stokes_idx = 0; stokes_idx < 4; stokes_idx++) {
    
        for (int coord_idx = 0; coord_idx < 4 ; coord_idx++) {
    
            Stokes_Basis_Pol_vec[stokes_idx] += inv_Stokes_Tetrad[stokes_idx][coord_idx] * Coord_Basis_Pol_vec[coord_idx];

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

    double Polarized_Intensity_after_mapping = sqrt(Stokes_Vector[Q] * Stokes_Vector[Q] +
        Stokes_Vector[U] * Stokes_Vector[U] +
        Stokes_Vector[V] * Stokes_Vector[V]);

    if (!isnan(1. / Polarized_Intensity_after_mapping) && !isinf(1.0 / Polarized_Intensity_after_mapping)) {

        for (int idx = 1; idx < 4; idx++) {

            Stokes_Vector[idx] *= Polarized_Intensity / Polarized_Intensity_after_mapping;

        }

    }

}

void static Map_Stokes_to_Polarization_Vector(const double* const Stokes_Vector,
                                              const double Stokes_Tetrad[4][4],
                                              std::complex<double>* const Coord_Basis_Pol_vec) {

    memset(Coord_Basis_Pol_vec, 0, 4 * sizeof(std::complex<double>));

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

    for (int coord_idx = 0; coord_idx < 4; coord_idx++) {

        for (int stokes_idx = 0; stokes_idx < 4; stokes_idx++) {

            Coord_Basis_Pol_vec[coord_idx] += Stokes_Tetrad[stokes_idx][coord_idx] * Stokes_Basis_Pol_vec[stokes_idx];

        }

    }

}

bool static Is_inside_emission_medium(const Simulation_Context_type* const p_Sim_Context,
                                      const double* const State_Vector) {

    Emission_medium_state_type s_Hotspot_state{};
    Emission_medium_state_type s_Disk_state{};

    Return_Values Hotspot_velocity_OK = p_Sim_Context->p_Emission_Model->get_plasma_velocity(State_Vector,
                                                                                             p_Sim_Context,
                                                                                             p_Sim_Context->p_Init_Conditions->Hotspot_params.Velocity_profile_type,
                                                                                             p_Sim_Context->p_Init_Conditions->Hotspot_params.Radial_velocity_fraction,
                                                                                             s_Hotspot_state.Plasma_Velocity);

    bool In_hotspot = false;
    
    if (OK == Hotspot_velocity_OK) {

        In_hotspot = p_Sim_Context->p_Emission_Model->p_Hotspot_Model->is_inside_hotspot(State_Vector,
                                                                                         &s_Hotspot_state);
    }

    const bool In_disk = p_Sim_Context->p_Emission_Model->p_Disk_Model->is_inside_disk(State_Vector,
                                                                                       p_Sim_Context->p_Emission_Model->p_Disk_Model->s_Disk_params.e_Disk_model,
                                                                                       &s_Disk_state);

    return In_hotspot || In_disk;

}

bool static Evaluate_Equatorial_Disk(const Simulation_Context_type* const p_Sim_Context,
                                     Results_type* const p_Ray_results,
                                     const double* const State_vector,
                                     const double* const Old_state,
                                     const int N_theta_turning_points) {

    /* ------------ The number of components is e_State_Number - 1 because we do not include the integration step. */
    double Crossing_State[e_Dynamic_state_size]{};
    double& R_throat = p_Sim_Context->p_Init_Conditions->Metric_parameters.R_throat;

    if (interpolate_equatorial_crossing(State_vector, Old_state, Crossing_State)) {

        if (Wormhole == p_Sim_Context->p_Init_Conditions->Metric_parameters.e_Spacetime) {

            // The wormhole metric uses the global coordinate ell^2 = r^2 - r_throat^2
            // Here I convert back to the r coordinate for the NT model evaluation

            Crossing_State[e_r] = sqrt(Crossing_State[e_r] * Crossing_State[e_r] - R_throat * R_throat);

        }

        double& r_in = p_Sim_Context->p_Init_Conditions->Disk_params.Novikov_Thorne_params.r_in;
        double& r_out = p_Sim_Context->p_Init_Conditions->Disk_params.Novikov_Thorne_params.r_out;

        if (Crossing_State[e_r] < r_out && Crossing_State[e_r] > r_in && !p_Ray_results->NT_Disk_found) {

            p_Ray_results->Redshift_NT = get_redshift(Crossing_State, p_Sim_Context->p_NT_model->get_disk_velocity_vector(Crossing_State), p_Sim_Context->p_Observer);
            p_Ray_results->Flux_NT = p_Sim_Context->p_NT_model->get_flux(Crossing_State);

            double* Polarization_vector_coord = p_Sim_Context->p_NT_model->Construct_coord_polarization_vector(Crossing_State);

            /* ====================== Parallel transport the polarization vector back to the observer ====================== */

            for (int log_idx = p_Ray_results->Ray_log_struct.Log_offset; log_idx > 0; log_idx--) {

                double* Logged_State = &p_Ray_results->Ray_log_struct.Ray_path_log[log_idx * e_Full_state_size];

                Parallel_Transport_Vector(Logged_State, p_Sim_Context->p_Spacetime, Contravariant, Polarization_vector_coord);

            }

            double Polarization_vector_ZAMO[4]{};
            Contravariant_coord_to_ZAMO(&p_Sim_Context->p_Init_Conditions->Init_metric, Polarization_vector_coord, Polarization_vector_ZAMO);

            p_Ray_results->Projected_polarization_vector[e_x] = Polarization_vector_ZAMO[e_phi];
            p_Ray_results->Projected_polarization_vector[e_y] = Polarization_vector_ZAMO[e_theta];

            /* ============================================================================================================= */

            p_Ray_results->NT_Disk_found = true;

        }

        memcpy(p_Ray_results->Thin_Disk_State_Vector, Crossing_State, e_Dynamic_state_size * sizeof(double));

        return true;
    }

    return false;
}

void static Propagate_forward_emission(const Simulation_Context_type* const p_Sim_Context, 
                                       Results_type* const p_Ray_results,
                                       int const N_theta_turning_points) {

    int Current_theta_turning_points = N_theta_turning_points;
    int Current_order = compute_image_order(N_theta_turning_points, p_Sim_Context->p_Init_Conditions);

    double Stokes_Vector[e_Stokes_param_num]{};

    std::complex<double> Coord_Basis_Pol_vec[4]{};
    // TODO: Propagate this aswell
    double Optical_Depth{};   

    for (int log_index = p_Ray_results->Ray_log_struct.Log_length; log_index >= 0; log_index--) {

        /* =============== Pick out the ray position / momenta from the Log, at the given log index =============== */

        double* Logged_ray_path = &(p_Ray_results->Ray_log_struct.Ray_path_log[log_index * e_Full_state_size]);
        
        if (log_index > 0) { Current_theta_turning_points -= Check_for_theta_turning_point(Logged_ray_path, Logged_ray_path - e_Full_state_size); };

        Current_order = compute_image_order(Current_theta_turning_points, p_Sim_Context->p_Init_Conditions);

        log_ray_emission(Stokes_Vector, Optical_Depth, p_Ray_results, log_index, p_Sim_Context->p_Init_Conditions->Simulation_mode);

        if (Is_inside_emission_medium(p_Sim_Context, Logged_ray_path)) {

            double Tetrad[4][4]{};
            double inv_Tetrad[4][4]{};

            if (p_Sim_Context->p_Init_Conditions->Observer_params.include_polarization) {

                switch (Construct_Stokes_Tetrad(Tetrad, inv_Tetrad, p_Sim_Context, false, Logged_ray_path)) {

                case ERROR:

                    throw std::runtime_error("Could not construct the Stokes basis at the current ray point!\n");
                    break;

                case NOT_IN_EMISSION_MEDIUM:
                    continue;

                default:

                    Map_Polarization_Vector_to_Stokes(std::as_const(inv_Tetrad), std::as_const(Coord_Basis_Pol_vec), Stokes_Vector);
                    break;

                }
            }

            /* ====================================== Propagate the radiative transfer equations ====================================== */

            if (Current_order >= p_Sim_Context->p_Init_Conditions->Min_order && Current_order <= p_Sim_Context->p_Init_Conditions->Max_order) {

                Propagate_Stokes_vector(p_Sim_Context->p_Init_Conditions->Integrator_params.e_Radiative_transfer_integrator, p_Sim_Context, Logged_ray_path, Stokes_Vector);

            }

            if (p_Sim_Context->p_Init_Conditions->Observer_params.include_polarization) {

                Map_Stokes_to_Polarization_Vector(std::as_const(Stokes_Vector), std::as_const(Tetrad), Coord_Basis_Pol_vec);

            }
        }

        /* ======================================================================================================================== */

        /* ====================================== Parallel transport the polarization vector ====================================== */

        if (p_Sim_Context->p_Init_Conditions->Observer_params.include_polarization) {

                Parallel_Transport_Vector(Logged_ray_path, p_Sim_Context->p_Spacetime, Contravariant, Coord_Basis_Pol_vec);

        }

        /* ======================================================================================================================== */

    }

    /* =============== The final mapping of the polarization vector to Stokes parameters at the observer ===================== */

    if (p_Sim_Context->p_Init_Conditions->Observer_params.include_polarization && Stokes_Vector[I] > 0) {

        double Observer_Tetrad[4][4]{};
        double Observer_inv_Tetrad[4][4]{};

        if (OK != Construct_Stokes_Tetrad(Observer_Tetrad, Observer_inv_Tetrad, p_Sim_Context, true, p_Ray_results->Ray_log_struct.Ray_path_log)) {

            throw std::runtime_error("Could not construct the Stokes basis at the observer! \n");

        }

        Map_Polarization_Vector_to_Stokes(std::as_const(Observer_inv_Tetrad), Coord_Basis_Pol_vec, Stokes_Vector);

        /*std::complex<double> ZAMO_pol_vec[4]{};

        Contravariant_coord_to_ZAMO(&p_Sim_Context->p_Init_Conditions->Init_metric, Coord_Basis_Pol_vec, ZAMO_pol_vec);

        double Polarized_Intensity = sqrt(Stokes_Vector[Q] * Stokes_Vector[Q] +
            Stokes_Vector[U] * Stokes_Vector[U] +
            Stokes_Vector[V] * Stokes_Vector[V]);

        Stokes_Vector[Q] = -(Polarized_Intensity * (ZAMO_pol_vec[e_phi] * std::conj(ZAMO_pol_vec[e_phi]) -
            ZAMO_pol_vec[e_theta] * std::conj(ZAMO_pol_vec[e_theta]))).real();

        Stokes_Vector[U] = (Polarized_Intensity * (ZAMO_pol_vec[e_phi] * std::conj(ZAMO_pol_vec[e_theta]) +
            ZAMO_pol_vec[e_theta] * std::conj(ZAMO_pol_vec[e_phi]))).real();

        Stokes_Vector[V] = (Polarized_Intensity * (std::conj(ZAMO_pol_vec[e_phi]) * ZAMO_pol_vec[e_theta] -
            ZAMO_pol_vec[e_phi] * std::conj(ZAMO_pol_vec[e_theta]))).imag();*/

    }

    /* ====================================================================================================================== */

    memcpy(p_Ray_results->Intensity, Stokes_Vector, 4 * sizeof(double));

}

void Propagate_ray(const Simulation_Context_type* const p_Sim_Context, Results_type* const p_Ray_results) {

    int N_theta_turning_points{}, Current_order{};

    // Calculate the image coordinates from the initial conditions
    get_image_coordinates(p_Sim_Context->p_Init_Conditions, p_Ray_results->Image_Coords);

    Integrator_class Geodesic_Integrator(p_Sim_Context, p_Ray_results);

    p_Ray_results->NT_Disk_found = false;

    while (!Geodesic_Integrator.integration_complete) {

        Geodesic_Integrator.Propagate_ray();

        if (Geodesic_Integrator.continue_integration) {

            N_theta_turning_points += Check_for_theta_turning_point(Geodesic_Integrator.get_current_State_Vector(), Geodesic_Integrator.get_previous_State_Vector());

            Current_order = compute_image_order(N_theta_turning_points, p_Sim_Context->p_Init_Conditions);

            /* ======================================== Evaluate the thin disk models ======================================== */

            if (e_Novikov_Thorne == p_Sim_Context->p_Init_Conditions->Disk_params.e_Disk_model && 
                Current_order >= p_Sim_Context->p_Init_Conditions->Min_order && 
                Current_order <= p_Sim_Context->p_Init_Conditions->Max_order) {

                Evaluate_Equatorial_Disk(p_Sim_Context, p_Ray_results, Geodesic_Integrator.get_current_State_Vector(), Geodesic_Integrator.get_previous_State_Vector(), N_theta_turning_points);

            }

            /* ============================================================================================================== */

        }

    }

    p_Ray_results->Ray_log_struct.Log_length = p_Ray_results->Ray_log_struct.Log_offset;
    p_Ray_results->Metric_parameters      = p_Sim_Context->p_Init_Conditions->Metric_parameters;

    interpolate_celestial_sphere_crossing(Geodesic_Integrator.get_current_State_Vector(),
                                          Geodesic_Integrator.get_previous_State_Vector(), 
                                          p_Sim_Context->p_Init_Conditions->Metric_parameters.Scattering_radius, 
                                          p_Ray_results->Celestial_sphere_crossing_coords);

    memcpy(p_Ray_results->Final_State_Vector, Geodesic_Integrator.get_current_State_Vector(), e_Full_state_size * sizeof(double));

    /* =========== Integrate the radiative transfer equations forward along the ray for the RIAF models =========== */

    if (e_Novikov_Thorne != p_Sim_Context->p_Init_Conditions->Disk_params.e_Disk_model) {

        Propagate_forward_emission(p_Sim_Context, p_Ray_results, N_theta_turning_points);

    }

    /* ============================================================================================================ */

}