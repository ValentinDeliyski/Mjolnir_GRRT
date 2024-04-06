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
#include "Radiative_Transfer.h"

#include "Lensing.h"
#include <iostream>
#include <complex>

void static log_ray_path(double State_Vector[], Results_type* s_Ray_Results, Step_controller Controller){

    int &log_offset = s_Ray_Results->Ray_log_struct.Log_offset;

    // The loop continues up untill e_path_log_number - 2 in order to exclude the log step from the loop,
    // because its not parat of the photon state vector - I take care of it after the loop "by hand".

    for (int index = e_r; index <= e_path_log_number - 2; index++) {

        s_Ray_Results->Ray_log_struct.Ray_path_log[index + log_offset * e_path_log_number] = State_Vector[index];

        // The wormhole metric works with a "global" radial coordinate, that goes negative on the other side of the throat.
        // The emission model can't work with this coordinate, so I log the normal spherical radial coordinate instead. 

        if (e_metric == Wormhole && index == e_r) {

            s_Ray_Results->Ray_log_struct.Ray_path_log[e_r + log_offset * e_path_log_number] = sqrt(State_Vector[e_r] * State_Vector[e_r] + WH_R_THROAT * WH_R_THROAT);

        }

    }

    s_Ray_Results->Ray_log_struct.Ray_path_log[e_path_log_step + log_offset * e_path_log_number] = Controller.previous_step;

}

void static log_ray_emission(double Stokes_Vector[STOKES_PARAM_NUM], double Optical_depth, Results_type* s_Ray_Results, int log_index) {
    
    for (int stokes_idx = 0; stokes_idx <= STOKES_PARAM_NUM - 1; stokes_idx++) {

        s_Ray_Results->Ray_log_struct.Ray_emission_log[0 + 2 * log_index][stokes_idx] = Stokes_Vector[stokes_idx];
        s_Ray_Results->Ray_log_struct.Ray_emission_log[1 + 2 * log_index][stokes_idx] = Optical_depth;
    }

}

void static Evaluate_Equatorial_Disk(Initial_conditions_type* const s_Initial_Conditions,
                                  Results_type* const s_Ray_results,
                                  double State_vector[],
                                  double Old_state[], 
                                  int N_theta_turning_points) {

    double crossing_coords[3]{}, crossing_momenta[3]{};

    if (interpolate_crossing(State_vector, Old_state, crossing_coords, crossing_momenta, s_Initial_Conditions->NT_model)) {


        int Image_Order = compute_image_order(N_theta_turning_points, s_Initial_Conditions);

        double r_crossing = std::sqrt(crossing_coords[x] * crossing_coords[x] + crossing_coords[y] * crossing_coords[y]);
        double state_crossing[2] = { r_crossing, M_PI_2 };

        if (Evaluate_NT_disk){

            s_Ray_results->Redshift_NT[Image_Order] = s_Initial_Conditions->NT_model->Redshift(State_vector[e_p_phi], 
                                                                                               state_crossing, 
                                                                                               r_obs, 
                                                                                               theta_obs, 
                                                                                               s_Initial_Conditions->Spacetimes);

            if (s_Ray_results->Redshift_NT[Image_Order] > std::numeric_limits<double>::min()) {

                s_Ray_results->Flux_NT[Image_Order] = s_Initial_Conditions->NT_model->get_flux(r_crossing, s_Initial_Conditions->Spacetimes);

            }

        }

        s_Ray_results->Source_Coords[e_r][Image_Order] = r_crossing;
        s_Ray_results->Source_Coords[e_phi][Image_Order] = State_vector[e_phi];

        s_Ray_results->Three_Momentum[e_r][Image_Order] = crossing_momenta[e_r];
        s_Ray_results->Three_Momentum[e_theta][Image_Order] = crossing_momenta[e_theta];
    }

}

void static Seperate_Image_into_orders(int const Max_theta_turning_points, 
                                int const N_theta_turning_points,
                                Initial_conditions_type* const s_Initial_Conditions, 
                                Results_type* const p_Ray_results,
                                double const Intensity[STOKES_PARAM_NUM]) {

        int Image_Order = compute_image_order(Max_theta_turning_points - N_theta_turning_points, s_Initial_Conditions);

        for (int stokes_idx = 0; stokes_idx <= STOKES_PARAM_NUM - 1; stokes_idx++) {

            p_Ray_results->Intensity[Image_Order][stokes_idx] = Intensity[stokes_idx];

            for (int order_scan = Image_Order + 1; order_scan <= compute_image_order(Max_theta_turning_points, s_Initial_Conditions); order_scan++) {

                p_Ray_results->Intensity[Image_Order][stokes_idx] -= p_Ray_results->Intensity[order_scan][stokes_idx];

            }

        }

}

void static Propagate_Stokes_vector(Radiative_Transfer_Integrator Integrator,
                             double const emission_functions[2][STOKES_PARAM_NUM], 
                             double const absorbtion_functions[2][STOKES_PARAM_NUM],
                             double const faradey_functions[2][STOKES_PARAM_NUM], 
                             double const step, 
                             double Intensity[STOKES_PARAM_NUM]) {

    switch (Integrator) {

    case Analytic:

        Analytic_Radiative_Transfer(emission_functions[Current], absorbtion_functions[Current], faradey_functions[Current], step / 2, Intensity);
        Analytic_Radiative_Transfer(emission_functions[Next], absorbtion_functions[Next], faradey_functions[Next], step / 2, Intensity);

        break;

    case Implicit_Trapezoid:

        Implicit_Trapezoid_Radiative_Transfer(emission_functions[Current], absorbtion_functions[Current], faradey_functions[Current], step / 2, Intensity);
        Implicit_Trapezoid_Radiative_Transfer(emission_functions[Next], absorbtion_functions[Next], faradey_functions[Next], step / 2, Intensity);

        break;

    case RK4:

        RK4_Radiative_Transfer(emission_functions, absorbtion_functions, faradey_functions, step, Intensity);

        break;

    default:

        std::cout << "Integration method not supported for the radiative transfer equations!" << '\n';

        break;

    }

}

Return_Values static Construct_Stokes_Tetrad(double Tetrad[4][4],
                                             double inv_Tetrad[4][4],
                                             Initial_conditions_type* s_Initial_Conditions, 
                                             double State_vector[e_State_Number]) {

    /* The reference of this implementation is the second RAPTOR paper: https://arxiv.org/pdf/2007.03045.pdf 
       Expressions 10 - 11. Note that their definition of 9d is wrong... the "g" in the denominator should 
       be omega, as defined in 9c. */

     /* --------------------- Get the three 4-vectors from which we will construct the tetrad --------------------- */

     double* Plasma_velocity_contravariant = s_Initial_Conditions->OTT_model->get_disk_velocity(State_vector, s_Initial_Conditions);

     double B_field_contravariant[4]{};
     s_Initial_Conditions->OTT_model->get_magnetic_field(B_field_contravariant, State_vector, s_Initial_Conditions);

     double Wave_Vector_covariant[4] = { -1, State_vector[e_p_r], State_vector[e_p_theta], State_vector[e_p_phi] };

     /* --------------------- Evaluate the inner products (expressions 9a - 9d), and compute the contravariant Wave-Vector --------------------- */

     double Wave_vec_dot_Plasma_vel{}, Plasma_vel_dot_B_field{}, B_field_norm_squared{}, Wave_vec_dot_B_field{};

     Metric_type s_Metric = s_Initial_Conditions->Spacetimes[e_metric]->get_metric(State_vector);

     for (int left_idx = 0; left_idx <= 3; left_idx++) {

         Wave_vec_dot_Plasma_vel += Wave_Vector_covariant[left_idx] * Plasma_velocity_contravariant[left_idx];
         Wave_vec_dot_B_field    += Wave_Vector_covariant[left_idx] * B_field_contravariant[left_idx];

         for (int right_idx = 0; right_idx <= 3; right_idx++) {

             Plasma_vel_dot_B_field += s_Metric.Metric[left_idx][right_idx] * Plasma_velocity_contravariant[left_idx] * B_field_contravariant[right_idx];
             B_field_norm_squared   += s_Metric.Metric[left_idx][right_idx] * B_field_contravariant[left_idx] * B_field_contravariant[right_idx];

         }
         
     }

     double metric_determinant = get_metric_det(s_Metric.Metric);
     double sqrt_determinant = sqrt(-metric_determinant);

     double C_coeff = -Wave_vec_dot_B_field / Wave_vec_dot_Plasma_vel - Plasma_vel_dot_B_field;
     double N_coeff = sqrt(B_field_norm_squared + Plasma_vel_dot_B_field * Plasma_vel_dot_B_field - C_coeff * C_coeff);

     /* --- Perform checks on these coefficients, because for very low plasma densities they blow up --- */

     if (isnan(N_coeff) || isinf(1.0 / N_coeff) || isnan(C_coeff) || isinf(1.0 / C_coeff)) {

         return ERROR;

     }

     /* --- Raise / Lower indicies on the three main 4-vectors - this is needed for computing the final tetrad vector --- */

     double inv_Metric[4][4]{};
     invert_metric(inv_Metric, s_Metric.Metric);

     double Wave_vector_contravariant[4]{};
     double B_field_covariant[4]{};
     double Plasma_velocity_covavriant[4]{};

     for (int left_idx = 0; left_idx <= 3; left_idx++) {

         for (int right_idx = 0; right_idx <= 3; right_idx++) {

             Wave_vector_contravariant[left_idx]  +=      inv_Metric[left_idx][right_idx] * Wave_Vector_covariant[right_idx];
             B_field_covariant[left_idx]          += s_Metric.Metric[left_idx][right_idx] * B_field_contravariant[right_idx];
             Plasma_velocity_covavriant[left_idx] += s_Metric.Metric[left_idx][right_idx] * Plasma_velocity_contravariant[right_idx];
         }

     }

     /* --------------------- Compute the first two tetrad basis vectors --------------------- */

     // Wave_vec_dot_Plasma_vel is safe to divide by as it can never go to zero - Plasma_vel is never null.

     for (int index = 0; index <= 3; index++) {

         Tetrad[e_t_coord][index]   = Plasma_velocity_contravariant[index];
         Tetrad[e_phi_coord][index] = -Wave_vector_contravariant[index] / Wave_vec_dot_Plasma_vel - Plasma_velocity_contravariant[index];

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

         Tetrad[e_theta_coord][index] = (B_field_contravariant[index] + Plasma_vel_dot_B_field * Plasma_velocity_contravariant[index] - C_coeff * Tetrad[e_phi_coord][index]) / N_coeff;

         for (int i = 0; i <= 3; i++) {

             for (int j = 0; j <= 3; j++) {

                 for (int k = 0; k <= 3; k++) {

                     Tetrad[e_r_coord][index] += -Levi_Cevita_tensor[index][i][j][k] *
                                                      Plasma_velocity_covavriant[i] *
                                                      Wave_Vector_covariant[j] *
                                                      B_field_covariant[k] / Wave_vec_dot_Plasma_vel / N_coeff;

                 }

             }

         }

     }


     /* --------------------- Compute the inverse tetrad --------------------- */

     double Minkowski_Metric[4][4] = { {-1., 0., 0., 0.},
                                       { 0., 1., 0., 0.},
                                       { 0., 0., 1., 0.},
                                       { 0., 0., 0., 1.} };

     for (int left_idx = 0; left_idx <= 3; left_idx++) {

         for (int right_idx = 0; right_idx <= 3; right_idx++) {

             for (int m = 0; m <= 3; m++) {

                 for (int g = 0; g <= 3; g++) {

                     inv_Tetrad[left_idx][right_idx] += Minkowski_Metric[left_idx][m] * s_Metric.Metric[right_idx][g] * Tetrad[m][g];

                     if (isnan(inv_Tetrad[left_idx][right_idx]) || isinf(inv_Tetrad[left_idx][right_idx]) ||
                         isnan(Tetrad[left_idx][right_idx])     || isinf(Tetrad[left_idx][right_idx])) {

                         return ERROR;

                     }
                         
                 }
             }
         }
     }

     return OK;

}

void static Parallel_Transport_Polarization_Vector(double State_Vector[], 
                                                   Initial_conditions_type* const s_Initial_Conditions, 
                                                   std::complex<double> Polarization_Vector[]) {

    Metric_type s_Metric        = s_Initial_Conditions->Spacetimes[e_metric]->get_metric(State_Vector);
    Metric_type s_dr_Metric     = s_Initial_Conditions->Spacetimes[e_metric]->get_dr_metric(State_Vector);
    Metric_type s_dtheta_Metric = s_Initial_Conditions->Spacetimes[e_metric]->get_dtheta_metric(State_Vector);

    double inv_Metric[4][4]{};
    invert_metric(inv_Metric, s_Metric.Metric);
    
    double Connection_Coefficients[4][4][4]{};

    get_connection_coefficients(s_Metric, s_dr_Metric, s_dtheta_Metric, Connection_Coefficients);

    std::complex<double> Polarization_Vector_Derivative[4]{};

  /* ========================== Construct the full CONTRVARIANT photon wave vector ========================== */

    double& J = State_Vector[e_path_log_p_phi];

    double p_t_contravariant     = inv_Metric[e_t_coord][e_t_coord] * (-1) + inv_Metric[e_t_coord][e_phi_coord] * J;
    double p_r_contravariant     = inv_Metric[e_r_coord][e_r_coord] * State_Vector[e_path_log_p_r];
    double p_theta_contravariant = inv_Metric[e_theta_coord][e_theta_coord] * State_Vector[e_path_log_p_theta];
    double p_phi_contravariant   = inv_Metric[e_phi_coord][e_phi_coord] * J + inv_Metric[e_phi_coord][e_t_coord] * (-1);

    double photon_wave_vector[4] = {p_t_contravariant, p_r_contravariant, p_theta_contravariant, p_phi_contravariant};

  /* ========================== Compute the derivative of the polarization vector from the parallel transport ========================== */

    for (int derivative_index = 0; derivative_index <= STOKES_PARAM_NUM - 1; derivative_index++) {

        for (int polarization_index = 0; polarization_index <= STOKES_PARAM_NUM - 1; polarization_index++) {

            for (int wave_vector_index = 0; wave_vector_index <= STOKES_PARAM_NUM - 1; wave_vector_index++) {

                Polarization_Vector_Derivative[derivative_index] = -Connection_Coefficients[derivative_index][wave_vector_index][polarization_index] * photon_wave_vector[wave_vector_index] * Polarization_Vector[polarization_index];

            }

        }

    }

    for (int index = 0; index <= 3; index++) {

        Polarization_Vector[index] += Polarization_Vector_Derivative[index] * State_Vector[e_path_log_step];

    }

}

void static Map_Polarization_Vector_to_Stokes(const double inv_Stokes_Tetrad[4][4],
                                              std::complex<double> Coord_Basis_Pol_vec[STOKES_PARAM_NUM],
                                              double Stokes_Vector[STOKES_PARAM_NUM]) {

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

void static Map_Stokes_to_Polarization_Vector(const double Stokes_Vector[STOKES_PARAM_NUM],
                                              const double Stokes_Tetrad[4][4],
                                              std::complex<double> Coord_Basis_Pol_vec[STOKES_PARAM_NUM]) {

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

void static Propagate_forward_emission(Initial_conditions_type* const s_Initial_Conditions, Results_type* const s_Ray_results, int* const N_theta_turning_points) {

    int const Max_theta_turning_points = *N_theta_turning_points;
    *N_theta_turning_points = 0;

    int Image_Order = compute_image_order(Max_theta_turning_points, s_Initial_Conditions);

    double Stokes_Vector[STOKES_PARAM_NUM]{};

    std::complex<double> Coord_Basis_Pol_vec[4] = { std::complex<double>(0,0),
                                                    std::complex<double>(1,0), 
                                                    std::complex<double>(0,0), 
                                                    std::complex<double>(0,0) };

    double* Logged_ray_path[INTERPOLATION_NUM]{};

    // TODO: Propagate this aswell
    double Optical_Depth{};

    for (int log_index = s_Ray_results->Ray_log_struct.Log_length; log_index > 0; log_index--) {

        /* ================== Pick out the ray position / momenta from the Log, at the given log index ================== */

        Logged_ray_path[Current] = &(s_Ray_results->Ray_log_struct.Ray_path_log[ log_index      * e_path_log_number]);
        Logged_ray_path[Next]    = &(s_Ray_results->Ray_log_struct.Ray_path_log[(log_index - 1) * e_path_log_number]);

        *N_theta_turning_points += Increment_theta_turning_points(Logged_ray_path[Current], Logged_ray_path[Next]);

        double step = Logged_ray_path[Current][e_path_log_step] * MASS_TO_CM;

        if (INCLUDE_POLARIZATION) {

            /* ======================================== Parallel transport the polarization vector ======================================== */

            Parallel_Transport_Polarization_Vector(Logged_ray_path[Current], s_Initial_Conditions, Coord_Basis_Pol_vec);

            /* =========================================================================================================================== */
        }

        double Normalized_plasma_density = s_Initial_Conditions->OTT_model->get_disk_density_profile(Logged_ray_path[Current]);

        if (Normalized_plasma_density > 1e-6){

            double Tetrad[4][4]{};
            double inv_Tetrad[4][4]{};

            if (INCLUDE_POLARIZATION) {

                int Tetrad_OK = Construct_Stokes_Tetrad(Tetrad, inv_Tetrad, s_Initial_Conditions, Logged_ray_path[Current]);

                if (OK == Tetrad_OK) {

                    Map_Polarization_Vector_to_Stokes(std::as_const(inv_Tetrad), Coord_Basis_Pol_vec, Stokes_Vector);

                }

            }

            /* ======================================== Propagate thhe radiative transfer equations ======================================== */

            double* U_source_coord      = s_Initial_Conditions->OTT_model->get_disk_velocity(Logged_ray_path[Current], s_Initial_Conditions);
            double* U_source_coord_next = s_Initial_Conditions->OTT_model->get_disk_velocity(Logged_ray_path[Next], s_Initial_Conditions);

            double redshift[2]{};
            redshift[Current] = Redshift(Logged_ray_path[Current], U_source_coord);
            redshift[Next]    = Redshift(Logged_ray_path[Next], U_source_coord_next);

            double emission_functions[INTERPOLATION_NUM][STOKES_PARAM_NUM]{}, 
                    faradey_functions[INTERPOLATION_NUM][STOKES_PARAM_NUM]{}, 
                 absorbtion_functions[INTERPOLATION_NUM][STOKES_PARAM_NUM]{};

            if (!isinf(1.0 / redshift[Current]) && !isinf(1.0 / redshift[Next])) {

                s_Initial_Conditions->OTT_model->get_synchotron_transfer_functions(Logged_ray_path[Current], s_Initial_Conditions, emission_functions[Current], faradey_functions[Current], absorbtion_functions[Current]);
                s_Initial_Conditions->OTT_model->get_synchotron_transfer_functions(Logged_ray_path[Next], s_Initial_Conditions, emission_functions[Next], faradey_functions[Next], absorbtion_functions[Next]);

                /* ------ Account for the relativistic doppler effet via the redshift ------ */

                for (int index = 0; index <= 1; index++) {

                    for (int stokes_idx = 0; stokes_idx <= 3; stokes_idx++) {

                        emission_functions[index][stokes_idx]   *= redshift[index] * redshift[index];
                        faradey_functions[index][stokes_idx]    /= redshift[index];
                        absorbtion_functions[index][stokes_idx] /= redshift[index];

                    }

                }

                Propagate_Stokes_vector(Implicit_Trapezoid, emission_functions, absorbtion_functions, faradey_functions, step, Stokes_Vector);

            }

            if (INCLUDE_POLARIZATION){

                /* =================================== Convert the stokes vector into a polarization vector =================================== */

                Map_Stokes_to_Polarization_Vector(std::as_const(Stokes_Vector), std::as_const(Tetrad), Coord_Basis_Pol_vec);

            }

        }

        double Polarized_Intensity = sqrt(Stokes_Vector[Q] * Stokes_Vector[Q] +
                                          Stokes_Vector[U] * Stokes_Vector[U] +
                                          Stokes_Vector[V] * Stokes_Vector[V]);

        double fractional_polarization = Polarized_Intensity / Stokes_Vector[I];

        if (fractional_polarization > 1.001 && Stokes_Vector[I] > 1e-10) {

            exit(ERROR);

        }

        log_ray_emission(Stokes_Vector, Optical_Depth, s_Ray_results, log_index);

        Seperate_Image_into_orders(Max_theta_turning_points, *N_theta_turning_points, s_Initial_Conditions, s_Ray_results, Stokes_Vector);
         
    }

}

Results_type* Propagate_ray(Initial_conditions_type* s_Initial_Conditions) {

    /* ===============================================================================================
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
    ================================================================================================ */

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

    s_Ray_results.Parameters = s_Initial_Conditions->Spacetimes[e_metric]->get_parameters();

    // Initialize the State Vector
    double State_vector[] = {r_obs, theta_obs, phi_obs, J, p_theta_0, p_r_0};

    // Initialize a vector that stores the old state
    double Old_state[e_State_Number]{};

    // Set the old State Vector and the Test State Vector to the Initial State Vector
    for (int vector_indexer = e_r; vector_indexer <= e_p_r; vector_indexer += 1) {

        Old_state[vector_indexer] = State_vector[vector_indexer];

    }

    // Initialize counters for the Number Of Integration Steps and the Number Of Turning points of the Polar Coordinate
    int integration_count{}, N_theta_turning_points{};

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

            Evaluate_Equatorial_Disk(s_Initial_Conditions, &s_Ray_results, State_vector, Old_state, N_theta_turning_points);

           /* ============= Evaluate logical flags for terminating the integration ============= */

            if (controller.integration_complete || integration_count >= MAX_INTEGRATION_COUNT) { 

                if (integration_count >= MAX_INTEGRATION_COUNT) {

                    std::cout << "Max iterations reached!" << '\n';

                }

                s_Ray_results.Ray_log_struct.Log_length = integration_count;

              /* =========== Integrate the radiative transfer equations forward along the ray =========== */

                Propagate_forward_emission(s_Initial_Conditions, &s_Ray_results, &N_theta_turning_points);

              /* ======================================================================================== */

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