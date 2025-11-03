#pragma once
#include "Enumerations.h"
#include "Spacetimes.h"
#include "General_GR_functions.h"
#include <complex>

template<typename Vec_type>
void static Parallel_Transport_RHS(const double* const State_Vector,
                                   Vec_type* Vector_to_transport,
                                   Spacetime_Base_Class* const p_Spacetime,
                                   Tensor_type_enums e_Vec_type,
                                   Vec_type* const Vector_derivative) {

    Metric_type s_Metric = p_Spacetime->get_metric(State_Vector);
    Metric_type s_dr_Metric = p_Spacetime->get_dr_metric(State_Vector);
    Metric_type s_dtheta_Metric = p_Spacetime->get_dtheta_metric(State_Vector);

    double inv_Metric[4][4]{};
    invert_metric(inv_Metric, s_Metric.Metric);

    double Connection_Coefficients[4][4][4]{};
    get_connection_coefficients(s_Metric, s_dr_Metric, s_dtheta_Metric, Connection_Coefficients);

    double Photon_momentum_contravariant[4]{};
    Manipulate_index(&s_Metric, State_Vector + e_p_t, Photon_momentum_contravariant, Raise_index);

    /* ========================== Compute the derivative of the polarization vector from the parallel transport ========================== */

    for (int derivative_index = 0; derivative_index < e_Stokes_param_num; derivative_index++) {

        for (int polarization_index = 0; polarization_index < e_Stokes_param_num; polarization_index++) {

            for (int wave_vector_index = 0; wave_vector_index < e_Stokes_param_num; wave_vector_index++) {

                if (e_Vec_type == Contravariant) {

                    Vector_derivative[derivative_index] += -Connection_Coefficients[derivative_index][wave_vector_index][polarization_index] * Photon_momentum_contravariant[wave_vector_index] * Vector_to_transport[polarization_index];

                }
                else {

                    Vector_derivative[derivative_index] += Connection_Coefficients[polarization_index][derivative_index][wave_vector_index] * Photon_momentum_contravariant[wave_vector_index] * Vector_to_transport[polarization_index];

                }

            }

        }

    }

}

template<typename Vec_type>
void static Parallel_Transport_Vector(const double* const State_Vector,
                                      Spacetime_Base_Class* const p_Spacetime,
                                      Tensor_type_enums const e_Vec_type,
                                      Vec_type* const Vector_to_transport) {

    Vec_type RHS[Nyström_size * e_Stokes_param_num]{};
    double EOM[Nyström_size * e_Dynamic_state_size]{};

    Vec_type Temp_Vector[e_Stokes_param_num]{};
    double Temp_State_Vector[e_Dynamic_state_size]{};

    for (int RK5_stage = 0; RK5_stage < Nyström_size; RK5_stage++) {

        memcpy(Temp_Vector, Vector_to_transport, e_Stokes_param_num * sizeof(Vec_type));
        memcpy(Temp_State_Vector, State_Vector, e_Dynamic_state_size * sizeof(double));

        for (int derivative_indexer = 0; derivative_indexer < RK5_stage; derivative_indexer++) {

            for (int idx = 0; idx < 4; idx++) {

                Temp_Vector[idx] += Nyström_Deriv_coeffs[RK5_stage][derivative_indexer] * RHS[idx + derivative_indexer * e_Stokes_param_num] * State_Vector[e_step];

            }

            for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

                Temp_State_Vector[idx] += Nyström_Deriv_coeffs[RK5_stage][derivative_indexer] * EOM[idx + derivative_indexer * e_Dynamic_state_size] * State_Vector[e_step];

            }
        }

        Parallel_Transport_RHS(Temp_State_Vector, Temp_Vector, p_Spacetime, e_Vec_type, RHS + RK5_stage * e_Stokes_param_num);
        p_Spacetime->get_EOM(Temp_State_Vector, EOM + RK5_stage * e_Dynamic_state_size);

    } 

    for (int idx = 0; idx < e_Stokes_param_num; idx++) {

        for (int deriv_idx = 0; deriv_idx < Nyström_size; deriv_idx++) {

            Vector_to_transport[idx] += State_Vector[e_step] * Nyström_Coeff_sol[deriv_idx] * RHS[idx + deriv_idx * e_Stokes_param_num];

        }
    }
}