#include "Lensing.h"

//! Runs one iteration of the Dormond - Prince adaptive integrator.
/*! Runs one iteration of the Dormond - Prince adaptive integrator, and updates the State Vector and Step Controller instance accordingly.
*
*   \param [out] State_Vector - Pointer to the array that holds the photon State Vector.
*   \param [out] p_Controller - Pointer to the Step Controller class instance.
*   \param [in] p_Sim_context - Pointer to the Simulation Context struct.
*   \return Nothing
*/
void RK45(double* const State_Vector, Step_controller* const p_Controller, const Simulation_Context_type* const p_Sim_context) {

    // Initialize the iteration counter
    int iteration = 0;

    // Initialize the state errors.
    double state_error[e_Dynamic_state_size]{};
    double state_rel_err[e_Dynamic_state_size]{};

    // Initialize the array that holds the intermediate EOM RHS evaluations.
    double Derivatives[RK45_size * e_Dynamic_state_size]{};

    // Initialize the array that holds the intermediate State Vectors.
    double inter_State_vector[RK45_size * e_Dynamic_state_size]{};

    // Initialize the array that holds the two new solutions that the DP54 method computes.
    double New_State_vector_O5[e_Dynamic_state_size]{};
    double New_State_vector_O4[e_Dynamic_state_size]{};

    // Runs trough the EOM evaluations in-between t and t + step.
    while (iteration <= RK45_size - 1) { 

        // Runs trough the state vector components.
        for (int vector_indexer = 0; vector_indexer <= e_Dynamic_state_size - 1; vector_indexer += 1) {

            inter_State_vector[vector_indexer + iteration * e_Dynamic_state_size] = State_Vector[vector_indexer];

            // Runs trough tough the Dormand-Prince coeficients matrix and adds on the contributions from the derivatives at the points between t and t + step.
            for (int derivative_indexer = 0; derivative_indexer <= iteration - 1; derivative_indexer += 1) { 

                inter_State_vector[vector_indexer + iteration * e_Dynamic_state_size] += -p_Controller->step * Coeff_deriv[iteration][derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_Dynamic_state_size];

            }
        }

        p_Sim_context->p_Spacetime->get_EOM(&inter_State_vector[iteration * e_Dynamic_state_size], &Derivatives[iteration * e_Dynamic_state_size]);

        iteration += 1;

    }

    // Compute the new state vectors.
    for (int vector_indexer = 0; vector_indexer <= e_Dynamic_state_size - 1; vector_indexer += 1) {

        New_State_vector_O5[vector_indexer] = State_Vector[vector_indexer];
        New_State_vector_O4[vector_indexer] = State_Vector[vector_indexer];

        for (int derivative_indexer = 0; derivative_indexer <= RK45_size - 1; derivative_indexer += 1) {

            New_State_vector_O5[vector_indexer] += -p_Controller->step * Coeff_sol[derivative_indexer]      * Derivatives[vector_indexer + derivative_indexer * e_Dynamic_state_size];
            New_State_vector_O4[vector_indexer] += -p_Controller->step * Coeff_test_sol[derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_Dynamic_state_size];

        }

        state_error[vector_indexer] = New_State_vector_O5[vector_indexer] - New_State_vector_O4[vector_indexer];
       
    }

    p_Controller->integration_complete = p_Sim_context->p_Spacetime->terminate_integration(State_Vector);

    // The integrator might jump pass surfaces that are singular for the EOM (like the JNW singularity at 2 / gamma)
    // In this case the whole state vector becomes a NaN. I check for this and update the integration step by hand,
    // then set the continue_integration flag to "false" to force the integrator to redo the current iteration with a smaller step.
    if (isnan(New_State_vector_O5[e_r])) {

        p_Controller->continue_integration = false;
        p_Controller->step /= 10.0;

        return;

    }

    // Update the state errors
    p_Controller->previous_step = p_Controller->step;
    p_Controller->sec_prev_err  = p_Controller->prev_err;
    p_Controller->prev_err      = p_Controller->current_err;
    p_Controller->current_err   = get_max_element(state_error, e_Dynamic_state_size);

    // Update the controller step
    p_Controller->update_step(std::as_const(State_Vector));

    if (p_Controller->step > 10) {

        //p_Controller->step = 10;

    }

    if (p_Controller->continue_integration) {

        // Update the state vector
        memcpy(State_Vector, New_State_vector_O5, e_Dynamic_state_size * sizeof(double));

        // Using the "previous_step" here, because the step was updated by the above call to "update_step()" and we want the step that got us to this point.
        State_Vector[e_step] = p_Controller->previous_step;
        State_Vector[e_affine_param] -= p_Controller->previous_step;

        // For the JNW Naked Singularity, certain photons scatter from very close to the singularity.
        // Close enough that it requires "manual" scattering, by flipping the p_r sign.
        // Otherwise the photons never reach the turning point and the integration grinds to a halt.
        if (p_Sim_context->p_Init_Conditions->Metric_parameters.e_Spacetime == Janis_Newman_Winicour && p_Sim_context->p_Init_Conditions->Metric_parameters.JNW_Gamma_Parameter < 0.5) {

            if (State_Vector[e_r] - 2.0 / p_Sim_context->p_Init_Conditions->Metric_parameters.JNW_Gamma_Parameter < p_Sim_context->p_Init_Conditions->Metric_parameters.Min_distance_to_singular_point) {

                State_Vector[e_p_r] *= -1.0;

            }

        }
    }
}

Step_controller::Step_controller(const Integrator_parameters_type Integrator_parameters) {

    this->Parameters = Integrator_parameters;

    this->step = Integrator_parameters.Init_stepzie;
    this->previous_step = Integrator_parameters.Init_stepzie;

    this->current_err  = Integrator_parameters.RK_45_accuracy;
    this->prev_err     = Integrator_parameters.RK_45_accuracy;
    this->sec_prev_err = Integrator_parameters.RK_45_accuracy;

    this->continue_integration = false;
    this->integration_complete = false;

}

//! Updates the integration step, based on the previous State Error estimates, and the current State Vector.
/*! Updates the integration step, based on the previous State Error estimates, and the current State Vector.
 *   Currently the following step controllers are implemented. The reference is https://arxiv.org/pdf/1806.08693:
 *      1) PID controller
 *      2) Gustafsson controller
 *
 *   \param [in] State_Vector - Pointer to the array that holds the photon State Vector.
 *   \return Nothing
 */
void Step_controller::update_step(const double* const State_Vector) {

    double Rel_step_increase{};
    
    double Error_threshold = this->Parameters.RK_45_accuracy;

    switch (this->Parameters.Controller_type) {

    case PID:

        Rel_step_increase = this->Parameters.Safety_1 * pow(Error_threshold / (this->current_err + this->Parameters.Safety_2), this->Parameters.PID_gain_I) *
                                                        pow(Error_threshold / (this->prev_err + this->Parameters.Safety_2), this->Parameters.PID_gain_P) *
                                                        pow(Error_threshold / (this->sec_prev_err + this->Parameters.Safety_2), this->Parameters.PID_gain_D);

        break;

    default:

        Rel_step_increase = this->Parameters.Safety_1 * pow(Error_threshold / (this->current_err + this->Parameters.Safety_2), this->Parameters.Gustafsson_k1) *
                                                        pow(Error_threshold / (this->prev_err + this->Parameters.Safety_2), this->Parameters.Gustafsson_k2);

        break;

    }

    Rel_step_increase = std::min(this->Parameters.Max_rel_step_increase, std::max(this->Parameters.Min_rel_step_increase, Rel_step_increase));

    this->step = Rel_step_increase * this->step;

    if (this->current_err < Error_threshold)
    {
        this->continue_integration = true;
    }
    else
    {
        this->continue_integration = false;
    }
}