#include "Lensing.h"

//! Runs one iteration of the Dormond - Prince adaptive integrator.
/*! Runs one iteration of the Dormond - Prince adaptive integrator, and updates the State Vector and Step Controller instance accordingly.
*
*   \param [out] State_Vector - Pointer to the array that holds the photon State Vector.
*   \param [out] Controller - Pointer to the Step Controller class instance.
*   \param [in] p_Sim_context - Pointer to the Simulation Context struct.
*   \return Nothing
*/
void RK45(double* const State_Vector, Step_controller* const controller, const Simulation_Context_type* const p_Sim_context) {

    // Initialize the iteration counter
    int iteration = 0;

    // Initialize the state errors.
    double state_error[e_State_Number]{};
    double state_rel_err[e_State_Number]{};

    // Initialize the array that holds the intermediate EOM RHS evaluations.
    double Derivatives[RK45_size * e_State_Number]{};

    // Initialize the array that holds the intermediate State Vectors.
    double inter_State_vector[RK45_size * e_State_Number]{};

    // Initialize the array that holds the two new solutions that the DP54 method computes.
    double New_State_vector_O5[e_State_Number]{};
    double New_State_vector_O4[e_State_Number]{};

    // Runs trough the EOM evaluations in-between t and t + step.
    while (iteration <= RK45_size - 1) { 

        // Runs trough the state vector components.
        for (int vector_indexer = 0; vector_indexer <= e_State_Number - 2; vector_indexer += 1) { 

            inter_State_vector[vector_indexer + iteration * e_State_Number] = State_Vector[vector_indexer];

            // Runs trough tough the Dormand-Prince coeficients matrix and adds on the contributions from the derivatives at the points between t and t + step.
            for (int derivative_indexer = 0; derivative_indexer <= iteration - 1; derivative_indexer += 1) { 

                inter_State_vector[vector_indexer + iteration * e_State_Number] += -controller->step * Coeff_deriv[iteration][derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_State_Number];

            }
        }

        p_Sim_context->p_Spacetime->get_EOM(&inter_State_vector[iteration * e_State_Number], &Derivatives[iteration * e_State_Number]);

        iteration += 1;

    }

    // Compute the new state vectors.
    for (int vector_indexer = 0; vector_indexer <= e_State_Number - 2; vector_indexer += 1) {

        New_State_vector_O5[vector_indexer] = State_Vector[vector_indexer];
        New_State_vector_O4[vector_indexer] = State_Vector[vector_indexer];

        for (int derivative_indexer = 0; derivative_indexer <= iteration - 1; derivative_indexer += 1) {

            New_State_vector_O5[vector_indexer] += -controller->step * Coeff_sol[derivative_indexer]      * Derivatives[vector_indexer + derivative_indexer * e_State_Number];
            New_State_vector_O4[vector_indexer] += -controller->step * Coeff_test_sol[derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_State_Number];

        }

        state_error[vector_indexer] = New_State_vector_O5[vector_indexer] - New_State_vector_O4[vector_indexer];
       
    }

    // The integrator might jump pass surfaces that are singular for the EOM (like the JNW singularity at 2 / gamma)
    // In this case the whole state vector becomes a NaN. I check for this and update the integration step by hand,
    // then set the continue_integration flag to "false" to force the integrator to redo the current iteration with a smaller step.
    if (isnan(New_State_vector_O5[e_r])) {

        controller->continue_integration = false;
        controller->step /= 10.0;

        return;

    }

    // Update the state errors
    controller->previous_step = controller->step;
    controller->sec_prev_err  = controller->prev_err;
    controller->prev_err      = controller->current_err;
    controller->current_err   = get_max_element(state_error, e_State_Number - 1);

    // Update the controller step
    controller->update_step(std::as_const(State_Vector));

    if (controller->continue_integration) {

        // Update the state vector
        for (int vector_indexer = 0; vector_indexer <= e_State_Number - 1; vector_indexer += 1) {

            State_Vector[vector_indexer] = New_State_vector_O5[vector_indexer];

        }

        controller->integration_complete = p_Sim_context->p_Spacetime->terminate_integration(New_State_vector_O5, Derivatives);

        // For the JNW Naked Singularity, certain photons scatter from very close to the singularity.
        // Close enough that it requires "manual" scattering, by flipping the p_r sign.
        // Otherwise the photons never reach the turning point and the integration grinds to a halt.
        if (p_Sim_context->p_Init_Conditions->Metric_params.e_Spacetime == Janis_Newman_Winicour && p_Sim_context->p_Init_Conditions->Metric_params.JNW_Gamma_Parameter < 0.5) {

            if (State_Vector[e_r] - 2 / p_Sim_context->p_Init_Conditions->Metric_params.JNW_Gamma_Parameter < 1e-8) {

                State_Vector[e_p_r] *= -1;

            }

        }
    }
}

Step_controller::Step_controller(const Integrator_parameters_type Integrator_parameters) {

    this->Controller_type = Integrator_parameters.Controller_type;

    this->Gustafsson_k_1 = Integrator_parameters.Gustafsson_k1;
    this->Gustafsson_k_2 = Integrator_parameters.Gustafsson_k2;

    this->Max_rel_step_increase = Integrator_parameters.Max_rel_step_increase;
    this->Min_rel_step_increase = Integrator_parameters.Min_rel_step_increase;

    this->Gain_I = Integrator_parameters.PID_gain_I;
    this->Gain_P = Integrator_parameters.PID_gain_P;
    this->Gain_D = Integrator_parameters.PID_gain_D;

    this->step = Integrator_parameters.Init_stepzie;
    this->previous_step = Integrator_parameters.Init_stepzie;

    this->Max_absolute_err = Integrator_parameters.RK_45_accuracy;

    this->current_err  = this->Max_absolute_err;
    this->prev_err     = this->Max_absolute_err;
    this->sec_prev_err = this->Max_absolute_err;

    this->Safety_1 = Integrator_parameters.Safety_1;
    this->Safety_2 = Integrator_parameters.Safety_2;

    this->Max_integration_count = Integrator_parameters.Max_integration_count;

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
    
    double Error_threshold = this->Max_absolute_err * (1 + get_max_element(State_Vector, e_State_Number - 1));

    switch (this->Controller_type) {

    case PID:

        Rel_step_increase = this->Safety_1 * pow(Error_threshold / (current_err + this->Safety_2), this->Gain_I) *
                                             pow(Error_threshold / (prev_err + this->Safety_2), this->Gain_P) *
                                             pow(Error_threshold / (sec_prev_err + this->Safety_2), this->Gain_D);

        Rel_step_increase = std::min(this->Max_rel_step_increase, std::max(this->Min_rel_step_increase, Rel_step_increase));

        this->step = Rel_step_increase * this->step;

        break;

    default:

        Rel_step_increase = this->Safety_1 * pow(Error_threshold / (current_err + this->Safety_2), this->Gustafsson_k_1) *
                                             pow(current_err / (prev_err + this->Safety_2), this->Gustafsson_k_2);

        Rel_step_increase = std::min(this->Max_rel_step_increase, std::max(this->Min_rel_step_increase, Rel_step_increase));

        this->step = Rel_step_increase * this->step;

        break;

    }

    if (current_err < Error_threshold)
    {
        this->continue_integration = true;
    }
    else
    {
        this->continue_integration = false;
    }
}