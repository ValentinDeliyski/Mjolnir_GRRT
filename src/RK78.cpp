#include "Lensing.h"
#include "General_GR_functions.h"

void RK78(double* const State_Vector, Step_controller* const p_Controller, const Simulation_Context_type* const p_Sim_context) {

    int iteration = 0;

    double state_error[e_Dynamic_state_size]{};

    double Derivatives[RK78_size * e_Dynamic_state_size]{};

    double inter_State_vector[e_Dynamic_state_size]{};

    double New_State_vector_O8[e_Dynamic_state_size]{};
    double New_State_vector_O9[e_Dynamic_state_size]{};

    while (iteration < RK78_size) { 

        memcpy(inter_State_vector, State_Vector, e_Dynamic_state_size * sizeof(double));

        for (int vector_indexer = 0; vector_indexer < e_Dynamic_state_size; vector_indexer++) {

            for (int derivative_indexer = 0; derivative_indexer < iteration; derivative_indexer++) { 

                inter_State_vector[vector_indexer] += -p_Controller->step * RK78_Coeff_deriv[iteration][derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_Dynamic_state_size];

            }
        }

        p_Sim_context->p_Spacetime->get_EOM(inter_State_vector, &Derivatives[iteration * e_Dynamic_state_size]);

        iteration += 1;

    }

    for (int vector_indexer = 0; vector_indexer < e_Dynamic_state_size; vector_indexer++) {

        New_State_vector_O8[vector_indexer] = State_Vector[vector_indexer];
        New_State_vector_O9[vector_indexer] = State_Vector[vector_indexer];

        for (int derivative_indexer = 0; derivative_indexer < RK78_size; derivative_indexer++) {

            New_State_vector_O8[vector_indexer] += -p_Controller->step * RK78_Coeff_sol[derivative_indexer]      * Derivatives[vector_indexer + derivative_indexer * e_Dynamic_state_size];
            New_State_vector_O9[vector_indexer] += -p_Controller->step * RK78_Coeff_test_sol[derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_Dynamic_state_size];

        }

        state_error[vector_indexer] = New_State_vector_O8[vector_indexer] - New_State_vector_O9[vector_indexer];
       
    }

    // The integrator might jump pass surfaces that are singular for the EOM (like the JNW singularity at 2 / gamma)
    // In this case the whole state vector becomes a NaN. I check for this and update the integration step by hand,
    // then set the continue_integration flag to "false" to force the integrator to redo the current iteration with a smaller step.
    if (isnan(New_State_vector_O8[e_r])) {

        p_Controller->continue_integration = false;
        p_Controller->step /= 10.0;

        return;

    }

    p_Controller->update_state_errors(New_State_vector_O8, state_error);
    p_Controller->update_step(New_State_vector_O8);

    if (p_Controller->continue_integration) {

        // Update the state vector
        memcpy(State_Vector, New_State_vector_O8, e_Dynamic_state_size * sizeof(double));

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

        p_Controller->integration_complete = p_Sim_context->p_Spacetime->terminate_integration(State_Vector);

    }

}

void Step_controller::update_state_errors(const double* State_Vector, const double* State_Error_Vector) {

    this->sec_prev_err = this->prev_err;
    this->prev_err     = this->current_err;

    double Error_scale[e_Dynamic_state_size]{};

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        Error_scale[idx] = this->Parameters.RK_78_accuracy + std::fabs(State_Vector[idx]) * this->Parameters.RK_78_accuracy;

    }

    double Total_State_Error{};

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        Total_State_Error += (State_Error_Vector[idx] / Error_scale[idx]) * (State_Error_Vector[idx] / Error_scale[idx]);

    }

    Total_State_Error /= e_Dynamic_state_size;

    this->current_err = std::sqrt(Total_State_Error) + this->Parameters.Safety_2;

}

Step_controller::Step_controller(const Integrator_parameters_type Integrator_parameters) {

    this->Parameters = Integrator_parameters;

    this->step = this->Parameters.Init_stepzie;
    this->previous_step = this->Parameters.Init_stepzie;

    this->current_err  = this->Parameters.Safety_2;
    this->prev_err     = this->Parameters.Safety_2;
    this->sec_prev_err = this->Parameters.Safety_2;

    this->continue_integration = false;
    this->integration_complete = false;

}

void Step_controller::update_step(const double* const State_Vector) {

    if (!this->Parameters.Use_adaptive_step) {

        this->continue_integration = true;

        return;

    }

    this->previous_step = this->step;

    double Rel_step_increase{};

    switch (this->Parameters.Controller_type) {

    case PID:

        Rel_step_increase = this->Parameters.Safety_1 * pow(this->current_err, this->Parameters.PID_gain_I) *
                                                        pow(this->prev_err, this->Parameters.PID_gain_P) *
                                                        pow(this->sec_prev_err, this->Parameters.PID_gain_D);

        break;

    default:

        Rel_step_increase = this->Parameters.Safety_1 * pow(this->current_err, this->Parameters.Gustafsson_k1) *
                                                        pow(this->current_err / this->prev_err, this->Parameters.Gustafsson_k2);

        break;

    }

    Rel_step_increase = std::min(this->Parameters.Max_rel_step_increase, std::max(this->Parameters.Min_rel_step_increase, Rel_step_increase));

    this->step *= Rel_step_increase;

    if (this->step > this->Parameters.Max_stepsize) { this->step = this->Parameters.Max_stepsize; };

    if (this->current_err < 1.0)
    {
        this->continue_integration = true;
    }
    else
    {
        this->continue_integration = false;
    }
}