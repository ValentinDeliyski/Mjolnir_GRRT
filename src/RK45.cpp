#pragma once

#include "Enumerations.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "Disk_Models.h"
#include "General_GR_functions.h"
#include "General_math_functions.h"

#include <cmath>
#include <vector>
#include "Lensing.h"

void RK45(double* const State_Vector, Step_controller* const controller, const Simulation_Context_type* const p_Sim_context) {

    /*************************************************************************************************
    |                                                                                                |
    |   @ Description: Perfomrs one iteration of numerical integration, using the Dormand - Prince   |
    |     method, then updates the photon State Vector and the Step Controller properties.           |
    |                                                                                                |
    |   @ Inputs:                                                                                    |
    |     * State_Vector: Pointer to an array that holds the photon State Vector to be updated       |
    |     * Derivatives: Pointer to an array that holds the evaluation of the E.O.M.                 |
    |     * controller: Pointer to class instance of the integrator step controller                  |
    |                                                                                                |
    |   @ Ouput: None                                                                                |
    |                                                                                                |
    *************************************************************************************************/

    int iteration = 0;

    double state_error[e_State_Number]{};
    double state_rel_err[e_State_Number]{};
    double New_State_vector_O5[e_State_Number]{};
    double New_State_vector_O4[e_State_Number]{};
    double Derivatives[RK45_size * e_State_Number]{};
    double inter_State_vector[RK45_size * e_State_Number]{};

    while (iteration <= RK45_size - 1) { //runs trough the EOM evaluations in-between t and t + step

        for (int vector_indexer = 0; vector_indexer <= e_State_Number - 2; vector_indexer += 1) { //runs trough the state vector components

            inter_State_vector[vector_indexer + iteration * e_State_Number] = State_Vector[vector_indexer];

            for (int derivative_indexer = 0; derivative_indexer <= iteration - 1; derivative_indexer += 1) { //runs trough tough the Dormand-Prince coeficients matrix and adds on the contributions from the derivatives at the points between t and t + step;

                inter_State_vector[vector_indexer + iteration * e_State_Number] += -controller->step * Coeff_deriv[iteration][derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_State_Number];

            }
        }

        p_Sim_context->p_Spacetime->get_EOM(&inter_State_vector[iteration * e_State_Number], &Derivatives[iteration * e_State_Number]);

        iteration += 1;

    }

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
    // In this case the whole state vector becomes a NaN. I check for this and update the step with some huge error 
    // to force the integrator to redo this interation with a smaller step untill it succeeds.

    if (isnan(New_State_vector_O5[e_r])) {

        controller->update_step(std::as_const(State_Vector));

        return;

    }

    controller->previous_step = controller->step;
    controller->sec_prev_err  = controller->prev_err;
    controller->prev_err      = controller->current_err;

    controller->current_err = my_max(state_error, e_State_Number);

    controller->update_step(std::as_const(State_Vector));

    if (controller->continue_integration) {

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

    this->Gain_I = Integrator_parameters.Gain_I;
    this->Gain_P = Integrator_parameters.Gain_P;
    this->Gain_D = Integrator_parameters.Gain_D;

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

void Step_controller::update_step(const double const* State_Vector) {

    double const Error_threshold = this->Max_absolute_err + this->Max_absolute_err * my_max(State_Vector, e_State_Number);

    if (current_err < Error_threshold)
    {
        this->step = pow(Error_threshold / (current_err  + this->Safety_2), this->Gain_I) *
                     pow(Error_threshold / (prev_err     + this->Safety_2), this->Gain_P) *
                     pow(Error_threshold / (sec_prev_err + this->Safety_2), this->Gain_D) * this->step;

        this->continue_integration = true;
    }
    else
    {

        step = this->Safety_1 * this->step * pow(Error_threshold / (current_err + this->Safety_2), 0.25);

        this->continue_integration = false;
    }
}