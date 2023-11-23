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

void get_Radiative_Transfer(double State_Vector[], double Derivatives[], Initial_conditions_type* s_Initial_Conditions) {

    /************************************************************************************************
    |                                                                                               |
    |   @ Description: Evaluate the radiative transfer equations, set up as two first order ODEs.   |
    |     Intensity as a function of optical depth, and optical depth as a function of the          |
    |     affine parameter.                                                                         |
    |                                                                                               |
    |   @ Inputs:                                                                                   |
    |     * State_Vector: Pointer to an array that holds the ray / photon State Vector              |
    |     * Derivatives: Pointer to an array that holds the evaluation of the E.O.M                 |
    |                                                                                               |
    |   @ Ouput: None                                                                               |
    |                                                                                               |
    ************************************************************************************************/

    double* temp_State_Vector = State_Vector;

    if (e_metric == Wormhole) {

        temp_State_Vector[e_r] = sqrt(State_Vector[e_r] * State_Vector[e_r] + WH_R_THROAT * WH_R_THROAT);

    }

    /* Get Disk Cooridinate Velocity */

    double* U_source_coord = s_Initial_Conditions->OTT_model->get_disk_velocity(temp_State_Vector, s_Initial_Conditions);

    /* Get The Redshift */

    double redshift = Redshift(temp_State_Vector, U_source_coord);

    double Emission_function{}, Absorbtion_function{};

    switch (e_emission) {

    case Synchotron_exact:

        Emission_function   = s_Initial_Conditions->OTT_model->get_emission_function_synchotron_exact(temp_State_Vector, s_Initial_Conditions);
        Absorbtion_function = s_Initial_Conditions->OTT_model->get_absorbtion_function(Emission_function, temp_State_Vector, redshift, OBS_FREQUENCY_CGS / redshift);

        break;

    case Synchotron_phenomenological:

        Emission_function   = s_Initial_Conditions->OTT_model->get_emission_function_synchotron_phenomenological(temp_State_Vector, s_Initial_Conditions);
        Absorbtion_function = s_Initial_Conditions->OTT_model->get_absorbtion_function(Emission_function, temp_State_Vector, redshift, OBS_FREQUENCY_CGS / redshift);
    
        break;

    default:

        std::cout << "Unsupported emission model!" << "\n";

    }

    /* Fill in radiative transfer derivatives */

    *(Derivatives + e_Intensity    ) = -redshift * redshift * Emission_function * exp(-State_Vector[e_Optical_Depth]) * MASS_TO_CM;
    *(Derivatives + e_Optical_Depth) = -Absorbtion_function / redshift * MASS_TO_CM;

}

void RK45(double State_Vector[], double Derivatives[], Step_controller* controller, Initial_conditions_type* s_Initial_Conditions) {

    /*************************************************************************************************
    |                                                                                                |
    |   @ Description: Perfomrs one iteration of numerical integration, using the Dormand - Prince   |
    |     method, then updates the photon State Vector and the Step Controller properties.           |
    |                                                                                                |
    |   @ Inputs:                                                                                    |
    |     * State_Vector: Pointer to an array that holds the photon State Vector to be updated       |
    |     * Derivatives: Pointer to an array that holds the evaluation of the E.O.M                  |
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
    double inter_State_vector[RK45_size * e_State_Number]{};

    while (iteration <= RK45_size - 1) { //runs trough the EOM evaluations in-between t and t + step

        for (int vector_indexer = 0; vector_indexer <= e_State_Number - 1; vector_indexer += 1) { //runs trough the state vector components

            inter_State_vector[vector_indexer + iteration * e_State_Number] = State_Vector[vector_indexer];

            for (int derivative_indexer = 0; derivative_indexer <= iteration - 1; derivative_indexer += 1) { //runs trough tough the Dormand-Prince coeficients matrix and adds on the contributions from the derivatives at the points between t and t + step;

                inter_State_vector[vector_indexer + iteration * e_State_Number] += -controller->step * Coeff_deriv[iteration][derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_State_Number];

            }

        }

        s_Initial_Conditions->Spacetimes[e_metric]->get_EOM(&inter_State_vector[iteration * e_State_Number], &Derivatives[iteration * e_State_Number]);
        get_Radiative_Transfer(&inter_State_vector[iteration * e_State_Number], &Derivatives[iteration * e_State_Number], s_Initial_Conditions);

        iteration += 1;

        s_Initial_Conditions->Spacetimes[e_metric]->reset_eval_bitmask();

    }

    for (int vector_indexer = 0; vector_indexer <= e_State_Number - 1; vector_indexer += 1) {

        New_State_vector_O5[vector_indexer] = State_Vector[vector_indexer];
        New_State_vector_O4[vector_indexer] = State_Vector[vector_indexer];

        for (int derivative_indexer = 0; derivative_indexer <= iteration - 1; derivative_indexer += 1) {

            New_State_vector_O5[vector_indexer] += -controller->step * Coeff_sol[derivative_indexer]      * Derivatives[vector_indexer + derivative_indexer * e_State_Number];
            New_State_vector_O4[vector_indexer] += -controller->step * Coeff_test_sol[derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_State_Number];

        }

        state_error[vector_indexer] = New_State_vector_O5[vector_indexer] - New_State_vector_O4[vector_indexer];
       
    }

    controller->sec_prev_err = controller->prev_err;
    controller->prev_err     = controller->current_err;
    controller->current_err  = my_max(state_error, e_State_Number);
    controller->update_step();

    // For the Novikov-Thorne disk, depenging on where the equatorial crossing happens, the linear interpolation
    // of the crossing point might not be accurate enough, so halving the step is required.

    double z = State_Vector[e_r] * cos(State_Vector[e_theta]);
    bool near_NT_disk = z * z < 0.5 * 0.5 &&
                        State_Vector[e_r] * State_Vector[e_r] < s_Initial_Conditions->NT_model->get_r_out() * s_Initial_Conditions->NT_model->get_r_out() &&
                        State_Vector[e_r] * State_Vector[e_r] > s_Initial_Conditions->NT_model->get_r_in() *  s_Initial_Conditions->NT_model->get_r_in();

    if (Evaluate_NT_disk && near_NT_disk) {

        controller->step /= 2;

    }

    // For the Gauss-Bonnet Naked Singularity, certain photons scatter from very close to the singularity.
    // Generating the correct form of the images produced from those photons requires a finer step.
    // Otherwise the images look wonky.

    if (e_metric == Gauss_Bonnet && State_Vector[e_r] < 2.2) {

        controller->step /= 3;

    }

    // For the JNW Naked Singularity, certain photons scatter from very close to the singularity.
    // Close enough that it requires "manual" scattering, by flipping the p_r sign.
    // Otherwise the photons never reach the turning point and the integration grinds to a halt.

    if (e_metric == Naked_Singularity) {

        if (State_Vector[e_r] - JNW_R_SINGULARITY < 5e-8) {

            State_Vector[e_p_r] *= -1;

        }

    }

    if (controller->continue_integration) {

        for (int vector_indexer = 0; vector_indexer <= e_State_Number - 1; vector_indexer += 1) {

            State_Vector[vector_indexer] = New_State_vector_O5[vector_indexer];

        }

    }

}

Step_controller::Step_controller(double init_stepsize) {

    Gain_I =  0.58 / 5;
    Gain_P = -0.21 / 5;
    Gain_D =  0.10 / 5;

    step = init_stepsize;

    current_err     = RK45_ACCURACY;
    prev_err        = RK45_ACCURACY;
    sec_prev_err    = RK45_ACCURACY;

    continue_integration = false;

}

void Step_controller::update_step() {

    if (current_err < RK45_ACCURACY )
    {
        step = pow(RK45_ACCURACY / (current_err  + SAFETY_2), Gain_I) *
               pow(RK45_ACCURACY / (prev_err     + SAFETY_2), Gain_P) *
               pow(RK45_ACCURACY / (sec_prev_err + SAFETY_2), Gain_D) * step;

        continue_integration = true;

    }
    else
    {

        step = SAFETY_1 * step * pow(RK45_ACCURACY / (current_err + SAFETY_2), 0.25);

        continue_integration = false;

    }

}