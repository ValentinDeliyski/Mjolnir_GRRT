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

extern Spacetime_Base_Class* Spacetimes[];
extern Optically_Thin_Toroidal_Model OTT_Model;

void get_Radiative_Transfer(double State_Vector[], double Derivatives[], int iteration) {

    /************************************************************************************************
    |                                                                                               |
    |   @ Description: Evaluate the radiative transfer equations, set up as two first order ODEs.   |
    |     Intensity as a function of optical depth, and optical depth as a function of the          |
    |     affine parameter.                                                                         |
    |                                                                                               |
    |   @ Inputs:                                                                                   |
    |     * State_Vector: Pointer to an array that holds the ray / photon State Vector              |
    |     * Derivatives: Pointer to an array that holds the evaluation of the E.O.M                 |
    |     * iteration: The current RK45 iteration (0 to 7 - 1)                                      |
    |                                                                                               |
    |   @ Ouput: None                                                                               |
    |                                                                                               |
    ************************************************************************************************/

    double temp_State_Vector[e_State_Number]{};

    for (int index = e_r; index < e_State_Number - 1; index++) {

        temp_State_Vector[index] = State_Vector[index + iteration * e_State_Number];

    }

    if (e_metric == Wormhole) {

        temp_State_Vector[e_r] = sqrt(State_Vector[e_r + iteration * e_State_Number] * State_Vector[e_r + iteration * e_State_Number] + WH_R_THROAT * WH_R_THROAT);

    }

    /* Get Disk Cooridinate Velocity */

    double* U_source_coord = OTT_Model.get_disk_velocity(temp_State_Vector);

    /* Get The Redshift */

    double redshift = Redshift(temp_State_Vector, U_source_coord);

    double Emission_function{}, Absorbtion_function{};

    switch (e_emission) {

    case Synchotron_exact:

        Emission_function = OTT_Model.get_emission_function_synchotron_exact(temp_State_Vector);
        Absorbtion_function = OTT_Model.get_absorbtion_function(Emission_function, temp_State_Vector, redshift, OBS_FREQUENCY_CGS / redshift);

        break;

    case Synchotron_phenomenological:

        Emission_function = OTT_Model.get_emission_function_synchotron_phenomenological(temp_State_Vector);
        Absorbtion_function = OTT_Model.get_absorbtion_function(Emission_function, temp_State_Vector, redshift, OBS_FREQUENCY_CGS / redshift);
    
        break;

    }

    /* Fill in radiative transfer derivatives */

    Derivatives[e_Intensity     + iteration * e_State_Number] = -redshift * redshift * Emission_function * exp(-State_Vector[e_Optical_Depth]) * MASS_TO_CM;
    Derivatives[e_Optical_Depth + iteration * e_State_Number] = -Absorbtion_function / redshift * MASS_TO_CM;

}

void RK45(double State_Vector[], double Derivatives[], Step_controller* controller) {

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

        Spacetimes[e_metric]->get_EOM(inter_State_vector, Derivatives, iteration);
        get_Radiative_Transfer(inter_State_vector, Derivatives, iteration);

        iteration += 1;

        Spacetimes[e_metric]->reset_eval_bitmask();

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

    controller->current_err = my_max(state_error, e_State_Number);
    controller->update_step();

    //bool near_NT_disk = State_Vector[e_r] * State_Vector[e_r] * cos(State_Vector[e_theta]) * cos(State_Vector[e_theta]) < 0.5 * 0.5 &&
    //                    State_Vector[e_r] * State_Vector[e_r] < NT_Model.get_r_out() * NT_Model.get_r_out();

    //if (near_NT_disk) {

    //    controller->step /= 2;

    //}

    if (e_metric == Gauss_Bonnet && State_Vector[e_r] < 2.2) {

        controller->step /= 3;

    }

    if (controller->continue_integration) {

        for (int vector_indexer = 0; vector_indexer <= e_State_Number - 1; vector_indexer += 1) {

            State_Vector[vector_indexer] = New_State_vector_O5[vector_indexer];

        }

    }
    
    controller->sec_prev_err = controller->prev_err;
    controller->prev_err     = controller->current_err;

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