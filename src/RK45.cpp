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

extern std::vector<c_Spacetime_Base*> Spacetimes;
extern Novikov_Thorne_Model NT_Model;
extern Optically_Thin_Toroidal_Model OTT_Model;

void get_Radiative_Transfer(double State_Vector[], double Derivatives[], int iteration, double J) {

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
    |     * J: Covariant, azimuthal component of the ray / photon 4 - momentum                      |
    |                                                                                               |
    |   @ Ouput: None                                                                               |
    |                                                                                               |
    ************************************************************************************************/

    double r = State_Vector[e_r];

    if (e_metric == Wormhole) {

        r = sqrt(State_Vector[e_r] * State_Vector[e_r] + WH_R_THROAT * WH_R_THROAT);

    }

    /* Get Disk Cooridinate Velocity */

    double U_source_coord[4]{};

    OTT_Model.get_disk_velocity(U_source_coord, State_Vector, Spacetimes);

    /* Get The Redshift */

    double redshift = Redshift(J, State_Vector, U_source_coord);

    double Emission_function{}, Absorbtion_function{};

    switch (e_emission) {

    case Synchotron_exact:

        Emission_function = OTT_Model.get_emission_fucntion_synchotron_exact(State_Vector, J, Spacetimes);
        Absorbtion_function = OTT_Model.get_absorbtion_fucntion(Emission_function, State_Vector, redshift, OBS_FREQUENCY_CGS / redshift, OTT_Model.get_disk_temperature(State_Vector));

        break;

    case Synchotron_phenomenological:

        Emission_function = OTT_Model.get_emission_fucntion_synchotron_phenomenological(State_Vector, J, Spacetimes);
        Absorbtion_function = OTT_Model.get_absorbtion_fucntion(Emission_function, State_Vector, redshift, OBS_FREQUENCY_CGS / redshift, OTT_Model.get_disk_temperature(State_Vector));
    
        break;

    }

    /* Fill in radiative transfer derivatives */

    Derivatives[e_Intensity     + iteration * e_State_Number] = -redshift * redshift * Emission_function * exp(-State_Vector[e_Optical_Depth]) * MASS_TO_CM;
    Derivatives[e_Optical_Depth + iteration * e_State_Number] = -Absorbtion_function / redshift * MASS_TO_CM;

}

void RK45(double State_Vector[], double Derivatives[], double J, Step_controller* controller) {

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

        Spacetimes[e_metric]->get_EOM(inter_State_vector, J, Derivatives, iteration);
            
        double current_iteration[e_State_Number] = { inter_State_vector[e_r             + iteration * e_State_Number],
                                                     inter_State_vector[e_theta         + iteration * e_State_Number],
                                                     inter_State_vector[e_phi           + iteration * e_State_Number],
                                                     inter_State_vector[e_phi_FD        + iteration * e_State_Number],
                                                     inter_State_vector[e_p_theta       + iteration * e_State_Number],
                                                     inter_State_vector[e_p_r           + iteration * e_State_Number],
                                                     inter_State_vector[e_Intensity     + iteration * e_State_Number],
                                                     inter_State_vector[e_Optical_Depth + iteration * e_State_Number] };

        get_Radiative_Transfer(current_iteration, Derivatives, iteration, J);

        iteration += 1;

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

    controller->current_err = my_max(state_error);

    controller->update_step();

    bool near_NT_disk = State_Vector[e_r] * State_Vector[e_r] * cos(State_Vector[e_theta]) * cos(State_Vector[e_theta]) < 0.5 * 0.5 &&
                        State_Vector[e_r] * State_Vector[e_r] < NT_Model.get_r_out() * NT_Model.get_r_out();

    if (near_NT_disk) {

        controller->step /= 2;

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

    current_err  = RK45_ACCURACY;
    prev_err     = RK45_ACCURACY;
    sec_prev_err = RK45_ACCURACY;

    continue_integration = false;

}

void Step_controller::update_step() {

    if (current_err < RK45_ACCURACY)
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