#pragma once

#include "Enumerations.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "Disk_Models.h"
#include "General_GR_functions.h"
#include "General_math_functions.h"

#include <cmath>
#include <vector>

extern e_Spacetimes e_metric;
extern std::vector<c_Spacetime_Base*> Spacetimes;
extern Optically_Thin_Toroidal_Model OTT_Model;

int get_Radiative_Transfer(double State_Vector[], double Derivatives[], int iteration, double J) {

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
        Absorbtion_function = OTT_Model.get_absorbtion_fucntion(Emission_function, OBS_FREQUENCY_CGS / redshift, OTT_Model.get_disk_temperature(State_Vector));

        break;

    case Synchotron_phenomenological:

        Emission_function = OTT_Model.get_emission_fucntion_synchotron_phenomenological(State_Vector, J, Spacetimes);
        Absorbtion_function = OTT_Model.get_absorbtion_fucntion(Emission_function, OBS_FREQUENCY_CGS / redshift, OTT_Model.get_disk_temperature(State_Vector));
    
        break;

    }

    /* Fill in radiative transfer derivatives */

    Derivatives[e_Intensity     + iteration * e_State_Number] = -redshift * redshift * Emission_function * exp(-State_Vector[e_Optical_Depth]) * MASS_TO_CM;
    Derivatives[e_Optical_Depth + iteration * e_State_Number] = -Absorbtion_function / redshift * MASS_TO_CM;

    return  OK;

}

Return_Value_enums RK45(double State_Vector[], double Derivatives[], double* step, double J, bool* continue_integration) {

    int iteration = 0;
    bool EOM_Status = OK;

    double state_error[e_State_Number]{};
    double State_vector_test[e_State_Number]{};
    double inter_State_vector[RK45_size * e_State_Number]{};

    double integration_error{};

    while (iteration <= RK45_size - 1 && EOM_Status == OK) { //runs trough the EOM evaluations in-between t and t + step

        for (int vector_indexer = 0; vector_indexer <= e_State_Number - 1; vector_indexer += 1) { //runs trough the state vector components

            inter_State_vector[vector_indexer + iteration * e_State_Number] = State_Vector[vector_indexer];

            for (int derivative_indexer = 0; derivative_indexer <= iteration - 1; derivative_indexer += 1) { //runs trough tough the Dormand-Prince coeficients matrix and adds on the contributions from the derivatives at the points between t and t + step;

                inter_State_vector[vector_indexer + iteration * e_State_Number] += -*step * Coeff_deriv[iteration][derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_State_Number];

            }

        }

        EOM_Status = Spacetimes[e_metric]->get_EOM(inter_State_vector, J, Derivatives, iteration);

        //if (State_Vector[e_r + iteration * e_State_Number] > Spacetimes[e_metric]->get_ISCO(Prograde)) {

            
        double current_iteration[e_State_Number] = { inter_State_vector[e_r + iteration * e_State_Number],
                                                    inter_State_vector[e_theta + iteration * e_State_Number],
                                                    inter_State_vector[e_phi + iteration * e_State_Number],
                                                    inter_State_vector[e_phi_FD + iteration * e_State_Number],
                                                    inter_State_vector[e_p_theta + iteration * e_State_Number],
                                                    inter_State_vector[e_p_r + iteration * e_State_Number],
                                                    inter_State_vector[e_Intensity + iteration * e_State_Number],
                                                    inter_State_vector[e_Optical_Depth + iteration * e_State_Number] };

        get_Radiative_Transfer(current_iteration, Derivatives, iteration, J);



        //}

        iteration += 1;

    }

    for (int vector_indexer = 0; vector_indexer <= e_State_Number - 1; vector_indexer += 1) {

        State_vector_test[vector_indexer] = State_Vector[vector_indexer];

        for (int derivative_indexer = 0; derivative_indexer <= iteration - 1; derivative_indexer += 1) {

            State_Vector[vector_indexer] += -*step * Coeff_sol[derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_State_Number];
            State_vector_test[vector_indexer] += -*step * Coeff_test_sol[derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_State_Number];

        }

        state_error[vector_indexer] = State_Vector[vector_indexer] - State_vector_test[vector_indexer];

    }

    integration_error = my_max(state_error);


    if (integration_error < RK45_ACCURACY)
    {
        *step = SAFETY_1 * *step * pow(RK45_ACCURACY / (integration_error + SAFETY_2), 0.2);

        *continue_integration = true;

    }
    else
    {

        *step = SAFETY_1 * *step * pow(RK45_ACCURACY / (integration_error + SAFETY_2), 0.25);

        *continue_integration = false;

    }

    return OK;

}