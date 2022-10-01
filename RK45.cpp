#pragma once

#include "Enumerations.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "Disk_Models.h"
#include "General_functions.h"

#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>
#include <vector>

extern e_Spacetimes e_metric;

Return_Value_enums RK45(double State_Vector[], double Derivatives[], double* step, double J, bool* continue_integration,
                        c_Observer Observer_class, Optically_Thin_Toroidal_Model OTT_Model, std::vector<c_Spacetime_Base*> Spacetimes) {

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


        get_Radiative_Transfer(inter_State_vector, Derivatives, iteration, J, OTT_Model, Observer_class, Spacetimes);


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

        double z = State_Vector[e_r] * cos(State_Vector[e_theta]);

        bool near_disk = z * z < 0.1 * 0.1 &&
            State_Vector[e_r] * State_Vector[e_r] < 1.1 * 1.1 * 50 * 50;

        if (near_disk)
        {
            *step = SAFETY_1 / 10 * *step * pow(RK45_ACCURACY / (integration_error + SAFETY_2), 0.2);
        }

        *continue_integration = true;

    }
    else
    {

        *step = SAFETY_1 * *step * pow(RK45_ACCURACY / (integration_error + SAFETY_2), 0.25);

        *continue_integration = false;

    }

    return OK;

}