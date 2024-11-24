#pragma once
#include "General_math_functions.h"
#include "Enumerations.h"
#include "Spacetimes.h"
#include "Constants.h"
#include "Structs.h"
#include <cmath>

class Step_controller {

public:

    Step_controller(const Integrator_parameters_type Integrator_parameters);

    //! Updates the integration step, based on the previous State Error estimates, and the current State Vector.
    /*! Updates the integration step, based on the previous State Error estimates, and the current State Vector.
     *   Currently the following step controllers are implemented. The reference is https://arxiv.org/pdf/1806.08693:
     *      1) PID controller
     *      2) Gustafsson controller
     *
     *   \param [in] State_Vector - Pointer to the array that holds the photon State Vector.
     *   \return Nothing
     */
    void update_step(const double* State_Vector);

    Step_controller_type_enums Controller_type;

    double Gustafsson_k_1;
    double Gustafsson_k_2;

    double Gain_I;
    double Gain_P;
    double Gain_D;

    double Max_rel_step_increase;
    double Min_rel_step_increase;

    double step;
    double previous_step;

    double current_err;
    double prev_err;
    double sec_prev_err;

    double Max_absolute_err;

    double Safety_1;
    double Safety_2;

    int Max_integration_count;

    bool continue_integration;
    bool integration_complete;

};

//! Runs one iteration of the Dormond - Prince adaptive integrator.
/*! Runs one iteration of the Dormond - Prince adaptive integrator, and updates the State Vector and Step Controller instance accordingly.
*
*   \param [out] State_Vector - Pointer to the array that holds the photon State Vector.
*   \param [out] Controller - Pointer to the Step Controller class instance.
*   \param [in] p_Sim_context - Pointer to the Simulation Context struct.
*   \return Nothing
*/
void RK45(double* const State_Vector, Step_controller* const Controller, const Simulation_Context_type* const p_Sim_context);

void Propagate_ray(const Simulation_Context_type* s_Sim_Context, Results_type* const Ray_results);
