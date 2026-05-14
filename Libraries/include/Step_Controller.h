#pragma once
#include "Enumerations.h"
#include "Structs.h"

class Step_controller_class {

public:

    Step_controller_class(const Step_Controller_parameters_type Controller_parameters);

    //! Updates the integration step, based on the previous State Error estimates, and the current State Vector.
    /*! Updates the integration step, based on the previous State Error estimates, and the current State Vector.
     *   Currently the following step controllers are implemented. The reference is https://arxiv.org/pdf/1806.08693:
     *      1) PID controller
     *      2) Gustafsson controller
     *
     */
    void update_step(Integrator_enums e_Active_integrator, const double r);

    double get_max_step(const double r) const;

    template<typename Vec_type>
    void update_state_errors(const Vec_type* State_Vector, const Vec_type* State_Error_Vector, Integrator_enums e_Active_integrator, int State_size) {

        this->sec_prev_err = this->prev_err;
        this->prev_err = this->current_err;

        double Abs_tol{}, Rel_tol{};

        if (e_Active_integrator >= Radiative_only_integrators) {

            throw std::runtime_error("Wrong active integrator in Step_controller_class::update_state_errors()!");

        }

        switch (e_Active_integrator) {

        default:

            Abs_tol = this->Parameters.RK_abs_accuracy;
            Rel_tol = this->Parameters.RK_rel_accuracy;

            break;

        case ESDIRK54:

            Abs_tol = this->Parameters.ESDIRK54_abs_accuracy;
            Rel_tol = this->Parameters.ESDIRK54_rel_accuracy;

        }

        double Total_State_Error{};

        for (int idx = 0; idx < State_size; idx++) {

            Total_State_Error += std::pow(std::abs(State_Error_Vector[idx]) / (Abs_tol + std::abs(State_Vector[idx]) * Rel_tol), 2);

        }

        Total_State_Error /= (State_size - 1);

        this->current_err = std::sqrt(Total_State_Error) + this->Parameters.Safety_2;

    }

    Step_Controller_parameters_type Parameters;

    double step;
    double previous_step;

    double current_err;
    double prev_err;
    double sec_prev_err;

};