#pragma once

#ifndef LENSING

    #define LENSING

    class Step_controller {

    public:

        Step_controller(const Integrator_parameters_type Integrator_parameters);

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

    void RK45(double* const State_Vector, Step_controller* const controller, const Simulation_Context_type* const p_Sim_context);

    void Propagate_ray(const Simulation_Context_type* s_Sim_Context, Results_type* const Ray_results);

#endif 