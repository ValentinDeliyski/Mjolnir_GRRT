#pragma once

#ifndef LENSING

    #define LENSING

    class Step_controller {

    public:

        Step_controller(Integrator_parameters_type Integrator_parameters);

        void update_step(double const State_Vector[]);

        double Gain_I;
        double Gain_P;
        double Gain_D;

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

    void RK45(double State_Vector[], Step_controller* controller, Simulation_Context_type* p_Sim_context);

    void Propagate_ray(Simulation_Context_type* s_Sim_Context, Results_type* Ray_results);

#endif 