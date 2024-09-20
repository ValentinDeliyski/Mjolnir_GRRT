#pragma once

#ifndef LENSING

    #define LENSING

    class Step_controller {

    public:

        Step_controller(double const init_stepsize);

        void update_step(double const State_Vector[]);

        double Gain_I;
        double Gain_P;
        double Gain_D;

        double step;
        double previous_step;

        double current_err;
        double prev_err;
        double sec_prev_err;

        bool continue_integration;
        bool integration_complete;

    };

    void RK45(double State_Vector[], Step_controller* controller, Spacetime_Base_Class* p_Spacetime);

    Results_type* Propagate_ray(Simulation_Context_type* s_Sim_Context);

#endif 