#include "Spacetimes.h"

Observer_class::Observer_class(Simulation_Context_type* p_Sim_Context) {

    // Copy the observer parameters into the class variable for the sake of convenience
    memcpy(&this->obs_params, &p_Sim_Context->p_Init_Conditions->Observer_params, sizeof(this->obs_params));

    double obs_position[4] = {0, this->obs_params.distance, this->obs_params.inclination, this->obs_params.azimuth };

    Metric_type s_init_Metric = p_Sim_Context->p_Spacetime->get_metric(obs_position);
    p_Sim_Context->p_Init_Conditions->init_metric_Redshift_func = s_init_Metric.Lapse_function;
    p_Sim_Context->p_Init_Conditions->init_metric_Shitft_func   = s_init_Metric.Shift_function;

    /*

    Obs_velocity is given in contravatiant components

    */

    obs_velocity[0] = 1.0 / s_init_Metric.Lapse_function;
    obs_velocity[1] = 0;
    obs_velocity[2] = 0;
    obs_velocity[3] = s_init_Metric.Shift_function / s_init_Metric.Lapse_function;

}

Observer_parameters_type Observer_class::get_parameters() { return this->obs_params; }

double* Observer_class::get_obs_velocity() {

    return this->obs_velocity;

}