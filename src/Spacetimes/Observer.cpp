#include "Spacetimes.h"

Observer_class::Observer_class(Simulation_Context_type* p_Sim_Context) {

    // Copy the observer parameters into the class variable for the sake of convenience
    memcpy(&this->obs_params, &p_Sim_Context->p_Init_Conditions->Observer_params, sizeof(this->obs_params));

    double obs_position[4] = {0, this->obs_params.distance, this->obs_params.inclination, this->obs_params.azimuth };

    Metric_type s_init_Metric = p_Sim_Context->p_Spacetime->get_global_metric(obs_position);
    p_Sim_Context->p_Init_Conditions->Init_metric.Lapse_function = s_init_Metric.Lapse_function;
    p_Sim_Context->p_Init_Conditions->Init_metric.Shift_function = s_init_Metric.Shift_function;

    /*

    The velocities are given in contravatiant components

    */

    this->obs_velocity[e_t]     = 1.0 / s_init_Metric.Lapse_function;
    this->obs_velocity[e_r]     = 0;
    this->obs_velocity[e_theta] = 0;
    this->obs_velocity[e_phi]   = s_init_Metric.Shift_function / s_init_Metric.Lapse_function;

    this->fiducial_obs_velocity[e_t]     = obs_velocity[e_t];
    this->fiducial_obs_velocity[e_r]     = obs_velocity[e_r];
    this->fiducial_obs_velocity[e_theta] = obs_velocity[e_theta];
    this->fiducial_obs_velocity[e_phi]   = obs_velocity[e_phi];

}

Observer_parameters_type Observer_class::get_parameters() const { return this->obs_params; }

const double* const Observer_class::get_obs_velocity() const {

    return this->obs_velocity;

}

double* const Observer_class::get_fiducial_obs_velocity(Metric_type* p_Metric){

    this->fiducial_obs_velocity[e_t]     = 1.0 / p_Metric->Lapse_function;
    this->fiducial_obs_velocity[e_r]     = 0;
    this->fiducial_obs_velocity[e_theta] = 0;
    this->fiducial_obs_velocity[e_phi]   = p_Metric->Shift_function / p_Metric->Lapse_function;

    return this->fiducial_obs_velocity;
}
