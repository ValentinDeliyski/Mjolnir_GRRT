#pragma once

#include <cmath>

#include "Structs.h"
#include "Spacetimes.h"
#include "Constants.h"

Observer_class::Observer_class(Initial_conditions_type* p_Init_Conditions) {

    r_obs     = p_Init_Conditions->init_Pos[e_r];
    theta_obs = p_Init_Conditions->init_Pos[e_theta];
    phi_obs   = p_Init_Conditions->init_Pos[e_phi];

    Metric_type obs_metric{};
    memcpy(obs_metric.Metric, p_Init_Conditions->init_metric, sizeof(obs_metric.Metric));
    obs_metric.Lapse_function = p_Init_Conditions->init_metric_Redshift_func;
    obs_metric.Shift_function = p_Init_Conditions->init_metric_Shitft_func;

    /*

    Obs_velocity is given in contravatiant components

    */

    obs_velocity[0] = 1.0 / obs_metric.Lapse_function;
    obs_velocity[1] = 0;
    obs_velocity[2] = 0;
    obs_velocity[3] = obs_metric.Shift_function / obs_metric.Lapse_function;

}

double Observer_class::get_r_obs()     { return r_obs; };
double Observer_class::get_theta_obs() { return theta_obs; };
double Observer_class::get_phi_obs()   { return phi_obs; };

int Observer_class::get_obs_velocity(double Obs_velocity[4]) {

    for (int index = 0; index <= 3; index++) {

        Obs_velocity[index] = obs_velocity[index];

    }

    return OK;

}