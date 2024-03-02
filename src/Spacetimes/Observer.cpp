#pragma once

#include <cmath>

#include "Structs.h"
#include "Spacetimes.h"
#include "Constants.h"

extern Spacetime_Base_Class* Spacetimes[];

Observer_class::Observer_class(double r, double theta, double phi) {

    r_obs     = r;
    theta_obs = theta;
    phi_obs   = phi;

    double Obs_State_Vector[3] = {r, theta, phi};

    Metric_type obs_metric = Spacetimes[e_metric]->get_metric(Obs_State_Vector);

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