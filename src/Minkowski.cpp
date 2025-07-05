#include "Spacetimes.h"

Metric_type Minkowski_class::get_metric(const double* const State_Vector) const {

    Metric_type s_Minkowski_metric{};

    const double& r = State_Vector[e_r];
    const double sin_theta = sin(State_Vector[e_theta]);

    s_Minkowski_metric.Metric[e_t][e_t] = -1;
    s_Minkowski_metric.Metric[e_r][e_r] = 1;
    s_Minkowski_metric.Metric[e_theta][e_theta] = r * r;
    s_Minkowski_metric.Metric[e_phi][e_phi] = r * r * sin_theta * sin_theta;

    s_Minkowski_metric.Lapse_function = 1;
    s_Minkowski_metric.Shift_function = 0;

    return s_Minkowski_metric;

};

Metric_type Minkowski_class::get_dr_metric(const double* const State_Vector) const {

    Metric_type s_dr_Minkowski_metric{};

    const double& r = State_Vector[e_r];
    const double sin_theta = sin(State_Vector[e_theta]);

    s_dr_Minkowski_metric.Metric[e_t][e_t] = 0;
    s_dr_Minkowski_metric.Metric[e_r][e_r] = 0;
    s_dr_Minkowski_metric.Metric[e_theta][e_theta] = 2 * r;
    s_dr_Minkowski_metric.Metric[e_phi][e_phi] = 2 * r * sin_theta * sin_theta;

    s_dr_Minkowski_metric.Lapse_function = 0;
    s_dr_Minkowski_metric.Shift_function = 0;

    return s_dr_Minkowski_metric;
}

Metric_type Minkowski_class::get_dtheta_metric(const double* const State_Vector) const {

    Metric_type s_dtheta_Minkowski_metric{};

    const double& r = State_Vector[e_r];
    const double sin_theta = sin(State_Vector[e_theta]);
    const double cos_theta = cos(State_Vector[e_theta]);

    s_dtheta_Minkowski_metric.Metric[e_t][e_t] = 0;
    s_dtheta_Minkowski_metric.Metric[e_r][e_r] = 0;
    s_dtheta_Minkowski_metric.Metric[e_theta][e_theta] = 0;
    s_dtheta_Minkowski_metric.Metric[e_phi][e_phi] = 2 * r * r * sin_theta * cos_theta;

    s_dtheta_Minkowski_metric.Lapse_function = 0;
    s_dtheta_Minkowski_metric.Shift_function = 0;

    return s_dtheta_Minkowski_metric;
}

Metric_type Minkowski_class::get_d2r_metric(const double* const State_Vector) const {

    Metric_type s_d2r_Minkowski_metric{};

    const double& r = State_Vector[e_r];
    const double sin_theta = sin(State_Vector[e_theta]);

    s_d2r_Minkowski_metric.Metric[e_t][e_t] = 0;
    s_d2r_Minkowski_metric.Metric[e_r][e_r] = 0;
    s_d2r_Minkowski_metric.Metric[e_theta][e_theta] = 2;
    s_d2r_Minkowski_metric.Metric[e_phi][e_phi] = 2 * sin_theta * sin_theta;

    s_d2r_Minkowski_metric.Lapse_function = 0;
    s_d2r_Minkowski_metric.Shift_function = 0;

    return s_d2r_Minkowski_metric;
}

int Minkowski_class::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

    double& r_obs = p_Initial_Conditions->Observer_params.distance;
    double& theta_obs = p_Initial_Conditions->Observer_params.inclination;

    p_Initial_Conditions->Init_Momentum[e_t] = -1;
    p_Initial_Conditions->Init_Momentum[e_phi] = -J_data[photon] * sin(theta_obs);
    p_Initial_Conditions->Init_Momentum[e_theta] = p_theta_data[photon];

    double& J = p_Initial_Conditions->Init_Momentum[e_phi];
    double& p_theta = p_Initial_Conditions->Init_Momentum[e_theta];

    double rad_potential = 1 - J * J / (r_obs * sin(theta_obs)) / (r_obs * sin(theta_obs));

    p_Initial_Conditions->Init_Momentum[e_r] = sqrt(rad_potential);

    return OK;
}

void Minkowski_class::get_EOM(const double* const State_vector, double* const Derivatives) const {

    const double& r = State_vector[e_r];

    double sin_theta = sin(State_vector[e_theta]);
    double cos_theta = cos(State_vector[e_theta]);

    const double& p_t     = State_vector[e_p_t];
    const double& p_r     = State_vector[e_p_r];
    const double& p_theta = State_vector[e_p_theta];
    const double& p_phi   = State_vector[e_p_phi];

    Derivatives[e_t]     = -p_t;
    Derivatives[e_r]     = p_r;
    Derivatives[e_theta] = p_theta / r / r;
    Derivatives[e_phi]   = p_phi / (r * sin_theta) / (r * sin_theta);

    Derivatives[e_p_t]     = 0.0;
    Derivatives[e_p_phi]   = 0.0;

    Derivatives[e_p_theta] = cos_theta / (sin_theta * sin_theta * sin_theta) / (r * r) * p_phi * p_phi;
    Derivatives[e_p_r]     = (p_theta * p_theta + (p_phi / sin_theta) * (p_phi / sin_theta)) / (r * r * r);

}

bool Minkowski_class::terminate_integration(const double* const State_vector, const double* const Derivatives) {

    bool scatter = State_vector[e_r] > 30 && Derivatives[e_r] < 0;

    return scatter;

};

Return_Values Minkowski_class::load_parameters(const Metric_parameters_type* const Metric_Parameters) {

    return OK;

}

