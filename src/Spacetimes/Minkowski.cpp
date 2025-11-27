#include "Spacetimes.h"

Minkowski_class::Minkowski_class(const Metric_parameters_type* const p_Metric_Parameters){

    if (isnan(p_Metric_Parameters->Scattering_radius) || isinf(p_Metric_Parameters->Scattering_radius) || p_Metric_Parameters->Scattering_radius < 0) {

        throw std::runtime_error(std::format("Invalid value for the scattering radius: {}", p_Metric_Parameters->Scattering_radius));

    }

    this->Scattering_radius = p_Metric_Parameters->Scattering_radius;

}

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

void Minkowski_class::get_EOM(const double* const State_vector, double* const Derivatives) {

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

bool Minkowski_class::terminate_integration(const double* const State_vector) {

    return State_vector[e_r] > this->Scattering_radius && State_vector[e_p_r] < 0;

};
