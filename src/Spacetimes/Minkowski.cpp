#pragma once
#include "Spacetimes.h"
#include <format>

Minkowski_class::Minkowski_class(const Metric_parameters_type* const p_Metric_Parameters){

    if (isnan(p_Metric_Parameters->Scattering_radius) || isinf(p_Metric_Parameters->Scattering_radius) || p_Metric_Parameters->Scattering_radius < 0) {

        throw std::runtime_error(std::format("Invalid value for the scattering radius: {}", p_Metric_Parameters->Scattering_radius));

    }

    this->Scattering_radius = p_Metric_Parameters->Scattering_radius;

}

Metric_type Minkowski_class::get_local_metric(const double* const Local_State_Vector) const {

    Metric_type s_Minkowski_metric{};

    const double& r = Local_State_Vector[e_r];
    const double sin_theta = sin(Local_State_Vector[e_theta]);

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_Minkowski_metric.Metric[e_t][e_t] = -1;
    s_Minkowski_metric.Metric[e_r][e_r] = 1;
    s_Minkowski_metric.Metric[e_theta][e_theta] = r * r;
    s_Minkowski_metric.Metric[e_phi][e_phi] = r * r * sin_theta * sin_theta;

    s_Minkowski_metric.Lapse_function = 1;
    s_Minkowski_metric.Shift_function = 0;

    return s_Minkowski_metric;

};

Metric_type Minkowski_class::get_global_metric(const double* const Global_State_Vector) const {

    return this->get_local_metric(Global_State_Vector);

};

Metric_type Minkowski_class::get_dr_local_metric(const double* const Local_State_Vector) const {

    Metric_type s_dr_Minkowski_metric{};

    const double& r = Local_State_Vector[e_r];
    const double sin_theta = sin(Local_State_Vector[e_theta]);

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_dr_Minkowski_metric.Metric[e_theta][e_theta] = 2 * r;
    s_dr_Minkowski_metric.Metric[e_phi][e_phi] = 2 * r * sin_theta * sin_theta;

    s_dr_Minkowski_metric.Lapse_function = 0;
    s_dr_Minkowski_metric.Shift_function = 0;

    return s_dr_Minkowski_metric;
}

Metric_type Minkowski_class::get_dr_global_metric(const double* const Global_State_Vector) const {

    return this->get_dr_local_metric(Global_State_Vector);

};

Metric_type Minkowski_class::get_dtheta_local_metric(const double* const Local_State_Vector) const {

    Metric_type s_dtheta_Minkowski_metric{};

    const double& r = Local_State_Vector[e_r];
    const double sin_theta = sin(Local_State_Vector[e_theta]);
    const double cos_theta = cos(Local_State_Vector[e_theta]);

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_dtheta_Minkowski_metric.Metric[e_phi][e_phi] = 2 * r * r * sin_theta * cos_theta;

    s_dtheta_Minkowski_metric.Lapse_function = 0;
    s_dtheta_Minkowski_metric.Shift_function = 0;

    return s_dtheta_Minkowski_metric;
}

Metric_type Minkowski_class::get_dtheta_global_metric(const double* const Global_State_Vector) const {

    return this->get_dtheta_local_metric(Global_State_Vector);

};

Metric_type Minkowski_class::get_d2r_local_metric(const double* const Local_State_Vector) const {

    Metric_type s_d2r_Minkowski_metric{};

    const double& r = Local_State_Vector[e_r];
    const double sin_theta = sin(Local_State_Vector[e_theta]);

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

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

void Minkowski_class::Convert_global_to_local_coords(const double* const State_Vector_Global, const double* const Global_Vec_to_Convert, double* Local_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

    switch (Entry_to_convert) {

    case e_Full_State_Vector:

        memcpy(Local_Vec_to_Convert, Global_Vec_to_Convert, e_Full_state_size * sizeof(double));
        break;

    case e_Contravariant_vector:

        memcpy(Local_Vec_to_Convert, Global_Vec_to_Convert, 4 * sizeof(double));
        break;

    case e_Covariant_vector:
        memcpy(Local_Vec_to_Convert, Global_Vec_to_Convert, 4 * sizeof(double));
        break;

    case e_Coordinates:
        memcpy(Local_Vec_to_Convert, Global_Vec_to_Convert, 4 * sizeof(double));
        break;

    default:

        throw std::runtime_error("Unsupported coordinate conversion type. Something Broke in Convert_global_to_local_coords()!");

    }

}

void Minkowski_class::Convert_local_to_global_coords(const double* const State_Vector_Local, const double* const Local_Vec_to_Convert, double* Global_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

    switch (Entry_to_convert) {

    case e_Full_State_Vector:

        memcpy(Global_Vec_to_Convert, Local_Vec_to_Convert, e_Full_state_size * sizeof(double));
        break;

    case e_Contravariant_vector:

        memcpy(Global_Vec_to_Convert, Local_Vec_to_Convert, 4 * sizeof(double));
        break;

    case e_Covariant_vector:
        memcpy(Global_Vec_to_Convert, Local_Vec_to_Convert, 4 * sizeof(double));
        break;

    case e_Coordinates:
        memcpy(Global_Vec_to_Convert, Local_Vec_to_Convert, 4 * sizeof(double));
        break;

    default:

        throw std::runtime_error("Unsupported coordinate conversion type. Something Broke in Convert_local_to_global_coords()!");

    }

}