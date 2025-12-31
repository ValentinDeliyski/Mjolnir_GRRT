#include "Spacetimes.h"

RBH_class::RBH_class(const Metric_parameters_type* const p_Metric_Parameters){

    if (isnan(p_Metric_Parameters->RBH_Parameter) || isinf(p_Metric_Parameters->RBH_Parameter) || p_Metric_Parameters->RBH_Parameter < 0) {

        throw std::runtime_error(std::format("Invalid value for the metric parameter: ", p_Metric_Parameters->RBH_Parameter));

    }

    if (isnan(p_Metric_Parameters->Scattering_radius) || isinf(p_Metric_Parameters->Scattering_radius) || p_Metric_Parameters->Scattering_radius < 0) {

        throw std::runtime_error(std::format("Invalid value for the scattering radius: ", p_Metric_Parameters->Scattering_radius));

    }

    if (isnan(p_Metric_Parameters->Min_distance_to_singular_point) || isinf(p_Metric_Parameters->Min_distance_to_singular_point) || p_Metric_Parameters->Min_distance_to_singular_point < 0) {

        throw std::runtime_error(std::format("Invalid value for the distance to the singular point: ", p_Metric_Parameters->Min_distance_to_singular_point));

    }

    this->Mass = p_Metric_Parameters->Mass;

    this->Horizon_radius = 0;

    if (this->Parameter > 2 * this->Mass) {

        this->Horizon_radius = sqrt(4 * this->Mass * this->Mass - this->Parameter * this->Parameter);

    }

    this->Parameter = p_Metric_Parameters->RBH_Parameter;
    this->Scattering_radius = p_Metric_Parameters->Scattering_radius;
    this->Min_distance_to_singular_point = p_Metric_Parameters->Min_distance_to_singular_point;

}

double* RBH_class::get_ISCO() {

    static double r_ISCO[2]{};

    r_ISCO[Inner] = sqrt(36 * this->Mass * this->Mass - this->Parameter * this->Parameter);
    r_ISCO[Outer] = r_ISCO[Inner];

    return r_ISCO;

}

double* RBH_class::get_Photon_Sphere() {

    double M = this->Mass;

    static double photon_orbit[2]{};

    photon_orbit[Inner] = sqrt(9 * M * M - this->Parameter * this->Parameter);
    photon_orbit[Outer] = photon_orbit[Inner];

    return photon_orbit;

}

Metric_type RBH_class::get_local_metric(const double* const Local_State_Vector) const {

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double rho = sqrt(r2 + this->Parameter * this->Parameter);

    Metric_type s_Metric{};

    /* --- Only the non-zero components are explicitly evaluated. --- */

    s_Metric.Metric[e_t][e_t]         = -(1 - 2 * this->Mass / rho);
    s_Metric.Metric[e_r][e_r]         = -1.0 / s_Metric.Metric[e_t][e_t];
    s_Metric.Metric[e_theta][e_theta] = rho * rho;
    s_Metric.Metric[e_phi][e_phi]     = s_Metric.Metric[e_theta][e_theta] * sin_theta * sin_theta;

    s_Metric.Lapse_function = sqrt(-s_Metric.Metric[e_t][e_t]);

    return s_Metric;
}

Metric_type RBH_class::get_global_metric(const double* const Global_State_Vector) const {

    return this->get_local_metric(Global_State_Vector);

}

Metric_type RBH_class::get_dr_local_metric(const double* const Local_State_Vector) const {

    Metric_type s_Metric = this->get_local_metric(Local_State_Vector);

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho = sqrt(r2 + this->Parameter * this->Parameter);
    double rho3 = rho * rho * rho;

    Metric_type s_dr_Metric{};

    /* --- Only the non-zero components are explicitly evaluated. --- */

    s_dr_Metric.Metric[e_t][e_t]         = -2 * this->Mass * r / rho3;
    s_dr_Metric.Metric[e_r][e_r]         = 1.0 / (s_Metric.Metric[e_t][e_t] * s_Metric.Metric[e_t][e_t]) * s_dr_Metric.Metric[e_t][e_t];
    s_dr_Metric.Metric[e_theta][e_theta] = 2 * r;
    s_dr_Metric.Metric[e_phi][e_phi]     = 2 * r * sin_theta * sin_theta;

    return s_dr_Metric;
}

Metric_type RBH_class::get_dr_global_metric(const double* const Global_State_Vector) const {

    return this->get_dr_local_metric(Global_State_Vector);

}

Metric_type RBH_class::get_dtheta_local_metric(const double* const Local_State_Vector) const {

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    Metric_type s_dtheta_Metric{};

    /* --- Only the non-zero components are explicitly evaluated. --- */

    s_dtheta_Metric.Metric[e_phi][e_phi] = 2 * r * r * sin_theta * cos_theta;

    return s_dtheta_Metric;
}


Metric_type RBH_class::get_dtheta_global_metric(const double* const Global_State_Vector) const {

    return this->get_dtheta_local_metric(Global_State_Vector);

}

Metric_type RBH_class::get_d2r_local_metric(const double* const Local_State_Vector) const {

    Metric_type s_Metric = this->get_local_metric(Local_State_Vector);
    Metric_type s_dr_Metric = this->get_dr_local_metric(Local_State_Vector);

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho = sqrt(r2 + this->Parameter * this->Parameter);
    double rho3 = rho * rho * rho;
    double rho5 = rho * rho * rho * rho * rho;

    Metric_type s_d2r_Metric{};

    /* --- Only the non-zero components are explicitly evaluated. --- */

    s_d2r_Metric.Metric[e_t][e_t]         = -2 * this->Mass / rho3 + 6 * this->Mass * r2 / (rho5);
    s_d2r_Metric.Metric[e_r][e_r]         = 1.0 / (s_Metric.Metric[e_t][e_t] * s_Metric.Metric[e_t][e_t]) * s_d2r_Metric.Metric[e_t][e_t] - 
                                            2.0 / (s_Metric.Metric[e_t][e_t] * s_Metric.Metric[e_t][e_t] * s_Metric.Metric[e_t][e_t]) * s_dr_Metric.Metric[e_t][e_t] * s_dr_Metric.Metric[e_t][e_t];
    s_d2r_Metric.Metric[e_theta][e_theta] = 2.0;
    s_d2r_Metric.Metric[e_phi][e_phi]     = 2 * sin_theta * sin_theta;

    return s_d2r_Metric;
}

void RBH_class::get_EOM(const double* const State_vector, double* const Derivatives) {

    double r = State_vector[e_r];
    double rho = sqrt(r * r + this->Parameter * this->Parameter);

    const double& J = State_vector[e_p_phi];

    double sin1 = sin(State_vector[e_theta]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(State_vector[e_theta]);
    double cos2 = cos1 * cos1;

    Derivatives[e_p_t] = - 1 / (1 - 2 * this->Mass / rho) * State_vector[e_p_t];
    Derivatives[e_r] = (1 - 2 * this->Mass / rho) * State_vector[e_p_r];
    Derivatives[e_theta] = 1.0 / (rho * rho) * State_vector[e_p_theta];
    Derivatives[e_phi] = J / (rho * rho * sin2);
    Derivatives[e_p_phi] = 0.0;
    Derivatives[e_p_theta] = cos1 / (rho * rho * sin1 * sin2) * J * J;
    Derivatives[e_p_t] = 0.0;

    double r_term_1 = -this->Mass * r / (rho * rho * rho) * (1.0 / ((1 - 2 * this->Mass / rho) * (1 - 2 * this->Mass / rho)) + State_vector[e_p_r] * State_vector[e_p_r]);
    double r_term_2 = r / (rho * rho * rho * rho) * (State_vector[e_p_theta] * State_vector[e_p_theta] + J * J / sin2);

    Derivatives[e_p_r] = r_term_1 + r_term_2;

}

bool RBH_class::terminate_integration(const double* const State_vector) {

    const bool scatter = State_vector[e_r] > this->Scattering_radius && State_vector[e_p_r] < 0;

    const bool hit_horizon = State_vector[e_r] - this->Horizon_radius < this->Min_distance_to_singular_point;

    return scatter || hit_horizon;

}

void RBH_class::Convert_global_to_local_coords(const double* const State_Vector_Global, const double* const Global_Vec_to_Convert, double* Local_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

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

void RBH_class::Convert_local_to_global_coords(const double* const State_Vector_Local, const double* const Local_Vec_to_Convert, double* Global_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

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