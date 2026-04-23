#include "Spacetimes.h"

JNW_class::JNW_class(const Metric_parameters_type* const p_Metric_Parameters) {

    if (isnan(p_Metric_Parameters->Scattering_radius) or isinf(p_Metric_Parameters->Scattering_radius) or p_Metric_Parameters->Scattering_radius < 0) {

        throw std::runtime_error(std::format("Invalid value for the scattering radius: {}", p_Metric_Parameters->Scattering_radius));

    }

    if (isnan(p_Metric_Parameters->Min_distance_to_singular_point) or isinf(p_Metric_Parameters->Min_distance_to_singular_point) or p_Metric_Parameters->Min_distance_to_singular_point < 0) {

        throw std::runtime_error(std::format("Invalid value for the distance to the singular point: {}", p_Metric_Parameters->Min_distance_to_singular_point));

    }

    if (isnan(p_Metric_Parameters->JNW_Gamma_Parameter) or isinf(p_Metric_Parameters->JNW_Gamma_Parameter) or p_Metric_Parameters->JNW_Gamma_Parameter < 0) {

        throw std::runtime_error(std::format("Invalid value for the gamma parameter: {}", p_Metric_Parameters->JNW_Gamma_Parameter));

    }

    this->Mass = p_Metric_Parameters->Mass;
    this->Gamma = p_Metric_Parameters->JNW_Gamma_Parameter;
    this->Scattering_radius = p_Metric_Parameters->Scattering_radius;
    this->Min_distance_to_singular_point = p_Metric_Parameters->Min_distance_to_singular_point;
    this->Horizon_radius = 0.0;

    if (this->Gamma < 1.0) {

        this->Horizon_radius = this->Mass + sqrt(this->Mass * this->Mass - this->Gamma * this->Gamma);

    }

    this->Scattered_off_singulariy = false;

}

double* JNW_class::get_ISCO() {

    static double r_ISCO[2]{};

    double r_singularity = 2 * this->Mass / this->Gamma;

    if (this->Gamma > 1.0 / 2) {

        r_ISCO[Inner] = 1.0 / this->Gamma * (3.0 * this->Gamma + 1.0 + sqrt(5 * this->Gamma * this->Gamma - 1));
        r_ISCO[Outer] = r_ISCO[Inner];


    }else if (this->Gamma > 1.0 / sqrt(5) and this->Gamma < 1.0 / 2) {

            r_ISCO[Inner] = 1.0 / this->Gamma * (3.0 * this->Gamma + 1.0 - sqrt(5 * this->Gamma * this->Gamma - 1));
            r_ISCO[Outer] = 1.0 / this->Gamma * (3.0 * this->Gamma + 1.0 + sqrt(5 * this->Gamma * this->Gamma - 1));

    }
    else {

        r_ISCO[Inner] = r_singularity;
        r_ISCO[Outer] = r_singularity;

    }

    return r_ISCO;

 };

double* JNW_class::get_Photon_Sphere() {

    static double photon_orbit[2]{};

    double r_singularity = 2 * this->Mass / this->Gamma;

    if (this->Gamma > 0.5) { // Weak naked singularity

        photon_orbit[Inner] = (2 * this->Gamma + 1) * r_singularity / 2;
       

    }
    else {

        photon_orbit[Inner] = r_singularity;

    }

    photon_orbit[Outer] = photon_orbit[Inner];

    return photon_orbit;

};


Metric_type JNW_class::get_local_metric(const double* const Local_State_Vector) const {

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double r_singularity = 2 * this->Mass / this->Gamma;

    Metric_type s_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_Metric.Metric[e_t][e_t]         = -pow(1 - r_singularity / r, this->Gamma);
    s_Metric.Metric[e_r][e_r]         = -1.0 / s_Metric.Metric[e_t][e_t];
    s_Metric.Metric[e_theta][e_theta] = pow(1 - r_singularity / r, 1 - this->Gamma) * r2;
    s_Metric.Metric[e_phi][e_phi]     = s_Metric.Metric[e_theta][e_theta] * sin_theta * sin_theta;

    s_Metric.Lapse_function = sqrt(-s_Metric.Metric[e_t][e_t]);
    s_Metric.Shift_function = 0.;

    return s_Metric;

}

Metric_type JNW_class::get_global_metric(const double* const Global_State_Vector) const {

    return this->get_local_metric(Global_State_Vector);

}

Metric_type JNW_class::get_dr_local_metric(const double* const Local_State_Vector) const {

    Metric_type s_Metric = this->get_local_metric(Local_State_Vector);

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double r_singularity = 2 * this->Mass / this->Gamma;

    Metric_type s_dr_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_dr_Metric.Metric[e_t][e_t]         = -this->Gamma * pow(1 - r_singularity / r, this->Gamma - 1) * r_singularity / r2;
    s_dr_Metric.Metric[e_r][e_r]         = 1.0 / (s_Metric.Metric[e_t][e_t] * s_Metric.Metric[e_t][e_t]) * s_dr_Metric.Metric[e_t][e_t];
    s_dr_Metric.Metric[e_theta][e_theta] = 2 * r * pow(1 - r_singularity / r, 1 - this->Gamma) + (1 - this->Gamma) * pow(1 - r_singularity / r, -this->Gamma) * r_singularity;
    s_dr_Metric.Metric[e_phi][e_phi]     = s_dr_Metric.Metric[e_theta][e_theta] * sin_theta * sin_theta;

    s_dr_Metric.Lapse_function = -1. / (2 * sqrt(-s_Metric.Metric[e_t][e_t])) * s_dr_Metric.Metric[e_t][e_t];

    return s_dr_Metric;

}

Metric_type JNW_class::get_dr_global_metric(const double* const Global_State_Vector) const {

    return this->get_dr_local_metric(Global_State_Vector);

}

Metric_type JNW_class::get_dtheta_local_metric(const double* const Local_State_Vector) const {

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    double r_singularity = 2 * this->Mass / this->Gamma;

    Metric_type s_dtheta_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_dtheta_Metric.Metric[e_phi][e_phi] = 2 * r * r * sin_theta * cos_theta * pow(1 - r_singularity / r, 1 - this->Gamma);

    return s_dtheta_Metric;
}

Metric_type JNW_class::get_dtheta_global_metric(const double* const Global_State_Vector) const {

    return this->get_dtheta_local_metric(Global_State_Vector);

}

Metric_type JNW_class::get_d2r_local_metric(const double* const Local_State_Vector) const {

    Metric_type s_Metric = this->get_local_metric(Local_State_Vector);
    Metric_type s_dr_Metric = this->get_dr_local_metric(Local_State_Vector);

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double r_singularity = 2 * this->Mass / this->Gamma;

    Metric_type s_d2r_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_d2r_Metric.Metric[e_t][e_t] = -this->Gamma * (this->Gamma - 1) * pow(1 - r_singularity / r, this->Gamma - 2) * r_singularity * r_singularity / r2 / r2
        + 2 * this->Gamma * pow(1 - r_singularity / r, this->Gamma - 1) * r_singularity / r2 / r;

    s_d2r_Metric.Metric[e_r][e_r] = 1.0 / (s_Metric.Metric[e_t][e_t] * s_Metric.Metric[e_t][e_t]) * s_d2r_Metric.Metric[e_t][e_t]
        - 2.0 / (s_Metric.Metric[e_t][e_t] * s_Metric.Metric[e_t][e_t] * s_Metric.Metric[e_t][e_t]) * s_dr_Metric.Metric[e_t][e_t] * s_dr_Metric.Metric[e_t][e_t];

    s_d2r_Metric.Metric[e_theta][e_theta] = 2 * pow(1 - r_singularity / r, 1 - this->Gamma) + 2 * (1 - this->Gamma) * pow(1 - r_singularity / r, -this->Gamma) * r_singularity / r
        - (1 - this->Gamma) * this->Gamma * pow(1 - r_singularity / r, -this->Gamma - 1) * r_singularity * r_singularity / r2;

    s_d2r_Metric.Metric[e_phi][e_phi] = s_d2r_Metric.Metric[e_theta][e_theta] * sin_theta * sin_theta;

    s_d2r_Metric.Lapse_function = -s_d2r_Metric.Metric[e_t][e_t];

    return s_d2r_Metric;

}

void JNW_class::get_EOM(const double* const State_vector, double* const Derivatives) {

    const double& r = State_vector[e_r];
    const double& J = State_vector[e_p_phi];

    double sin1 = sin(State_vector[e_theta]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(State_vector[e_theta]);

    double r_singularity = 2 * this->Mass / this->Gamma;

    double pow_gamma = pow(1 - r_singularity / r, this->Gamma);
    double pow_gamma_minus_1 = pow(1 - r_singularity / r, this->Gamma - 1);

    Derivatives[e_t] = - 1 / pow_gamma * State_vector[e_p_t];
    Derivatives[e_r] = pow_gamma * State_vector[e_p_r];
    Derivatives[e_theta] = pow_gamma_minus_1 / (r * r) * State_vector[e_p_theta];
    Derivatives[e_phi] = pow_gamma_minus_1 / (r * r * sin2) * J;
    Derivatives[e_p_phi] = 0.0;
    Derivatives[e_p_theta] = pow_gamma_minus_1 * cos1 / (r * r * sin1 * sin2) * J * J;
    Derivatives[e_p_t] = 0.0;

    double r_term_1 = -this->Gamma * r_singularity / 2 / r / r * pow_gamma_minus_1 * (1.0 / pow_gamma / pow_gamma
                    + State_vector[e_p_r] * State_vector[e_p_r]);
    double r_term_2 = 1.0 / r / r / r * pow_gamma_minus_1 * (1 - r_singularity / 2 / r * (this->Gamma - 1) / (1 - r_singularity / r))
                    * (State_vector[e_p_theta] * State_vector[e_p_theta] + J * J / sin2);

    Derivatives[e_p_r] = r_term_1 + r_term_2;

}

bool JNW_class::terminate_integration(const double* const State_vector) {

    const double r_singularity = 2.0 / this->Gamma;

    bool hit_horizon = false;

    if (r_singularity < this->Horizon_radius) {

        hit_horizon = State_vector[e_r] - this->Horizon_radius < this->Min_distance_to_singular_point;

    }

    const bool scatter = State_vector[e_r] > this->Scattering_radius and State_vector[e_p_r] < 0.0;

    if (scatter or hit_horizon) { this->Scattered_off_singulariy = false; }

    return scatter or hit_horizon;
};


void JNW_class::Convert_global_to_local_coords(const double* const, const double* const Global_Vec_to_Convert, double* Local_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

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

void JNW_class::Convert_local_to_global_coords(const double* const, const double* const Local_Vec_to_Convert, double* Global_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

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