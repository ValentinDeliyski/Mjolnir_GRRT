#include "Spacetimes.h"

Gauss_Bonnet_class::Gauss_Bonnet_class(const Metric_parameters_type* const p_Metric_Parameters) {

    if (isnan(p_Metric_Parameters->Scattering_radius) or isinf(p_Metric_Parameters->Scattering_radius) or p_Metric_Parameters->Scattering_radius < 0) {

        throw std::runtime_error(std::format("Invalid value for the scattering radius: {}", p_Metric_Parameters->Scattering_radius));

    }

    if (isnan(p_Metric_Parameters->Min_distance_to_singular_point) or isinf(p_Metric_Parameters->Min_distance_to_singular_point) or p_Metric_Parameters->Min_distance_to_singular_point < 0) {

        throw std::runtime_error(std::format("Invalid value for the distance to the singular point: {}", p_Metric_Parameters->Min_distance_to_singular_point));

    }

    if (isnan(p_Metric_Parameters->GB_Gamma_Parameter) or isinf(p_Metric_Parameters->GB_Gamma_Parameter) or p_Metric_Parameters->GB_Gamma_Parameter < 0) {

        throw std::runtime_error(std::format("Invalid value for the gamma parameter: {}", p_Metric_Parameters->GB_Gamma_Parameter));

    }

    this->Mass = p_Metric_Parameters->Mass;
    this->Gamma = p_Metric_Parameters->GB_Gamma_Parameter;
    this->Scattering_radius = p_Metric_Parameters->Scattering_radius;
    this->Min_distance_to_singular_point = p_Metric_Parameters->Min_distance_to_singular_point;
    this->Horizon_radius = 0.0;

    if (this->Gamma < 1.0) {

        this->Horizon_radius = this->Mass + sqrt(this->Mass * this->Mass - this->Gamma * this->Gamma);

    }

}

double* Gauss_Bonnet_class::get_ISCO() {

    /**************************************************************************
    |                                                                         |
    |   @ Description: Returns a pointer to the inner and outer ISCO radii.   |
    |     * The outer ISCO is the solution to the equation:                   |
    |       d2r_f(r) + 3 * dr_f(r) / r - 2 * dr_f(r)**2 / f(r) = 0            |
    |     * GAUSS_BONNET_GAMMA is in the range [0, 1.5]                       |
    |                                                                         |
    |   @ Inputs: None                                                        |
    |                                                                         |
    |   @ Ouput: Pointer to the ISCO radii                                    |
    |                                                                         |
    **************************************************************************/

    static double r_ISCO[2]{};

    double fit_coeffs[11] = { 5.99998915, -0.61042681, -0.11593137,  0.07275861, -0.46946788,
                              1.20693793, -1.99054947,  2.05041439, -1.29496979,  0.45787902, -0.07008574 };

    double Gamma2  = this->Gamma * this->Gamma;
    double Gamma4  = Gamma2 * Gamma2;
    double Gamma8  = Gamma4 * Gamma4;
    double Gamma10 = Gamma8 * Gamma2;

    r_ISCO[Outer] = fit_coeffs[0]  + 
                    fit_coeffs[1]  * this->Gamma +
                    fit_coeffs[2]  * Gamma2 +
                    fit_coeffs[3]  * Gamma2 * this->Gamma +
                    fit_coeffs[4]  * Gamma4 +
                    fit_coeffs[5]  * Gamma4 * this->Gamma +
                    fit_coeffs[6]  * Gamma4 * Gamma2 + 
                    fit_coeffs[7]  * Gamma8 / this->Gamma +
                    fit_coeffs[8]  * Gamma8 + 
                    fit_coeffs[9]  * Gamma8 * this->Gamma +
                    fit_coeffs[10] * Gamma10;

    r_ISCO[Inner] = pow(this->Gamma, 1.0 / 3);

    return r_ISCO;

};

double* Gauss_Bonnet_class::get_Photon_Sphere() {

    /* This expression is the root of a cubic equation */

    double q =  8 * this->Mass * this->Gamma;
    double p = -9 * this->Mass * this->Mass;

    static double photon_orbits[2]{};

    photon_orbits[Outer] = 2 * sqrt(-p / 3) * cos(1. / 3 * acos(3. / 2 * q / p * sqrt(-3. / p)));
    photon_orbits[Inner] = 2 * sqrt(-p / 3) * cos(1. / 3 * acos(3. / 2 * q / p * sqrt(-3. / p)) + 2. * M_PI / 3);

    return photon_orbits;

};

Metric_type Gauss_Bonnet_class::get_local_metric(const double* const Local_State_Vector) const {

    const double& M = this->Mass;
    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double f = 1. + r2 / this->Gamma / 2. * (1. - sqrt(1. + 8. * this->Gamma * M / r2 / r));

    Metric_type s_Metric{};

    /* --- Only the non-zero components are explicitly evaluated. --- */

    s_Metric.Metric[e_t][e_t]         = -f;
    s_Metric.Metric[e_r][e_r]         = 1. / f;
    s_Metric.Metric[e_theta][e_theta] = r2;
    s_Metric.Metric[e_phi][e_phi]     = r2 * sin_theta * sin_theta;

    s_Metric.Lapse_function = sqrt(-s_Metric.Metric[e_t][e_t]);

    return s_Metric;

}

Metric_type Gauss_Bonnet_class::get_global_metric(const double* const Global_State_Vector) const {

    return this->get_local_metric(Global_State_Vector);

}

Metric_type Gauss_Bonnet_class::get_dr_local_metric(const double* const Local_State_Vector) const {

    const double& M = this->Mass;
    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double f = 1. + r2 / this->Gamma / 2. * (1. - sqrt(1. + 8. * this->Gamma * M / r2 / r));
    double dr_f = 2. / r * (f - 1.) + 6. * M / sqrt(r2 * r2 + 8. * this->Gamma * M * r);

    Metric_type s_dr_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_dr_Metric.Metric[e_t][e_t]         = -dr_f;
    s_dr_Metric.Metric[e_r][e_r]         = -1. / f / f * dr_f;
    s_dr_Metric.Metric[e_theta][e_theta] = 2. * r;
    s_dr_Metric.Metric[e_phi][e_phi]     = 2. * r * sin_theta * sin_theta;

    s_dr_Metric.Lapse_function = 0;

    return s_dr_Metric;

}

Metric_type Gauss_Bonnet_class::get_dr_global_metric(const double* const Global_State_Vector) const {

    return this->get_dr_local_metric(Global_State_Vector);

}

Metric_type Gauss_Bonnet_class::get_dtheta_local_metric(const double* const Local_State_Vector) const {

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    Metric_type s_dtheta_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_dtheta_Metric.Metric[e_phi][e_phi] = 2 * r * r * sin_theta * cos_theta;

    return s_dtheta_Metric;
}

Metric_type Gauss_Bonnet_class::get_dtheta_global_metric(const double* const Global_State_Vector) const {

    return this->get_dtheta_local_metric(Global_State_Vector);

}

Metric_type Gauss_Bonnet_class::get_d2r_local_metric(const double* const Local_State_Vector) const {

    const double& M = this->Mass;
    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double root = sqrt(r2 * r2 + 8 * this->Gamma * M * r);

    double f = 1 + r2 / this->Gamma / 2. * (1 - sqrt(1. + 8. * this->Gamma * M / r2 / r));
    double dr_f = 2. / r * (f - 1.) + 6 * M / root;
    double d2r_f = -2. / r2 * (f - 1.) + 2. / r * dr_f - 12. * M / root / root / root * (r2 * r + 2. * this->Gamma * M);

    Metric_type s_d2r_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_d2r_Metric.Metric[e_t][e_t]         = -d2r_f;
    s_d2r_Metric.Metric[e_r][e_r]         = 2. / f / f / f * dr_f - 1. / f / f * d2r_f;
    s_d2r_Metric.Metric[e_theta][e_theta] = 2.;
    s_d2r_Metric.Metric[e_phi][e_phi]     = 2. * sin_theta * sin_theta;

    s_d2r_Metric.Lapse_function = -s_d2r_Metric.Metric[e_t][e_t];

    return s_d2r_Metric;

}

void Gauss_Bonnet_class::get_EOM(const double* const State_vector, double* const Derivatives) {

    const double& r = State_vector[e_r];
    const double& J = State_vector[e_p_phi];

    double sin1 = sin(State_vector[e_theta]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(State_vector[e_theta]);

    double root = sqrt(1. + 8. * this->Gamma * this->Mass / r / r / r);

    double f    = 1. + r * r / this->Gamma / 2. * (1. - root);
    double dr_f = 2. / r * (f - 1.) + 6. * this->Mass / root / r / r;

    Derivatives[e_t] = - 1. / f * State_vector[e_p_t];
    Derivatives[e_r] = f * State_vector[e_p_r];
    Derivatives[e_theta] = 1. / (r * r) * State_vector[e_p_theta];
    Derivatives[e_phi] = J / (r * r * sin2);
    Derivatives[e_p_phi] = 0.0;
    Derivatives[e_p_theta] = cos1 / (r * r * sin1 * sin2) * J * J;
    Derivatives[e_p_t] = 0.0;

    double r_term_1 = -1. / 2 * (1.0 / f / f + State_vector[e_p_r] * State_vector[e_p_r]) * dr_f;
    double r_term_2 = 1.0 / r / r / r * (State_vector[e_p_theta] * State_vector[e_p_theta] + J * J / sin2);

    Derivatives[e_p_r] = r_term_1 + r_term_2;
}

bool Gauss_Bonnet_class::terminate_integration(const double* const State_vector) {

    const bool scatter = State_vector[e_r] > this->Scattering_radius and State_vector[e_p_r] < 0;

    const bool hit_horizon = State_vector[e_r] - this->Horizon_radius < this->Min_distance_to_singular_point;

    return scatter or hit_horizon;

};

void Gauss_Bonnet_class::Convert_global_to_local_coords(const double* const, const double* const Global_Vec_to_Convert, double* Local_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

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

void Gauss_Bonnet_class::Convert_local_to_global_coords(const double* const, const double* const Local_Vec_to_Convert, double* Global_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

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