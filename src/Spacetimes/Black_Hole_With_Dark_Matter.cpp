#include "Spacetimes.h"

Black_Hole_w_Dark_Matter_Halo_class::Black_Hole_w_Dark_Matter_Halo_class(const Metric_parameters_type* const p_Metric_Parameters) {

    if (isnan(p_Metric_Parameters->Compactness) or isinf(p_Metric_Parameters->Compactness) or p_Metric_Parameters->Compactness < 0) {

        throw std::runtime_error(std::format("Invalid value for the compactness: {}", p_Metric_Parameters->Compactness));

    }

    if (isnan(p_Metric_Parameters->Halo_Mass) or isinf(p_Metric_Parameters->Halo_Mass) or p_Metric_Parameters->Halo_Mass < 0) {

        throw std::runtime_error(std::format("Invalid value for the halo mass: {}", p_Metric_Parameters->Halo_Mass));

    }

    if (isnan(p_Metric_Parameters->Scattering_radius) or isinf(p_Metric_Parameters->Scattering_radius) or p_Metric_Parameters->Scattering_radius < 0) {

        throw std::runtime_error(std::format("Invalid value for the scattering radius: {}", p_Metric_Parameters->Scattering_radius));

    }

    if (isnan(p_Metric_Parameters->Min_distance_to_singular_point) or isinf(p_Metric_Parameters->Min_distance_to_singular_point) or p_Metric_Parameters->Min_distance_to_singular_point < 0) {

        throw std::runtime_error(std::format("Invalid value for the distance to the singular point: {}", p_Metric_Parameters->Min_distance_to_singular_point));

    }

    this->Mass = p_Metric_Parameters->Mass;
    this->Compactness = p_Metric_Parameters->Compactness;
    this->Halo_Mass = p_Metric_Parameters->Halo_Mass;
    this->Scattering_radius = p_Metric_Parameters->Scattering_radius;
    this->Min_distance_to_singular_point = p_Metric_Parameters->Min_distance_to_singular_point;

}

double* Black_Hole_w_Dark_Matter_Halo_class::get_ISCO() {

    /**************************************************************************
    |                                                                         |
    |   @ Description: Returns a pointer to the inner and outer ISCO radii.   |
    |     * The outer ISCO is the solution to the equation:                   |
    |       d2r_f(r) + 3 * dr_f(r) / r - 2 * dr_f(r)**2 / f(r) = 0            |
    |     * Compactness is in the range [0, 1]                                |
    |                                                                         |
    |   @ Inputs: None                                                        |
    |                                                                         |
    |   @ Ouput: Pointer to the ISCO radii                                    |
    |                                                                         |
    **************************************************************************/

    static double r_ISCO[2]{};
    double fit_coeffs[11]{};

    if (this->Halo_Mass > 1e2 + 1) {

        fit_coeffs[0]  =  6.00000012;
        fit_coeffs[1]  = -1.83930084e-05;
        fit_coeffs[2]  = -1.85031991e-02;
        fit_coeffs[3]  =  6.48180463e-02;
        fit_coeffs[4]  = -1.86565972e-01;
        fit_coeffs[5]  =  4.07869236e-01;
        fit_coeffs[6]  = -6.40616716e-01;
        fit_coeffs[7]  =  6.89893238e-01;
        fit_coeffs[8]  = -4.79853789e-01;
        fit_coeffs[9]  =  1.93338491e-01;
        fit_coeffs[10] = -3.41859843e-02;

    }
    else {

        fit_coeffs[0]  =  6.00001128;
        fit_coeffs[1]  = -1.66482364e-03;
        fit_coeffs[2]  = -1.85739425e+00;
        fit_coeffs[3]  =  6.95943604e+00;
        fit_coeffs[4]  = -1.95778349e+01;
        fit_coeffs[5]  =  4.17612581e+01;
        fit_coeffs[6]  = -6.44415347e+01;
        fit_coeffs[7]  =  6.85608949e+01;
        fit_coeffs[8]  = -4.72854005e+01;
        fit_coeffs[9]  =  1.89361995e+01;
        fit_coeffs[10] = -3.33318571;
    }

    double Compactness2  = this->Compactness  * this->Compactness;
    double Compactness4  = Compactness2 * Compactness2;
    double Compactness8  = Compactness4 * Compactness4;
    double Compactness10 = Compactness8 * Compactness2;

    r_ISCO[Outer] = fit_coeffs[0]  +
                    fit_coeffs[1]  * this->Compactness +
                    fit_coeffs[2]  * Compactness2 +
                    fit_coeffs[3]  * Compactness2 * this->Compactness +
                    fit_coeffs[4]  * Compactness4 +
                    fit_coeffs[5]  * Compactness4 * this->Compactness +
                    fit_coeffs[6]  * Compactness4 * Compactness2 +
                    fit_coeffs[7]  * Compactness8 / this->Compactness +
                    fit_coeffs[8]  * Compactness8 +
                    fit_coeffs[9]  * Compactness8 * this->Compactness +
                    fit_coeffs[10] * Compactness10;

    r_ISCO[Inner] = r_ISCO[Outer];

    return r_ISCO;

};

Metric_type Black_Hole_w_Dark_Matter_Halo_class::get_local_metric(const double* const Local_State_Vector) const {

    const double& M = this->Mass;
    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double A_0 = this->Halo_Mass / this->Compactness;

    double ksi = 2 * A_0 - this->Halo_Mass + 4 * M;
    double Y = sqrt(this->Halo_Mass / ksi) * (2 * atan((r + A_0 + this->Halo_Mass) / sqrt(this->Halo_Mass * ksi)) - M_PI);
    double f = (1 - 2 * M / r) * exp(Y);
    double m = M + this->Halo_Mass * r2 / (A_0 + r) / (A_0 + r) * (1 - 2 * M / r) * (1 - 2 * M / r);

    Metric_type s_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_Metric.Metric[e_t][e_t]         = -f;
    s_Metric.Metric[e_r][e_r]         = 1. / (1 - 2 * m / r);
    s_Metric.Metric[e_theta][e_theta] = r2;
    s_Metric.Metric[e_phi][e_phi]     = r2 * sin_theta * sin_theta;

    s_Metric.Lapse_function = sqrt(- s_Metric.Metric[e_t][e_t]);

    return s_Metric;

}

Metric_type Black_Hole_w_Dark_Matter_Halo_class::get_global_metric(const double* const Global_State_Vector) const {

    return this->get_local_metric(Global_State_Vector);

}

Metric_type Black_Hole_w_Dark_Matter_Halo_class::get_dr_local_metric(const double* const Local_State_Vector) const {

    const double& M = this->Mass;
    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double A_0 = this->Halo_Mass / this->Compactness;

    double ksi = 2 * A_0 - this->Halo_Mass + 4 * M;
    double Y = sqrt(this->Halo_Mass / ksi) * (2 * atan((r + A_0 + this->Halo_Mass) / sqrt(this->Halo_Mass * ksi)) - M_PI);
    double dr_Y = 2 / ksi / (1 + (r + A_0 + this->Halo_Mass) * (r + A_0 + this->Halo_Mass) / M / ksi);

    double exp_y = exp(Y);

    double f = (1 - 2 * M / r) * exp_y;
    double dr_f = 2 * M / r2 * exp_y + f * dr_Y;
    double m = M + this->Halo_Mass * r2 / (A_0 + r) / (A_0 + r) * (1 - 2 * M / r) * (1 - 2 * M / r);
    double dr_m = 2 * (1 - 2 * M / r) * ((1 - 2 * M / r) * (1 - r / (r + A_0)) * r + 2 * M) * this->Halo_Mass / (r + A_0) / (r + A_0);

    Metric_type s_dr_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_dr_Metric.Metric[e_t][e_t]         = -dr_f;
    s_dr_Metric.Metric[e_r][e_r]         = -1. / (1 - 2 * m / r) / (1 - 2 * m / r) * (2 * m / r2 - 2 / r * dr_m);
    s_dr_Metric.Metric[e_theta][e_theta] = r2;
    s_dr_Metric.Metric[e_phi][e_phi]     = r2 * sin_theta * sin_theta;

    s_dr_Metric.Lapse_function = 0;

    return s_dr_Metric;

}

Metric_type Black_Hole_w_Dark_Matter_Halo_class::get_dr_global_metric(const double* const Global_State_Vector) const {

    return this->get_dr_local_metric(Global_State_Vector);

}

Metric_type Black_Hole_w_Dark_Matter_Halo_class::get_dtheta_local_metric(const double* const Local_State_Vector) const {

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    Metric_type s_dtheta_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_dtheta_Metric.Metric[e_phi][e_phi] = 2 * r * r * sin_theta * cos_theta;

    return s_dtheta_Metric;

}

Metric_type Black_Hole_w_Dark_Matter_Halo_class::get_dtheta_global_metric(const double* const Global_State_Vector) const {

    return this->get_dtheta_local_metric(Global_State_Vector);

}

void Black_Hole_w_Dark_Matter_Halo_class::get_EOM(const double* const Global_State_Vector, double* const Derivatives) {

    const double& r = Global_State_Vector[e_r];

    const double& J = Global_State_Vector[e_p_phi];

    double sin1 = sin(Global_State_Vector[e_theta]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(Global_State_Vector[e_theta]);
    double cos2 = cos1 * cos1;

    const double& M = this->Mass;
    double r2 = r * r;

    double A_0 = this->Halo_Mass / this->Compactness;

    double ksi   = 2 * A_0 - this->Halo_Mass + 4 * M;
    double Y     = sqrt(this->Halo_Mass / ksi) * (2 * atan((r + A_0 + this->Halo_Mass) / sqrt(this->Halo_Mass * ksi)) - M_PI);
    double dr_Y = 2 / ksi / (1 + (r + A_0 + this->Halo_Mass) * (r + A_0 + this->Halo_Mass) / M / ksi);

    double exp_y = exp(Y);

    double f    = (1 - 2 * M / r) * exp_y;
    double dr_f = 2 * M / r2 * exp_y + f * dr_Y;
    double m    = M + this->Halo_Mass * r2 / (A_0 + r) / (A_0 + r) * (1 - 2 * M / r) * (1 - 2 * M / r);
    double dr_m = 2 * (1 - 2 * M / r) * ((1 - 2 * M / r) * (1 - r / (r + A_0)) * r + 2 * M) * this->Halo_Mass / (r + A_0) / (r + A_0);

    Derivatives[e_t] = -1 / f * Global_State_Vector[e_p_t];
    Derivatives[e_r] = (1 - 2 * m / r) * Global_State_Vector[e_p_r];
    Derivatives[e_theta] = 1. / (r * r) * Global_State_Vector[e_p_theta];
    Derivatives[e_phi] = J / (r * r * sin2);
    Derivatives[e_p_phi] = 0.0;
    Derivatives[e_p_theta] = cos1 / (r * r * sin1 * sin2) * J * J;
    Derivatives[e_p_t] = 0.0;

    double r_term_1 = -1. / 2 / f / f * dr_f + (dr_m / r - m / r2) * Global_State_Vector[e_p_r] * Global_State_Vector[e_p_r];
    double r_term_2 = 1.0 / r / r / r * (Global_State_Vector[e_p_theta] * Global_State_Vector[e_p_theta] + J * J / sin2);

    Derivatives[e_p_r] = r_term_1 + r_term_2;

}

bool Black_Hole_w_Dark_Matter_Halo_class::terminate_integration(const double* const State_vector) {

    bool scatter     = State_vector[e_r] > this->Scattering_radius and State_vector[e_p_r] < 0;
    bool hit_horizon = State_vector[e_r] - 2 * this->Mass < this->Min_distance_to_singular_point;

    return scatter or hit_horizon;

};

void Black_Hole_w_Dark_Matter_Halo_class::Convert_global_to_local_coords(const double* const State_Vector_Global, const double* const Global_Vec_to_Convert, double* Local_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

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

void Black_Hole_w_Dark_Matter_Halo_class::Convert_local_to_global_coords(const double* const State_Vector_Local, const double* const Local_Vec_to_Convert, double* Global_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

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