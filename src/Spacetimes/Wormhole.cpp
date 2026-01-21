#include "Spacetimes.h"

Wormhole_class::Wormhole_class(const Metric_parameters_type* const p_Metric_Parameters){


    if (isnan(p_Metric_Parameters->Spin) || isinf(p_Metric_Parameters->Spin)) {

        throw std::runtime_error(std::format("Invalid value for the spin parameter: {}", p_Metric_Parameters->Spin));

    }

    if (true != p_Metric_Parameters->Stop_At_Throat && false != p_Metric_Parameters->Stop_At_Throat) {

        throw std::runtime_error(std::format("Invalid value for the \"Stop at throat\" flag: {}", p_Metric_Parameters->Stop_At_Throat));

    }

    if (isnan(p_Metric_Parameters->Redshift_Parameter) || isinf(p_Metric_Parameters->Redshift_Parameter) || p_Metric_Parameters->Redshift_Parameter < 0) {

        throw std::runtime_error(std::format("Invalid value for the redshift parameter: {}", p_Metric_Parameters->Redshift_Parameter));

    }

    if (isnan(p_Metric_Parameters->Scattering_radius) || isinf(p_Metric_Parameters->Scattering_radius) || p_Metric_Parameters->Scattering_radius < 0) {

        throw std::runtime_error(std::format("Invalid value for the scattering radius: {}", p_Metric_Parameters->Scattering_radius));

    }

    if (isnan(p_Metric_Parameters->Min_distance_to_singular_point) || isinf(p_Metric_Parameters->Min_distance_to_singular_point) || p_Metric_Parameters->Min_distance_to_singular_point < 0) {

        throw std::runtime_error(std::format("Invalid value for the distance to the throat: {}", p_Metric_Parameters->Min_distance_to_singular_point));

    }

    this->Mass = p_Metric_Parameters->Mass;
    this->R_Throat = p_Metric_Parameters->R_throat;
    this->Spin_Param = p_Metric_Parameters->Spin;
    this->Redshift_Param = p_Metric_Parameters->Redshift_Parameter;
    this->Stop_at_Throat = p_Metric_Parameters->Stop_At_Throat;
    this->Scattering_radius = p_Metric_Parameters->Scattering_radius;

    this->Affine_param_at_throat_corssing = 0.0;
    this->Crossed_throat = false;

    // I just reuse the "Min_distance_to_singular_point" for the min throat distance because it serves the same purpose.
    this->Min_distance_to_throat = p_Metric_Parameters->Min_distance_to_singular_point;

}

double* Wormhole_class::get_ISCO() {

    double M = this->Mass;

    static double r_ISCO[2]{};

    if (this->Spin_Param < 0.016) {

        r_ISCO[Inner] = 2 * M * (sqrt(4. / 9 * (6 * this->Redshift_Param + 1)) * cosh(1. / 3 * acosh((1 + 9 * this->Redshift_Param + 27. / 2 * this->Redshift_Param * this->Redshift_Param) / pow(6 * this->Redshift_Param + 1, 3. / 2))) + 1. / 3);
        r_ISCO[Outer] = r_ISCO[Inner];
    }
    else {

        r_ISCO[Inner] = this->R_Throat;
        r_ISCO[Outer] = r_ISCO[Inner];

    }

    return r_ISCO;

}

double* Wormhole_class::get_Photon_Sphere() {

    double M = this->Mass;
    double a = this->Spin_Param;

    static double photon_orbit[2]{};

    photon_orbit[Inner] = M / 2 * (1 + sqrt(1 + 8 * this->Redshift_Param));
    photon_orbit[Outer] = photon_orbit[Inner];

    return photon_orbit;

}

Metric_type Wormhole_class::get_local_metric(const double* const Local_State_Vector) const {

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double exponent = -this->Mass / r - this->Redshift_Param * this->Mass * this->Mass / r2;

    Metric_type s_Metric{};

    /* --- Only the non-zero components are explicitly evaluated. --- */

    s_Metric.Lapse_function = exp(exponent);
    s_Metric.Shift_function = 2 * this->Spin_Param * this->Mass * this->Mass / r2 / r;

    s_Metric.Metric[e_t][e_t] = -s_Metric.Lapse_function * s_Metric.Lapse_function +
                                r2 * s_Metric.Shift_function * s_Metric.Shift_function * sin_theta * sin_theta;

    s_Metric.Metric[e_t][e_phi]       = -r2 * sin_theta * sin_theta * s_Metric.Shift_function;
    s_Metric.Metric[e_phi][e_t]       = s_Metric.Metric[e_t][e_phi];

    s_Metric.Metric[e_r][e_r]         = 1 / (1 - this->R_Throat / r);
    s_Metric.Metric[e_theta][e_theta] = r2;
    s_Metric.Metric[e_phi][e_phi]     = r2 * sin_theta * sin_theta;

    return s_Metric;
}

Metric_type Wormhole_class::get_global_metric(const double* const Global_State_Vector) const {

    const double& Ell = Global_State_Vector[e_r];
    const double& theta = Global_State_Vector[e_theta];

    double r = sqrt(Ell * Ell + this->R_Throat * this->R_Throat);
    double r2 = r * r;
    double sin_theta = sin(theta);

    double exponent = -this->Mass / r - this->Redshift_Param * this->Mass * this->Mass / r2;

    Metric_type s_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_Metric.Lapse_function = exp(exponent);
    s_Metric.Shift_function = 2 * this->Spin_Param * this->Mass * this->Mass / r2 / r;

    s_Metric.Metric[e_t][e_t] = -s_Metric.Lapse_function * s_Metric.Lapse_function +
        r2 * s_Metric.Shift_function * s_Metric.Shift_function * sin_theta * sin_theta;

    s_Metric.Metric[e_t][e_phi] = -r2 * sin_theta * sin_theta * s_Metric.Shift_function;
    s_Metric.Metric[e_phi][e_t] = s_Metric.Metric[e_t][e_phi];

    s_Metric.Metric[e_r][e_r] = (1 + this->R_Throat / r);
    s_Metric.Metric[e_theta][e_theta] = r2;
    s_Metric.Metric[e_phi][e_phi] = r2 * sin_theta * sin_theta;

    return s_Metric;
}

Metric_type Wormhole_class::get_dr_local_metric(const double* const Local_State_Vector) const {

    Metric_type s_Metric = this->get_local_metric(Local_State_Vector);

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    Metric_type s_dr_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_dr_Metric.Lapse_function = s_Metric.Lapse_function * (1 / r2 + 2 * this->Redshift_Param / (r2 * r));
    s_dr_Metric.Shift_function = -3 * s_Metric.Shift_function / r;

    double& N = s_Metric.Lapse_function;
    double& dr_N = s_dr_Metric.Lapse_function;
    double& omega = s_Metric.Shift_function;
    double& dr_omega = s_dr_Metric.Shift_function;

    s_dr_Metric.Metric[e_t][e_t]         = -2 * N * dr_N + 2 * r * omega * (omega + r * dr_omega) * sin_theta * sin_theta;
    s_dr_Metric.Metric[e_t][e_phi]       = -r * (2 * omega + r * dr_omega) * sin_theta * sin_theta;
    s_dr_Metric.Metric[e_phi][e_t]       = s_dr_Metric.Metric[e_t][e_phi];
    s_dr_Metric.Metric[e_r][e_r]         = -1. / ((1 - this->R_Throat / r) * (1 - this->R_Throat / r)) * (this->R_Throat / r2);
    s_dr_Metric.Metric[e_theta][e_theta] = 2 * r;
    s_dr_Metric.Metric[e_phi][e_phi]     = 2 * r * sin_theta * sin_theta;

    return s_dr_Metric;
}

Metric_type Wormhole_class::get_dr_global_metric(const double* const Global_State_Vector) const {

    double Local_State_Vector[e_Full_state_size]{};
    memcpy(Local_State_Vector, Global_State_Vector, e_Full_state_size * sizeof(double));
    Local_State_Vector[e_r] = sqrt(Global_State_Vector[e_r] * Global_State_Vector[e_r] + this->R_Throat * this->R_Throat);

    /* ------------ This is called to get the lapse and shift functions. Technically its not nessicery to call the local metric here,
                    as the lapse and shift give the same value in both coordinates (no g_tr term). ------------ */
    Metric_type s_Metric = this->get_local_metric(Local_State_Vector);

    const double& Ell = Global_State_Vector[e_r];
    const double& theta = Global_State_Vector[e_theta];

    double& r = Local_State_Vector[e_r];
    double r2 = r * r;
    double sin_theta = sin(theta);

    double dr_dEll = Ell / r;

    Metric_type s_dEll_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    /* -------- The derivatives wrt Ell have an extra dr_dEll factor infront, which will be added on later. For now it is convenient to 
                use these as derivatives wrt r. ------------ */
    s_dEll_Metric.Lapse_function = s_Metric.Lapse_function * (1 / r2 + 2 * this->Redshift_Param / (r2 * r));
    s_dEll_Metric.Shift_function = -3 * s_Metric.Shift_function / r;

    double& N = s_Metric.Lapse_function;
    double& dr_N = s_dEll_Metric.Lapse_function;
    double& omega = s_Metric.Shift_function;
    double& dr_omega = s_dEll_Metric.Shift_function;

    s_dEll_Metric.Metric[e_t][e_t] = (-2 * N * dr_N + 2 * r * omega * (omega + r * dr_omega) * sin_theta * sin_theta) * dr_dEll;
    s_dEll_Metric.Metric[e_t][e_phi] = (-r * (2 * omega + r * dr_omega) * sin_theta * sin_theta)* dr_dEll;
    s_dEll_Metric.Metric[e_phi][e_t] = s_dEll_Metric.Metric[e_t][e_phi];

    s_dEll_Metric.Metric[e_r][e_r] = -this->R_Throat / r2 * dr_dEll;

    s_dEll_Metric.Metric[e_theta][e_theta] = 2 * r * dr_dEll;
    s_dEll_Metric.Metric[e_phi][e_phi] = 2 * r * sin_theta * sin_theta * dr_dEll;

    s_dEll_Metric.Lapse_function *= dr_dEll;
    s_dEll_Metric.Shift_function *= dr_dEll;

    return s_dEll_Metric;
}

Metric_type Wormhole_class::get_dtheta_local_metric(const double* const Local_State_Vector) const {

    Metric_type s_Metric = this->get_local_metric(Local_State_Vector);

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    Metric_type s_dtheta_Metric{};

    /* --- Only the non-zero components are explicitly evaluated. --- */

    s_dtheta_Metric.Metric[e_t][e_t]     = 2 * r2 * s_Metric.Shift_function * s_Metric.Shift_function * sin_theta * cos_theta;
    s_dtheta_Metric.Metric[e_t][e_phi]   = -2 * r2 * sin_theta * cos_theta * s_Metric.Shift_function;
    s_dtheta_Metric.Metric[e_phi][e_t]   = s_dtheta_Metric.Metric[e_t][e_phi];
    s_dtheta_Metric.Metric[e_phi][e_phi] = 2 * r2 * sin_theta * cos_theta;

    return s_dtheta_Metric;
}


Metric_type Wormhole_class::get_dtheta_global_metric(const double* const Global_State_Vector) const {

    double Local_State_Vector[e_Full_state_size]{};
    memcpy(Local_State_Vector, Global_State_Vector, e_Full_state_size * sizeof(double));
    Local_State_Vector[e_r] = sqrt(Global_State_Vector[e_r] * Global_State_Vector[e_r] + this->R_Throat * this->R_Throat);

    Metric_type s_Metric = this->get_local_metric(Local_State_Vector);

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    Metric_type s_dtheta_Metric{};

    /* ------------------------------------ Only the non-zero components are explicitly evaluated. ------------------------------------ */

    s_dtheta_Metric.Metric[e_t][e_t] = 2 * r2 * s_Metric.Shift_function * s_Metric.Shift_function * sin_theta * cos_theta;
    s_dtheta_Metric.Metric[e_t][e_phi] = -2 * r2 * sin_theta * cos_theta * s_Metric.Shift_function;
    s_dtheta_Metric.Metric[e_phi][e_t] = s_dtheta_Metric.Metric[e_t][e_phi];
    s_dtheta_Metric.Metric[e_phi][e_phi] = 2 * r2 * sin_theta * cos_theta;

    return s_dtheta_Metric;
}

Metric_type Wormhole_class::get_d2r_local_metric(const double* const State_Vector) const {

    Metric_type s_Metric = this->get_local_metric(State_Vector);
    Metric_type s_dr_Metric = this->get_dr_local_metric(State_Vector);

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double& N = s_Metric.Lapse_function;
    double& dr_N = s_dr_Metric.Lapse_function;
    double& omega = s_Metric.Shift_function;
    double& dr_omega = s_dr_Metric.Shift_function;

    Metric_type s_d2r_Metric{};

    /* --- Only the non-zero components are explicitly evaluated. --- */

    s_d2r_Metric.Lapse_function = dr_N * (1 / r2 + 2 * this->Redshift_Param / (r2 * r)) - N * (2. / (r2 * r) + 6 * this->Redshift_Param / (r2 * r2));
    s_d2r_Metric.Shift_function = -3 * dr_omega / r + 3 * omega / r2;

    s_d2r_Metric.Metric[e_t][e_t] = -2 * dr_N * dr_N - 2 * N * s_d2r_Metric.Lapse_function + 2 * ((omega + r * dr_omega) * (omega + r * dr_omega) +
                                    r * omega * (dr_omega + dr_omega + r * s_d2r_Metric.Shift_function)) * sin_theta * sin_theta;
    s_d2r_Metric.Metric[e_t][e_phi]       = -(2 * omega + r * dr_omega + r * (3 * dr_omega + r * s_d2r_Metric.Shift_function)) * sin_theta * sin_theta;
    s_d2r_Metric.Metric[e_phi][e_t]       = s_d2r_Metric.Metric[e_t][e_phi];
    s_d2r_Metric.Metric[e_r][e_r]         = 2 / ((1 - this->R_Throat / r) * (1 - this->R_Throat / r)) * ((this->R_Throat / r2) * (this->R_Throat / r2) / (1 - this->R_Throat / r) + this->R_Throat / (r2 * r));
    s_d2r_Metric.Metric[e_theta][e_theta] = 2.0;
    s_d2r_Metric.Metric[e_phi][e_phi]     = 2 * sin_theta * sin_theta;

    return s_d2r_Metric;
}

void Wormhole_class::get_EOM(const double* const State_Vector, double* const Derivatives) {

    double sqrt_r2 = sqrt(State_Vector[e_r] * State_Vector[e_r] + this->R_Throat * this->R_Throat);
    double d_ell_r = State_Vector[e_r] / sqrt_r2;

    const double& J = State_Vector[e_p_phi];

    double omega = 2 * this->Spin_Param * this->Mass * this->Mass / sqrt_r2 / sqrt_r2 / sqrt_r2;
    double d_ell_omega = -3 * omega / sqrt_r2 * d_ell_r;

    double exponent = -1 / sqrt_r2 - this->Redshift_Param / (sqrt_r2 * sqrt_r2);
    double N = exp(-this->Mass / sqrt_r2 - this->Redshift_Param * this->Mass * this->Mass / sqrt_r2 / sqrt_r2);
    double d_ell_N = N * (1 / (sqrt_r2 * sqrt_r2) + 2 * this->Redshift_Param / (sqrt_r2 * sqrt_r2 * sqrt_r2)) * d_ell_r;

    double N2 = N * N;

    double sin1 = sin(State_Vector[e_theta]);
    double sin2 = sin1 * sin1;

    Derivatives[e_t] = -1.0 / N / N * State_Vector[e_p_t] + omega * omega / N / N * State_Vector[e_p_phi];
    Derivatives[e_r] = 1.0 / (1 + this->R_Throat / sqrt_r2) * State_Vector[e_p_r];
    Derivatives[e_theta] = 1.0 / (sqrt_r2 * sqrt_r2) * State_Vector[e_p_theta];
    Derivatives[e_phi] = J / (sqrt_r2 * sqrt_r2 * sin2) + omega * (1 - omega * J) / N2;
    Derivatives[e_p_phi] = 0.0;
    Derivatives[e_p_theta] = (cos(State_Vector[e_theta]) / sin1) / (sqrt_r2 * sqrt_r2) * J * J / sin2;
    Derivatives[e_p_t] = 0.0;

    double term_1 = -1.0 / ((1 + this->R_Throat / sqrt_r2) * (1 + this->R_Throat / sqrt_r2)) * this->R_Throat * State_Vector[e_r] / (sqrt_r2 * sqrt_r2 * sqrt_r2) * State_Vector[e_p_r] * State_Vector[e_p_r] / 2;
    double term_2 = 1.0 / (sqrt_r2 * sqrt_r2 * sqrt_r2) * (State_Vector[e_p_theta] * State_Vector[e_p_theta] + J * J / sin2) * d_ell_r;
    double term_3 = -(1.0 / (N2 * N) * d_ell_N * ((1 - omega * J) * (1 - omega * J)) - 1.0 / N2 * (-d_ell_omega * (1 - omega * J) * J));

    Derivatives[e_p_r] = term_1 + term_2 + term_3;

}

bool Wormhole_class::terminate_integration(const double* const State_vector) {

    const double Scatter_radius_global_coords = sqrt(this->Scattering_radius * this->Scattering_radius - this->R_Throat * this->R_Throat);

    const bool scatter            = State_vector[e_r] > Scatter_radius_global_coords && State_vector[e_p_r] < 0;
    const bool scatter_other_side = State_vector[e_r] < -Scatter_radius_global_coords;
    const bool stop_at_throat     = State_vector[e_r] < this->Min_distance_to_throat;

    if (this->Stop_at_Throat) {

        return scatter || stop_at_throat;
    }
    else {

        return scatter || scatter_other_side;

    }
}

void Wormhole_class::Convert_global_to_local_coords(const double* const State_Vector_Global, const double* const Global_Vec_to_Convert, double* Local_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

    switch (Entry_to_convert) {

    case e_Full_State_Vector:

        if (State_Vector_Global[e_r] < 0) {

            int test{};

        }

        memcpy(Local_Vec_to_Convert, Global_Vec_to_Convert, e_Full_state_size * sizeof(double));
        Local_Vec_to_Convert[e_r] = sqrt(Global_Vec_to_Convert[e_r] * Global_Vec_to_Convert[e_r] + this->R_Throat * this->R_Throat);
        Local_Vec_to_Convert[e_p_r] *= Local_Vec_to_Convert[e_r] / Global_Vec_to_Convert[e_r];

        break;

    case e_Contravariant_vector:

        memcpy(Local_Vec_to_Convert, Global_Vec_to_Convert, 4 * sizeof(double));
        Local_Vec_to_Convert[e_r] *= State_Vector_Global[e_r] / sqrt(State_Vector_Global[e_r] * State_Vector_Global[e_r] + this->R_Throat * this->R_Throat);

        break;

    case e_Covariant_vector:

        memcpy(Local_Vec_to_Convert, Global_Vec_to_Convert, 4 * sizeof(double));
        Local_Vec_to_Convert[e_r] *= sqrt(State_Vector_Global[e_r] * State_Vector_Global[e_r] + this->R_Throat * this->R_Throat) / State_Vector_Global[e_r];

        break;

    case e_Coordinates:

        memcpy(Local_Vec_to_Convert, Global_Vec_to_Convert, 4 * sizeof(double));
        Local_Vec_to_Convert[e_r] = sqrt(Global_Vec_to_Convert[e_r] * Global_Vec_to_Convert[e_r] - this->R_Throat * this->R_Throat);

        break;

    default:

        throw std::runtime_error("Unsupported coordinate conversion type. Something Broke in Convert_global_to_local_coords()!");

    }

}

void Wormhole_class::Convert_local_to_global_coords(const double* const State_Vector_Global, const double* const Local_Vec_to_Convert, double* Global_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

    switch (Entry_to_convert) {

    case e_Contravariant_vector:

        memcpy(Global_Vec_to_Convert, Local_Vec_to_Convert, 4 * sizeof(double));
        Global_Vec_to_Convert[e_r] /= State_Vector_Global[e_r] / sqrt(State_Vector_Global[e_r] * State_Vector_Global[e_r] + this->R_Throat * this->R_Throat);

        break;

    case e_Covariant_vector:

        memcpy(Global_Vec_to_Convert, Local_Vec_to_Convert, 4 * sizeof(double));
        Global_Vec_to_Convert[e_r] /= sqrt(State_Vector_Global[e_r] * State_Vector_Global[e_r] + this->R_Throat * this->R_Throat) / State_Vector_Global[e_r];

        break;

    default:

        throw std::runtime_error("Unsupported coordinate conversion type. Something Broke in Convert_global_to_local_coords()!");

    }

}