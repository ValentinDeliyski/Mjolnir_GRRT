#include "Spacetimes.h"

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

Metric_type Wormhole_class::get_metric(const double* const State_Vector) const {

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double exponent = -this->Mass / r - this->Redshift_Param * this->Mass * this->Mass / r2;

    Metric_type s_Metric{};

    /* --- Only the non-zero components are exlicitly evaluated. --- */

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

Metric_type Wormhole_class::get_dr_metric(const double* const State_Vector) const {

    Metric_type s_Metric = this->get_metric(State_Vector);

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    Metric_type s_dr_Metric{};

    /* --- Only the non-zero components are exlicitly evaluated. --- */

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

Metric_type Wormhole_class::get_dtheta_metric(const double* const State_Vector) const {

    Metric_type s_Metric = this->get_metric(State_Vector);

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    double exponent = -this->Mass / r - this->Redshift_Param * this->Mass * this->Mass / r2;

    Metric_type s_dtheta_Metric{};

    /* --- Only the non-zero components are exlicitly evaluated. --- */

    s_dtheta_Metric.Metric[e_t][e_t]     = 2 * r2 * s_Metric.Shift_function * s_Metric.Shift_function * sin_theta * cos_theta;
    s_dtheta_Metric.Metric[e_t][e_phi]   = -2 * r2 * sin_theta * cos_theta * s_Metric.Shift_function;
    s_dtheta_Metric.Metric[e_phi][e_t]   = s_dtheta_Metric.Metric[e_t][e_phi];
    s_dtheta_Metric.Metric[e_phi][e_phi] = 2 * r2 * sin_theta * cos_theta;

    return s_dtheta_Metric;
}

Metric_type Wormhole_class::get_d2r_metric(const double* const State_Vector) const {

    Metric_type s_Metric = this->get_metric(State_Vector);
    Metric_type s_dr_Metric = this->get_dr_metric(State_Vector);

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double& N = s_Metric.Lapse_function;
    double& dr_N = s_dr_Metric.Lapse_function;
    double& omega = s_Metric.Shift_function;
    double& dr_omega = s_dr_Metric.Shift_function;

    Metric_type s_d2r_Metric{};

    /* --- Only the non-zero components are exlicitly evaluated. --- */

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

int Wormhole_class::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {


    double& r_obs = p_Initial_Conditions->Observer_params.distance;
    double& theta_obs = p_Initial_Conditions->Observer_params.inclination;

    p_Initial_Conditions->Init_Momentum[e_phi] = -J_data[photon] * sin(theta_obs);
    p_Initial_Conditions->Init_Momentum[e_theta] = p_theta_data[photon];

    double& p_theta = p_Initial_Conditions->Init_Momentum[e_theta];
    double& J = p_Initial_Conditions->Init_Momentum[e_phi];

    double& N = p_Initial_Conditions->Init_metric.Lapse_function;
    double& omega = p_Initial_Conditions->Init_metric.Shift_function;

    double rad_potential = -(p_theta * p_theta + J * J / sin(theta_obs) / sin(theta_obs)) * N * N / r_obs / r_obs + (1 - omega * J) * (1 - omega * J);

    double(*metric)[4] = p_Initial_Conditions->Init_metric.Metric;

    p_Initial_Conditions->Init_Momentum[e_r] = sqrt(rad_potential) / N / sqrt(metric[1][1]) * r_obs / sqrt(pow(r_obs, 2) - pow(1, 2)) * metric[1][1];

    return 0;
}

void Wormhole_class::get_EOM(const double* const State_Vector, double* const Derivatives) const{

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

    Derivatives[e_t] = 1.0 / N / N * State_Vector[e_p_r];
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

    const double Scatter_radius_global_coords = sqrt(this->Scattering_radius * this->Scattering_radius + this->R_Throat * this->R_Throat);

    const bool scatter            = State_vector[e_r] > Scatter_radius_global_coords && State_vector[e_p_r] < 0;
    const bool scatter_other_side = State_vector[e_r] < -Scatter_radius_global_coords;
    const bool stop_at_throat     = State_vector[e_r] < this->Min_distance_to_throat;

    if (this->Stop_at_Throat) {

        return scatter || stop_at_throat;
    }
    else {

        return scatter || scatter_other_side;

    }
};

Return_Values Wormhole_class::load_parameters(const Metric_parameters_type* const p_Metric_Parameters) {

    if (isnan(p_Metric_Parameters->Spin) || isinf(p_Metric_Parameters->Spin)) {

        std::cout << "Invalid value for the spin parameter: " << p_Metric_Parameters->Spin << "\n";

        return ERROR;
    }

    if (isnan(p_Metric_Parameters->Stop_At_Throat) || isinf(p_Metric_Parameters->Stop_At_Throat)) {

        std::cout << "Invalid value for the \"Stop at throat\" flag: " << p_Metric_Parameters->Stop_At_Throat << "\n";

        return ERROR;
    }

    if (isnan(p_Metric_Parameters->Redshift_Parameter) || isinf(p_Metric_Parameters->Redshift_Parameter) || p_Metric_Parameters->Redshift_Parameter < 0) {

        std::cout << "Invalid value for the redshift parameter: " << p_Metric_Parameters->Redshift_Parameter << "\n";

        return ERROR;
    }

    if (isnan(p_Metric_Parameters->Scattering_radius) || isinf(p_Metric_Parameters->Scattering_radius) || p_Metric_Parameters->Scattering_radius < 0) {

        std::cout << "Invalid value for the scattering radius: " << p_Metric_Parameters->Scattering_radius << "\n";

        return ERROR;
    }

    if (isnan(p_Metric_Parameters->Min_distance_to_singular_point) || isinf(p_Metric_Parameters->Min_distance_to_singular_point) || p_Metric_Parameters->Min_distance_to_singular_point < 0) {

        std::cout << "Invalid value for the distance to the throat: " << p_Metric_Parameters->Min_distance_to_singular_point << "\n";

        return ERROR;
    }
 
    this->Spin_Param = p_Metric_Parameters->Spin;
    this->Redshift_Param = p_Metric_Parameters->Redshift_Parameter;
    this->Stop_at_Throat = p_Metric_Parameters->Stop_At_Throat;
    this->Scattering_radius = p_Metric_Parameters->Scattering_radius;

    // I just reuse the "Min_distance_to_singular_point" for the min throat distance because it serves the same purpose.
    this->Min_distance_to_throat = p_Metric_Parameters->Min_distance_to_singular_point;

    return OK;


}
