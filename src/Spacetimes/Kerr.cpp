#include "Spacetimes.h"

double* Kerr_class::get_ISCO() {

    double Z_1 = 1 + pow(1 - this->Spin_Param * this->Spin_Param, 1. / 3) * (pow(1 + this->Spin_Param, 1. / 3) + pow(1 - this->Spin_Param, 1. / 3));
    double Z_2 = sqrt(3 * this->Spin_Param * this->Spin_Param + Z_1 * Z_1);

    static double r_ISCO[2]{};

    r_ISCO[Inner] = this->Mass * (3 + Z_2 - sqrt((3 - Z_1) * (3 + Z_1 + 2 * Z_2)));
    r_ISCO[Outer] = this->Mass * (3 + Z_2 + sqrt((3 - Z_1) * (3 + Z_1 + 2 * Z_2)));

    return r_ISCO;

}

double* Kerr_class::get_Photon_Sphere() {

    static double photon_orbit[2]{};

    photon_orbit[Inner] = 2 * this->Mass * (1 + cos(2.0 / 3 * acos(this->Spin_Param)));
    photon_orbit[Outer] = 2 * this->Mass * (1 + cos(2.0 / 3 * acos(-this->Spin_Param)));

    return photon_orbit;

}

Metric_type Kerr_class::get_metric(const double* const State_Vector) const {

    const double& M = this->Mass;
    const double& a = this->Spin_Param;

    const double& r     = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r2 + a * a * cos_theta * cos_theta;
    double delta = r2 - 2 * M * r + a * a;

    Metric_type s_Metric{};

    /* --- Only the non-zero components are exlicitly evaluated. --- */

    s_Metric.Metric[e_t][e_t]         = -(1 - 2 * M * r / rho2);
    s_Metric.Metric[e_t][e_phi]       = -2 * M * r * a * sin_theta * sin_theta / rho2;
    s_Metric.Metric[e_phi][e_t]       = s_Metric.Metric[e_t][e_phi];
    s_Metric.Metric[e_r][e_r]         = rho2 / delta;
    s_Metric.Metric[e_theta][e_theta] = rho2;
    s_Metric.Metric[e_phi][e_phi]     = (r2 + a * a + 2 * M * r * a * a / rho2 * sin_theta * sin_theta) * sin_theta * sin_theta;

    double sigma2 = (r2 + a * a) * (r2 + a * a) - a * a * delta * sin_theta * sin_theta;

    s_Metric.Lapse_function = sqrt(-s_Metric.Metric[e_t][e_t] + s_Metric.Metric[e_t][e_phi] * s_Metric.Metric[e_t][e_phi] / s_Metric.Metric[e_phi][e_phi]);
    s_Metric.Shift_function = 2 * a * r / sigma2;

    return s_Metric;

};

Metric_type Kerr_class::get_dr_metric(const double* const State_Vector) const {

    Metric_type s_Metric = this->get_metric(State_Vector);

    const double& M = this->Mass;
    const double& a = this->Spin_Param;

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r2 + a * a * cos_theta * cos_theta;
    double delta = r2 - 2 * M * r + a * a;

    Metric_type s_dr_Metric{};

    /* --- Only the non-zero components are exlicitly evaluated. --- */

    s_dr_Metric.Metric[e_t][e_t]         = -2 * M / rho2 * (2 * r2 / rho2 - 1);
    s_dr_Metric.Metric[e_t][e_phi]       = 2 * M * a * sin_theta * sin_theta / rho2 * (2 * r2 / rho2 - 1);
    s_dr_Metric.Metric[e_phi][e_t]       = s_dr_Metric.Metric[e_t][e_phi];
    s_dr_Metric.Metric[e_r][e_r]         = 2 * r / delta * (1 - rho2 / delta * (1 - M / r));
    s_dr_Metric.Metric[e_theta][e_theta] = 2 * r;
    s_dr_Metric.Metric[e_phi][e_phi]     = 2 * (r - M * a * a / rho2 * (2 * r2 / rho2 - 1) * sin_theta * sin_theta) * sin_theta * sin_theta;

    double sigma2 = rho2 * s_Metric.Metric[e_phi][e_phi] / sin_theta / sin_theta;
    double dr_sigma2 = 2 * r * s_Metric.Metric[e_phi][e_phi] + rho2 * s_dr_Metric.Metric[e_phi][e_phi];

    s_dr_Metric.Lapse_function = s_Metric.Lapse_function * (r / rho2 + (r - M) / delta - dr_sigma2 / 2 / sigma2);
    s_dr_Metric.Shift_function = s_Metric.Shift_function / r * (1 - r * dr_sigma2 / sigma2);

    return s_dr_Metric;
}

Metric_type Kerr_class::get_dtheta_metric(const double* const State_Vector) const {

    Metric_type s_Metric = this->get_metric(State_Vector);

    const double& M = this->Mass;
    const double& a = this->Spin_Param;

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r2 + a * a * cos_theta * cos_theta;
    double delta = r2 - 2 * M * r + a * a;

    Metric_type s_dtheta_Metric{};

    /* --- Only the non-zero components are exlicitly evaluated. --- */

    s_dtheta_Metric.Metric[e_t][e_t]         = 4 * M * r / rho2 / rho2 * (a * a * cos_theta * sin_theta);
    s_dtheta_Metric.Metric[e_t][e_phi]       = -4 * M * r * a * sin_theta * cos_theta / rho2 * (1 - a * a * sin_theta * sin_theta / rho2);
    s_dtheta_Metric.Metric[e_phi][e_t]       = s_dtheta_Metric.Metric[e_t][e_phi];
    s_dtheta_Metric.Metric[e_r][e_r]         = -2 * a * a * cos_theta * sin_theta / delta;
    s_dtheta_Metric.Metric[e_theta][e_theta] = -2 * a * a * cos_theta * sin_theta;
    s_dtheta_Metric.Metric[e_phi][e_phi]     = 4 * M * r * a * a * sin_theta * cos_theta / rho2 * (1 + a * a * sin_theta * sin_theta / rho2) * sin_theta * sin_theta +
                                               2 * s_Metric.Metric[e_phi][e_phi] / sin_theta * cos_theta;

    double sigma2 = rho2 * s_Metric.Metric[e_phi][e_phi] / sin_theta / sin_theta;
    double dtheta_sigma2 = -2 * a * a * delta * sin_theta * cos_theta;

    s_dtheta_Metric.Lapse_function = s_Metric.Lapse_function * a * a * sin_theta * cos_theta * (delta / sigma2 - 1 / rho2);
    s_dtheta_Metric.Shift_function = -s_Metric.Shift_function * dtheta_sigma2 / sigma2;

    return s_dtheta_Metric;
}

Metric_type Kerr_class::get_d2r_metric(const double* const State_Vector) const {

    Metric_type s_Metric = this->get_metric(State_Vector);
    Metric_type s_dr_Metric = this->get_dr_metric(State_Vector);

    const double& M = this->Mass;
    const double& a = this->Spin_Param;

    const double& r     = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r2 + a * a * cos_theta * cos_theta;
    double delta = r2 - 2 * M * r + a * a;

    Metric_type s_d2r_Metric{};

    /* --- Only the non-zero components are exlicitly evaluated. --- */

    s_d2r_Metric.Metric[e_t][e_t]         = 4 * M * r / rho2 / rho2 * (4 * r2 / rho2 - 3);
    s_d2r_Metric.Metric[e_t][e_phi]       = -4 * M * a * r * sin_theta * sin_theta / rho2 / rho2 * (4 * r2 / rho2 - 3);
    s_d2r_Metric.Metric[e_phi][e_t]       = s_d2r_Metric.Metric[e_t][e_phi];
    s_d2r_Metric.Metric[e_r][e_r]         = 2 / delta * (1 - 4 * (r2 - r * M) / delta + rho2 / delta * (4 * (r - M) * (r - M) / delta - 1));
    s_d2r_Metric.Metric[e_theta][e_theta] = 2.0;
    s_d2r_Metric.Metric[e_phi][e_phi]     = 2 * (1 + 2 * M * a * a * r / rho2 / rho2 * (4 * r2 / rho2 - 3) * sin_theta * sin_theta) * sin_theta * sin_theta;

    double sigma2 = rho2 * s_Metric.Metric[e_phi][e_phi] / sin_theta / sin_theta;
    double dr_sigma2 = (2 * r * s_Metric.Metric[e_phi][e_phi] + rho2 * s_dr_Metric.Metric[e_phi][e_phi]) / sin_theta / sin_theta;
    double d2r_sigma2 = (2 * s_Metric.Metric[e_phi][e_phi] + 4 * r * s_dr_Metric.Metric[e_phi][e_phi] + rho2 * s_d2r_Metric.Metric[e_phi][e_phi]) / sin_theta / sin_theta;

    double& N = s_Metric.Lapse_function;
    double& dr_N = s_dr_Metric.Lapse_function;
    s_d2r_Metric.Lapse_function = dr_N * dr_N / N + N / rho2 * (1 - 2 * r2 / rho2 + rho2 / delta * (1 - (r - M) * (r - M) / delta) - rho2 / sigma2 / 2 * (d2r_sigma2 - dr_sigma2 * dr_sigma2 / sigma2));

    double& omega = s_Metric.Shift_function;
    double& dr_omega = s_dr_Metric.Shift_function;
    s_d2r_Metric.Shift_function = -omega / r2 * (1 - r * dr_omega / omega + r * dr_sigma2 / sigma2) * (1 - r * dr_sigma2 / sigma2);

    return s_d2r_Metric;
}

int Kerr_class::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

    double M = this->Mass;
    double a = this->Spin_Param;

    double& r_obs = p_Initial_Conditions->Observer_params.distance;
    double& theta_obs = p_Initial_Conditions->Observer_params.inclination;

    p_Initial_Conditions->Init_Momentum[e_t] = -1;
    p_Initial_Conditions->Init_Momentum[e_phi] = -J_data[photon] * sin(theta_obs);
    p_Initial_Conditions->Init_Momentum[e_theta] = p_theta_data[photon];

    double& J = p_Initial_Conditions->Init_Momentum[e_phi];
    double& p_theta = p_Initial_Conditions->Init_Momentum[e_theta];

    double delta = pow(r_obs, 2) + pow(a, 2) - 2 * M * r_obs;
    double K = pow(p_theta, 2) + pow(cos(theta_obs), 2) * (pow(J / sin(theta_obs), 2) - pow(a, 2));

    double rad_potential = pow(r_obs * r_obs + a * a - a * J, 2) - delta * (pow(J - a, 2) + K);

    p_Initial_Conditions->Init_Momentum[e_r] = sqrt(rad_potential) / delta;

    return OK;
}

void Kerr_class::get_EOM(const double* const State_vector, double* const Derivatives) const {

    const double& r = State_vector[e_r];
    double r2 = r * r;

    const double& J = State_vector[e_p_phi];

    double sin1 = sin(State_vector[e_theta]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(State_vector[e_theta]);
    double cos2 = cos1 * cos1;

    double rho2 = r2 + this->Spin_Param * this->Spin_Param * cos2;

    double P = r2 + this->Spin_Param * this->Spin_Param - this->Spin_Param * J;
    double delta = r2 - 2 * this->Mass * r + this->Spin_Param * this->Spin_Param;
    double F = P * P - delta * ((J - this->Spin_Param) * (J - this->Spin_Param) + cos2 * (J * J / sin2 - this->Spin_Param * this->Spin_Param));

    const double& p_r     = State_vector[e_p_r];
    const double& p_theta = State_vector[e_p_theta];
    const double& p_t     = State_vector[e_p_t];

    Derivatives[e_t] = ((r2 + this->Spin_Param * this->Spin_Param) / delta * P - this->Spin_Param * (this->Spin_Param * sin2 - J)) / rho2;
    Derivatives[e_r] = delta / rho2 * p_r;
    Derivatives[e_theta] = 1.0 / rho2 * p_theta;
    Derivatives[e_phi] = 1.0 / (delta * rho2) * (P * this->Spin_Param + delta * (J / sin2 - this->Spin_Param));
    Derivatives[e_p_phi] = 0.0;
    Derivatives[e_p_t] = 0.0;

    double theta_term_1 = -(delta * p_r * p_r + p_theta * p_theta) * this->Spin_Param * this->Spin_Param * cos1 * sin1 / (rho2 * rho2);
    double theta_term_2 = F * this->Spin_Param * this->Spin_Param * cos1 * sin1 / (delta * rho2 * rho2) + (J * J * cos1 / (sin2 * sin1) - this->Spin_Param * this->Spin_Param * cos1 * sin1) / rho2;

    Derivatives[e_p_theta] = theta_term_1 + theta_term_2;

    double r_term_1 = p_r * p_r / (rho2) * (this->Mass - r * (1 - delta / rho2)) + p_theta * p_theta * r / (rho2 * rho2);
    double r_term_2 = (2 * P * r - (r - this->Mass) * ((J - this->Spin_Param) * (J - this->Spin_Param) + cos2 * (J * J / (sin2) - this->Spin_Param * this->Spin_Param))) / (delta * rho2)
                    - F * (rho2 * (r - this->Mass) + r * delta) / (delta * delta * rho2 * rho2);

    Derivatives[e_p_r] = r_term_1 + r_term_2;

}

bool Kerr_class::terminate_integration(const double* const State_vector, const double* const Derivatives) {

    bool scatter = State_vector[e_r] > 30 && Derivatives[e_r] < 0;

    double r_horizon = this->Mass * (1 + sqrt(1 - this->Spin_Param * this->Spin_Param));

    bool hit_horizon = State_vector[e_r] - r_horizon < 1e-4;

    return scatter || hit_horizon;

};

Return_Values Kerr_class::load_parameters(const Metric_parameters_type* const Metric_Parameters) {

    if (!isnan(Metric_Parameters->Spin)) {

        this->Spin_Param = Metric_Parameters->Spin;

        return OK;

    }

    return ERROR;

}

