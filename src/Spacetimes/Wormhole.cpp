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

Metric_type Wormhole_class::get_metric(const double* const State_Vector) {

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double exponent = -this->Mass / r - this->Redshift_Param * this->Mass * this->Mass / r2;

    memset(&this->s_Metric, 0, sizeof(this->s_Metric));

    this->s_Metric.Lapse_function = exp(exponent);
    this->s_Metric.Shift_function = 2 * this->Spin_Param * this->Mass * this->Mass / r2 / r;

    this->s_Metric.Metric[0][0] = -this->s_Metric.Lapse_function * this->s_Metric.Lapse_function +
        r2 * this->s_Metric.Shift_function * this->s_Metric.Shift_function * sin_theta * sin_theta;
    this->s_Metric.Metric[0][3] = -r2 * sin_theta * sin_theta * this->s_Metric.Shift_function;
    this->s_Metric.Metric[3][0] = this->s_Metric.Metric[0][3];
    this->s_Metric.Metric[1][1] = 1 / (1 - this->R_Throat / r);
    this->s_Metric.Metric[2][2] = r2;
    this->s_Metric.Metric[3][3] = r2 * sin_theta * sin_theta;

    return this->s_Metric;
}

Metric_type Wormhole_class::get_dr_metric(const double* const State_Vector) {

    Metric_type s_Metric = this->get_metric(State_Vector);

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    memset(&this->s_dr_Metric, 0, sizeof(this->s_dr_Metric));

    this->s_dr_Metric.Lapse_function = this->s_Metric.Lapse_function * (1 / r2 + 2 * this->Redshift_Param / (r2 * r));
    this->s_dr_Metric.Shift_function = -3 * this->s_Metric.Shift_function / r;

    double& N = this->s_Metric.Lapse_function;
    double& dr_N = this->s_dr_Metric.Lapse_function;
    double& omega = this->s_Metric.Shift_function;
    double& dr_omega = this->s_dr_Metric.Shift_function;

    this->s_dr_Metric.Metric[0][0] = -2 * N * dr_N + 2 * r * omega * (omega + r * dr_omega) * sin_theta * sin_theta;
    this->s_dr_Metric.Metric[0][3] = -r * (2 * omega + r * dr_omega) * sin_theta * sin_theta;
    this->s_dr_Metric.Metric[3][0] = this->s_dtheta_Metric.Metric[0][3];
    this->s_dr_Metric.Metric[1][1] = -1. / ((1 - this->R_Throat / r) * (1 - this->R_Throat / r)) * (this->R_Throat / r2);
    this->s_dr_Metric.Metric[2][2] = 2 * r;
    this->s_dr_Metric.Metric[3][3] = 2 * r * sin_theta * sin_theta;

    return this->s_dr_Metric;
}

Metric_type Wormhole_class::get_dtheta_metric(const double* const State_Vector) {

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    double exponent = -this->Mass / r - this->Redshift_Param * this->Mass * this->Mass / r2;

    memset(&this->s_dtheta_Metric, 0, sizeof(this->s_dtheta_Metric));

    this->s_dtheta_Metric.Lapse_function = 0.0;
    this->s_dtheta_Metric.Shift_function = 0.0;

    this->s_dtheta_Metric.Metric[0][0] = 2 * r2 * this->s_Metric.Shift_function * this->s_Metric.Shift_function * sin_theta * cos_theta;
    this->s_dtheta_Metric.Metric[0][3] = -2 * r2 * sin_theta * cos_theta * this->s_Metric.Shift_function;
    this->s_dtheta_Metric.Metric[3][0] = this->s_Metric.Metric[0][3];
    this->s_dtheta_Metric.Metric[1][1] = 0.0;
    this->s_dtheta_Metric.Metric[2][2] = 0.0;
    this->s_dtheta_Metric.Metric[3][3] = 2 * r2 * sin_theta * cos_theta;

    return this->s_Metric;
}

Metric_type Wormhole_class::get_d2r_metric(const double* const State_Vector) {

    Metric_type s_Metric = this->get_metric(State_Vector);
    Metric_type s_dr_Metric = this->get_dr_metric(State_Vector);

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double& N = this->s_Metric.Lapse_function;
    double& dr_N = this->s_dr_Metric.Lapse_function;
    double& omega = this->s_Metric.Shift_function;
    double& dr_omega = this->s_dr_Metric.Shift_function;

    memset(&this->s_d2r_Metric, 0, sizeof(this->s_d2r_Metric));

    this->s_d2r_Metric.Lapse_function = dr_N * (1 / r2 + 2 * this->Redshift_Param / (r2 * r)) - N * (2. / (r2 * r) + 6 * this->Redshift_Param / (r2 * r2));
    this->s_d2r_Metric.Shift_function = -3 * dr_omega / r + 3 * omega / r2;

    this->s_d2r_Metric.Metric[0][0] = -2 * dr_N * dr_N - 2 * N * this->s_d2r_Metric.Lapse_function + 2 * ((omega + r * dr_omega) * (omega + r * dr_omega) +
        r * omega * (dr_omega + dr_omega + r * this->s_d2r_Metric.Shift_function)) * sin_theta * sin_theta;
    this->s_d2r_Metric.Metric[0][3] = -(2 * omega + r * dr_omega + r * (3 * dr_omega + r * this->s_d2r_Metric.Shift_function)) * sin_theta * sin_theta;
    this->s_d2r_Metric.Metric[3][0] = this->s_d2r_Metric.Metric[0][3];
    this->s_d2r_Metric.Metric[1][1] = 2 / ((1 - this->R_Throat / r) * (1 - this->R_Throat / r)) * ((this->R_Throat / r2) * (this->R_Throat / r2) / (1 - this->R_Throat / r) + this->R_Throat / (r2 * r));
    this->s_d2r_Metric.Metric[2][2] = 2.0;
    this->s_d2r_Metric.Metric[3][3] = 2 * sin_theta * sin_theta;


    return this->s_d2r_Metric;
}

int Wormhole_class::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {


    double& r_obs = p_Initial_Conditions->Observer_params.distance;
    double& theta_obs = p_Initial_Conditions->Observer_params.inclination;

    p_Initial_Conditions->init_Three_Momentum[e_phi] = -J_data[photon] * sin(theta_obs);
    p_Initial_Conditions->init_Three_Momentum[e_theta] = p_theta_data[photon];

    double& p_theta = p_Initial_Conditions->init_Three_Momentum[e_theta];
    double& J = p_Initial_Conditions->init_Three_Momentum[e_phi];

    double& N = p_Initial_Conditions->init_metric_Redshift_func;
    double& omega = p_Initial_Conditions->init_metric_Shitft_func;

    double rad_potential = -(p_theta * p_theta + J * J / sin(theta_obs) / sin(theta_obs)) * N * N / r_obs / r_obs + (1 - omega * J) * (1 - omega * J);

    double(*metric)[4] = p_Initial_Conditions->init_metric;

    p_Initial_Conditions->init_Three_Momentum[e_r] = sqrt(rad_potential) / N / sqrt(metric[1][1]) * r_obs / sqrt(pow(r_obs, 2) - pow(1, 2)) * metric[1][1];

    return 0;
}

int Wormhole_class::get_EOM(double State_Vector[], double Derivatives[]) {

    double sqrt_r2 = sqrt(State_Vector[e_r] * State_Vector[e_r] + this->R_Throat * this->R_Throat);
    double d_ell_r = State_Vector[e_r] / sqrt_r2;

    double& J = State_Vector[e_p_phi];

    double omega = 2 * this->Spin_Param * this->Mass * this->Mass / sqrt_r2 / sqrt_r2 / sqrt_r2;
    double d_ell_omega = -3 * omega / sqrt_r2 * d_ell_r;

    double exponent = -1 / sqrt_r2 - this->Redshift_Param / (sqrt_r2 * sqrt_r2);
    double N = exp(-this->Mass / sqrt_r2 - this->Redshift_Param * this->Mass * this->Mass / sqrt_r2 / sqrt_r2);
    double d_ell_N = N * (1 / (sqrt_r2 * sqrt_r2) + 2 * this->Redshift_Param / (sqrt_r2 * sqrt_r2 * sqrt_r2)) * d_ell_r;

    double N2 = N * N;

    double sin1 = sin(State_Vector[e_theta]);
    double sin2 = sin1 * sin1;

    *(Derivatives + e_t) = -1.0 / N / N * State_Vector[e_p_r];
    *(Derivatives + e_r) = 1.0 / (1 + this->R_Throat / sqrt_r2) * State_Vector[e_p_r];
    *(Derivatives + e_theta) = 1.0 / (sqrt_r2 * sqrt_r2) * State_Vector[e_p_theta];
    *(Derivatives + e_phi) = J / (sqrt_r2 * sqrt_r2 * sin2) + omega * (1 - omega * J) / N2;
    *(Derivatives + e_p_phi) = 0.0;
    *(Derivatives + e_p_theta) = (cos(State_Vector[e_theta]) / sin1) / (sqrt_r2 * sqrt_r2) * J * J / sin2;
    *(Derivatives + e_p_t) = 0.0;

    double term_1 = -1.0 / ((1 + this->R_Throat / sqrt_r2) * (1 + this->R_Throat / sqrt_r2)) * this->R_Throat * State_Vector[e_r] / (sqrt_r2 * sqrt_r2 * sqrt_r2) * State_Vector[e_p_r] * State_Vector[e_p_r] / 2;
    double term_2 = 1.0 / (sqrt_r2 * sqrt_r2 * sqrt_r2) * (State_Vector[e_p_theta] * State_Vector[e_p_theta] + J * J / sin2) * d_ell_r;
    double term_3 = -(1.0 / (N2 * N) * d_ell_N * ((1 - omega * J) * (1 - omega * J)) - 1.0 / N2 * (-d_ell_omega * (1 - omega * J) * J));

    *(Derivatives + e_p_r) = term_1 + term_2 + term_3;

    return OK;
}

bool Wormhole_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool scatter            = State_vector[e_r] >  sqrt(100 * 100 + this->R_Throat * this->R_Throat) && Derivatives[e_r] < 0;
    bool scatter_other_side = State_vector[e_r] < -sqrt(100 * 100 + this->R_Throat * this->R_Throat);
    bool stop_at_throat     = State_vector[e_r] < 1e-5;

    if (this->Stop_at_Throat) {

        return scatter || stop_at_throat;
    }
    else {

        return scatter || scatter_other_side;

    }
};

bool Wormhole_class::load_parameters(Metric_parameters_type Metric_Parameters) {

    if (!isnan(Metric_Parameters.Spin) &&
        !isnan(Metric_Parameters.Redshift_Parameter)) {

        this->Spin_Param = Metric_Parameters.Spin;
        this->Redshift_Param = Metric_Parameters.Redshift_Parameter;
        this->Stop_at_Throat = Metric_Parameters.Stop_At_Throat;

        return true;

    }

    return false;

}

void Wormhole_class::update_parameters(double Param_value, Metric_Parameter_Selector Parameter) {


    if (Spin != Parameter && WH_Redshift != Parameter) {

        std::cout << "Wrong Parameter enum for sim mode 2! -> Can only be 'this->Spin_Param' or 'this->Redshift_Param'! Defaulting to 'this->Redshift_Param'..." << "\n";

        this->Redshift_Param = Param_value;

        return;

    }

    if (Spin == Parameter) {

        this->Spin_Param = Param_value;

    }
    else {

        this->Redshift_Param = Param_value;

    }

}