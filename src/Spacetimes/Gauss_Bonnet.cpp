#include "Spacetimes.h"

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

Metric_type Gauss_Bonnet_class::get_metric(const double* const State_Vector) const {

    const double& M = this->Mass;
    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double f = 1. + r2 / this->Gamma / 2. * (1. - sqrt(1. + 8. * this->Gamma * M / r2 / r));

    Metric_type s_Metric{};

    /* --- Only the non-zero components are exlicitly evaluated. --- */

    s_Metric.Metric[e_t][e_t]         = -f;
    s_Metric.Metric[e_r][e_r]         = 1. / f;
    s_Metric.Metric[e_theta][e_theta] = r2;
    s_Metric.Metric[e_phi][e_phi]     = r2 * sin_theta * sin_theta;

    s_Metric.Lapse_function = sqrt(-s_Metric.Metric[e_t][e_t]);

    return s_Metric;

}

Metric_type Gauss_Bonnet_class::get_dr_metric(const double* const State_Vector) const {

    const double& M = this->Mass;
    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double f = 1. + r2 / this->Gamma / 2. * (1. - sqrt(1. + 8. * this->Gamma * M / r2 / r));
    double dr_f = 2. / r * (f - 1.) + 6. * M / sqrt(r2 * r2 + 8. * this->Gamma * M * r);

    Metric_type s_dr_Metric{};

    /* --- Only the non-zero components are exlicitly evaluated. --- */

    s_dr_Metric.Metric[e_t][e_t]         = -dr_f;
    s_dr_Metric.Metric[e_r][e_r]         = -1. / f / f * dr_f;
    s_dr_Metric.Metric[e_theta][e_theta] = 2. * r;
    s_dr_Metric.Metric[e_phi][e_phi]     = 2. * r * sin_theta * sin_theta;

    s_dr_Metric.Lapse_function = 0;

    return s_dr_Metric;

}

Metric_type Gauss_Bonnet_class::get_dtheta_metric(const double* const State_Vector) const {

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    Metric_type s_dtheta_Metric{};

    /* --- Only the non-zero components are exlicitly evaluated. --- */

    s_dtheta_Metric.Metric[e_phi][e_phi] = 2 * r * r * sin_theta * cos_theta;

    return s_dtheta_Metric;
}

Metric_type Gauss_Bonnet_class::get_d2r_metric(const double* const State_Vector) const {

    const double& M = this->Mass;
    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double root = sqrt(r2 * r2 + 8 * this->Gamma * M * r);

    double f = 1 + r2 / this->Gamma / 2. * (1 - sqrt(1. + 8. * this->Gamma * M / r2 / r));
    double dr_f = 2. / r * (f - 1.) + 6 * M / root;
    double d2r_f = -2. / r2 * (f - 1.) + 2. / r * dr_f - 12. * M / root / root / root * (r2 * r + 2. * this->Gamma * M);

    Metric_type s_d2r_Metric{};

    /* --- Only the non-zero components are exlicitly evaluated. --- */

    s_d2r_Metric.Metric[e_t][e_t]         = -d2r_f;
    s_d2r_Metric.Metric[e_r][e_r]         = 2. / f / f / f * dr_f - 1. / f / f * d2r_f;
    s_d2r_Metric.Metric[e_theta][e_theta] = 2.;
    s_d2r_Metric.Metric[e_phi][e_phi]     = 2. * sin_theta * sin_theta;

    s_d2r_Metric.Lapse_function = -s_d2r_Metric.Metric[e_t][e_t];

    return s_d2r_Metric;

}

int Gauss_Bonnet_class::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

    double& r_obs     = p_Initial_Conditions->Observer_params.distance;
    double& theta_obs = p_Initial_Conditions->Observer_params.inclination;

    p_Initial_Conditions->Init_Momentum[e_phi] = -J_data[photon] * sin(theta_obs);
    p_Initial_Conditions->Init_Momentum[e_theta] = p_theta_data[photon];

    double& J = p_Initial_Conditions->Init_Momentum[e_phi];
    double  f = 1. + r_obs * r_obs / 2. / this->Gamma * (1. - sqrt(1. + 8. * this->Gamma * this->Mass / r_obs / r_obs / r_obs));

    double rad_potential = 1. - f * J * J / (r_obs * r_obs);

    double(*metric)[4] = p_Initial_Conditions->Init_metric.Metric;

    p_Initial_Conditions->Init_Momentum[e_r] = sqrt(rad_potential) * metric[1][1];

    return OK;

}

void Gauss_Bonnet_class::get_EOM(double State_vector[], double Derivatives[]) const{

    double& r = State_vector[e_r];
    double& J = State_vector[e_p_phi];

    double sin1 = sin(State_vector[e_theta]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(State_vector[e_theta]);
    double cos2 = cos1 * cos1;

    double root = sqrt(1. + 8. * this->Gamma * this->Mass / r / r / r);

    double f    = 1. + r * r / this->Gamma / 2. * (1. - root);
    double dr_f = 2. / r * (f - 1.) + 6. * this->Mass / root / r / r;

    *(Derivatives + e_t      ) = - 1 / f * State_vector[e_p_t];
    *(Derivatives + e_r      ) = f * State_vector[e_p_r];
    *(Derivatives + e_theta  ) = 1. / (r * r) * State_vector[e_p_theta];
    *(Derivatives + e_phi    ) = J / (r * r * sin2);
    *(Derivatives + e_p_phi  ) = 0.0;
    *(Derivatives + e_p_theta) = cos1 / (r * r * sin1 * sin2) * J * J;

    double r_term_1 = -1. / 2 * (1.0 / f / f + State_vector[e_p_r] * State_vector[e_p_r]) * dr_f;
    double r_term_2 = 1.0 / r / r / r * (State_vector[e_p_theta] * State_vector[e_p_theta] + J * J / sin2);

    *(Derivatives + e_p_r) = r_term_1 + r_term_2;
}

bool Gauss_Bonnet_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool scatter = State_vector[e_r] > 30 && Derivatives[e_r] < 0;

    bool too_high_order = State_vector[e_phi] * State_vector[e_phi] > 5 * M_PI * 5 * M_PI;

    double r_horizon{};
    bool hit_horizon = false;

    if (this->Gamma <= 1) {

        r_horizon = this->Mass + sqrt(this->Mass * this->Mass - this->Gamma * this->Gamma);
        hit_horizon = State_vector[e_r] - r_horizon < 1e-5;

    }

    return scatter || too_high_order || hit_horizon;

};

Return_Values Gauss_Bonnet_class::load_parameters(const Metric_parameters_type* const Metric_Parameters) {

    if (!isnan(Metric_Parameters->GB_Gamma_Parameter)) {

        this->Gamma = Metric_Parameters->GB_Gamma_Parameter;

        return OK;

    }

    return ERROR;

}
