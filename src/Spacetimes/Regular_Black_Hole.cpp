#pragma once

#include <cmath>

#include "Structs.h"
#include "Spacetimes.h"
#include "Constants.h"

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

Metric_type RBH_class::get_metric(const double* const State_Vector) {

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double rho = sqrt(r2 + this->Parameter * this->Parameter);

    memset(&this->s_Metric, 0., sizeof(this->s_Metric));

    this->s_Metric.Metric[0][0] = -(1 - 2 * this->Mass / rho);
    this->s_Metric.Metric[0][3] = 0.0;
    this->s_Metric.Metric[1][1] = -1.0 / this->s_Metric.Metric[0][0];
    this->s_Metric.Metric[2][2] = rho * rho;
    this->s_Metric.Metric[3][3] = this->s_Metric.Metric[2][2] * sin_theta * sin_theta;

    this->s_Metric.Lapse_function = -this->s_Metric.Metric[0][0];
    this->s_Metric.Shift_function = 0.0;

    return this->s_Metric;
}


Metric_type RBH_class::get_dr_metric(const double* const State_Vector) {

    Metric_type s_Metric = this->get_metric(State_Vector);

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho = sqrt(r2 + this->Parameter * this->Parameter);
    double rho3 = rho * rho * rho;

    memset(&this->s_dr_Metric, 0., sizeof(this->s_dr_Metric));

    this->s_dr_Metric.Metric[0][0] = -2 * this->Mass * r / rho3;
    this->s_dr_Metric.Metric[0][3] = 0.0;
    this->s_dr_Metric.Metric[1][1] = 1.0 / (this->s_Metric.Metric[0][0] * this->s_Metric.Metric[0][0]) * this->s_dr_Metric.Metric[0][0];
    this->s_dr_Metric.Metric[2][2] = 2 * r;
    this->s_dr_Metric.Metric[3][3] = 2 * r * sin_theta * sin_theta;

    this->s_dr_Metric.Lapse_function = -this->s_dr_Metric.Metric[0][0];
    this->s_dr_Metric.Shift_function = 0.0;

    return this->s_dr_Metric;
}

Metric_type RBH_class::get_dtheta_metric(const double* const State_Vector) {

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    memset(&this->s_dtheta_Metric, 0., sizeof(this->s_dtheta_Metric));

    this->s_dtheta_Metric.Metric[0][0] = 0.0;
    this->s_dtheta_Metric.Metric[0][3] = 0.0;
    this->s_dtheta_Metric.Metric[3][0] = this->s_dtheta_Metric.Metric[0][3];
    this->s_dtheta_Metric.Metric[1][1] = 0.0;
    this->s_dtheta_Metric.Metric[2][2] = 0.0;
    this->s_dtheta_Metric.Metric[3][3] = 2 * r * r * sin_theta * cos_theta;

    this->s_dtheta_Metric.Lapse_function = 0.0;
    this->s_dtheta_Metric.Shift_function = 0.0;

    return this->s_dtheta_Metric;
}

Metric_type RBH_class::get_d2r_metric(const double* const State_Vector) {

    Metric_type s_Metric = this->get_metric(State_Vector);
    Metric_type s_dr_Metric = this->get_dr_metric(State_Vector);

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho = sqrt(r2 + this->Parameter * this->Parameter);
    double rho3 = rho * rho * rho;
    double rho5 = rho * rho * rho * rho * rho;

    memset(&this->s_d2r_Metric, 0., sizeof(this->s_d2r_Metric));

    this->s_d2r_Metric.Metric[0][0] = -2 * this->Mass / rho3 + 6 * this->Mass * r2 / (rho5);
    this->s_d2r_Metric.Metric[0][3] = 0.0;
    this->s_d2r_Metric.Metric[1][1] = 1.0 / (this->s_Metric.Metric[0][0] * this->s_Metric.Metric[0][0]) * this->s_d2r_Metric.Metric[0][0] -
        2.0 / (this->s_Metric.Metric[0][0] * this->s_Metric.Metric[0][0] * this->s_Metric.Metric[0][0]) * this->s_dr_Metric.Metric[0][0] * this->s_dr_Metric.Metric[0][0];
    this->s_d2r_Metric.Metric[2][2] = 2.0;
    this->s_d2r_Metric.Metric[3][3] = 2 * sin_theta * sin_theta;

    this->s_d2r_Metric.Lapse_function = -this->s_d2r_Metric.Metric[0][0];
    this->s_d2r_Metric.Shift_function = 0.0;

    return this->s_d2r_Metric;
}

int RBH_class::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

    double& r_obs = p_Initial_Conditions->Observer_params.distance;
    double& theta_obs = p_Initial_Conditions->Observer_params.inclination;

    p_Initial_Conditions->init_Three_Momentum[e_phi] = -J_data[photon] * sin(theta_obs);
    p_Initial_Conditions->init_Three_Momentum[e_theta] = p_theta_data[photon];

    double& J = p_Initial_Conditions->init_Three_Momentum[e_phi];

    double rho = sqrt(r_obs * r_obs + this->Parameter * this->Parameter);
    double rad_potential = 1 - (1 - 2 * this->Mass / rho) * J * J / (rho * rho);

    double (*metric)[4] = p_Initial_Conditions->init_metric;

    p_Initial_Conditions->init_Three_Momentum[e_r] = sqrt(rad_potential) * metric[1][1];

    return OK;
}

int RBH_class::get_EOM(double State_vector[], double Derivatives[]) {

    double r = State_vector[e_r];
    double rho = sqrt(r * r + this->Parameter * this->Parameter);

    double& J = State_vector[e_p_phi];

    double sin1 = sin(State_vector[e_theta]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(State_vector[e_theta]);
    double cos2 = cos1 * cos1;

    *(Derivatives + e_p_t    ) = - 1 / (1 - 2 * this->Mass / rho) * State_vector[e_p_t];
    *(Derivatives + e_r      ) = (1 - 2 * this->Mass / rho) * State_vector[e_p_r];
    *(Derivatives + e_theta  ) = 1.0 / (rho * rho) * State_vector[e_p_theta];
    *(Derivatives + e_phi    ) = J / (rho * rho * sin2);
    *(Derivatives + e_p_phi  ) = 0.0;
    *(Derivatives + e_p_theta) = cos1 / (rho * rho * sin1 * sin2) * J * J;

    double r_term_1 = -this->Mass * r / (rho * rho * rho) * (1.0 / ((1 - 2 * this->Mass / rho) * (1 - 2 * this->Mass / rho)) + State_vector[e_p_r] * State_vector[e_p_r]);
    double r_term_2 = r / (rho * rho * rho * rho) * (State_vector[e_p_theta] * State_vector[e_p_theta] + J * J / sin2);

    *(Derivatives + e_p_r) = r_term_1 + r_term_2;

    return 0;
}

bool RBH_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool scatter = State_vector[e_r] > 100 && Derivatives[e_r] < 0;

    double r_horizon = sqrt(4 * this->Mass * this->Mass - this->Parameter * this->Parameter);

    bool hit_horizon_RBH = State_vector[e_r] - r_horizon < 0.05;

    return scatter || hit_horizon_RBH;

}

bool RBH_class::load_parameters(Metric_parameters_type Metric_Parameters) {

    if (!isnan(Metric_Parameters.RBH_Parameter)) {

        this->Parameter = Metric_Parameters.RBH_Parameter;

        return true;

    }

    return false;

}

void RBH_class::update_parameters(double Param_value, Metric_Parameter_Selector Parameter) {


    if (RBH_Param != Parameter) {

        std::cout << "Wrong Parameter enum for sim mode 2! -> Can only be 'RBH_Param'! Defaulting to it..." << "\n";

    }

    this->Parameter = Param_value;

}
