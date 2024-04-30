#pragma once

#include <cmath>

#include "Structs.h"
#include "Spacetimes.h"
#include "Constants.h"

double* JNW_class::get_ISCO() {

    static double r_ISCO[2]{};

    if (JNW_GAMMA > 1.0 / 2) {

        r_ISCO[Inner] = 1.0 / JNW_GAMMA * (3.0 * JNW_GAMMA + 1.0 + sqrt(5 * JNW_GAMMA * JNW_GAMMA - 1));
        r_ISCO[Outer] = r_ISCO[Inner];


    }else if (JNW_GAMMA > 1.0 / sqrt(5) && JNW_GAMMA < 1.0 / 2) {

            r_ISCO[Inner] = 1.0 / JNW_GAMMA * (3.0 * JNW_GAMMA + 1.0 - sqrt(5 * JNW_GAMMA * JNW_GAMMA - 1));
            r_ISCO[Outer] = 1.0 / JNW_GAMMA * (3.0 * JNW_GAMMA + 1.0 + sqrt(5 * JNW_GAMMA * JNW_GAMMA - 1));

    }
    else {

        r_ISCO[Inner] = JNW_R_SINGULARITY;
        r_ISCO[Outer] = JNW_R_SINGULARITY;

    }

    return r_ISCO;

 };

double* JNW_class::get_Photon_Sphere() {

    static double photon_orbit[2]{};

    if (JNW_GAMMA > 0.5) { // Weak naked singularity

        photon_orbit[Inner] = (2 * JNW_GAMMA + 1) * JNW_R_SINGULARITY / 2;
       

    }
    else {

        photon_orbit[Inner] = JNW_R_SINGULARITY;

    }

    photon_orbit[Outer] = photon_orbit[Inner];

    return photon_orbit;

};

int JNW_class::update_metric(double State_Vector[]) {

    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    this->s_Metric.Metric[0][0] = -pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA);
    this->s_Metric.Metric[1][1] = -1.0 / this->s_Metric.Metric[0][0];
    this->s_Metric.Metric[2][2] = pow(1 - JNW_R_SINGULARITY / r, 1 - JNW_GAMMA) * r2;
    this->s_Metric.Metric[3][3] = this->s_Metric.Metric[2][2] * sin_theta * sin_theta;
    this->s_Metric.Metric[0][3] = 0.;
    this->s_Metric.Metric[3][0] = this->s_Metric.Metric[0][3];

    this->s_Metric.Lapse_function = -this->s_Metric.Metric[0][0];
    this->s_Metric.Shift_function = 0.;

    return OK;

}

Metric_type JNW_class::get_metric(double State_Vector[]) {

    this->update_metric(State_Vector);

    return this->s_Metric;

}

int JNW_class::update_dr_metric(double State_Vector[]) {

    Metric_type Metric = this->get_metric(State_Vector);

    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    this->s_dr_Metric.Metric[0][0] = -JNW_GAMMA * pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA - 1) * JNW_R_SINGULARITY / r2;
    this->s_dr_Metric.Metric[1][1] = 1.0 / (this->s_Metric.Metric[0][0] * this->s_Metric.Metric[0][0]) * this->s_dr_Metric.Metric[0][0];;
    this->s_dr_Metric.Metric[2][2] = 2 * r * pow(1 - JNW_R_SINGULARITY / r, 1 - JNW_GAMMA) + (1 - JNW_GAMMA) * pow(1 - JNW_R_SINGULARITY / r, -JNW_GAMMA) * JNW_R_SINGULARITY;
    this->s_dr_Metric.Metric[3][3] = this->s_dr_Metric.Metric[2][2] * sin_theta * sin_theta;
    this->s_dr_Metric.Metric[0][3] = 0;
    this->s_dr_Metric.Metric[3][0] = this->s_dr_Metric.Metric[0][3];

    this->s_dr_Metric.Lapse_function = -this->s_dr_Metric.Metric[0][0];
    this->s_dr_Metric.Shift_function = 0.0;

    return OK;

};

Metric_type JNW_class::get_dr_metric(double State_Vector[]) {

    this->update_dr_metric(State_Vector);

    return this->s_dr_Metric;

}

int JNW_class::update_dtheta_metric(double State_Vector[]) {

    double& r = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    this->s_dtheta_Metric.Metric[0][0] = 0.0;
    this->s_dtheta_Metric.Metric[0][3] = 0.0;
    this->s_dtheta_Metric.Metric[3][0] = this->s_dtheta_Metric.Metric[0][3];
    this->s_dtheta_Metric.Metric[1][1] = 0.0;
    this->s_dtheta_Metric.Metric[2][2] = 0.0;
    this->s_dtheta_Metric.Metric[3][3] = 2 * r * r * sin_theta * cos_theta;

    this->s_dtheta_Metric.Lapse_function = 0.0;
    this->s_dtheta_Metric.Shift_function = 0.0;

    return OK;

}

Metric_type JNW_class::get_dtheta_metric(double State_Vector[]) {

    this->update_dtheta_metric(State_Vector);

    return this->s_dtheta_Metric;
}

int JNW_class::update_d2r_metric(double State_Vector[]) {

    Metric_type s_Metric = this->get_metric(State_Vector);
    Metric_type s_dr_Metric = this->get_dr_metric(State_Vector);

    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    this->s_d2r_Metric.Metric[0][0] = -JNW_GAMMA * (JNW_GAMMA - 1) * pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA - 2) * JNW_R_SINGULARITY * JNW_R_SINGULARITY / r2 / r2
                                    + 2 * JNW_GAMMA * pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA - 1) * JNW_R_SINGULARITY / r2 / r;

    this->s_d2r_Metric.Metric[0][3] = 0.0;

    this->s_d2r_Metric.Metric[1][1] = 1.0 / (this->s_Metric.Metric[0][0] * this->s_Metric.Metric[0][0]) * this->s_d2r_Metric.Metric[0][0]
                     - 2.0 / (this->s_Metric.Metric[0][0] * this->s_Metric.Metric[0][0] * this->s_Metric.Metric[0][0]) * this->s_dr_Metric.Metric[0][0] * this->s_dr_Metric.Metric[0][0];

    this->s_d2r_Metric.Metric[2][2] = 2 * pow(1 - JNW_R_SINGULARITY / r, 1 - JNW_GAMMA) + 2 * (1 - JNW_GAMMA) * pow(1 - JNW_R_SINGULARITY / r, -JNW_GAMMA) * JNW_R_SINGULARITY / r
                     - (1 - JNW_GAMMA) * JNW_GAMMA * pow(1 - JNW_R_SINGULARITY / r, -JNW_GAMMA - 1) * JNW_R_SINGULARITY * JNW_R_SINGULARITY / r2;

    this->s_d2r_Metric.Metric[3][3] = this->s_d2r_Metric.Metric[2][2] * sin_theta * sin_theta;

    this->s_d2r_Metric.Lapse_function = -this->s_d2r_Metric.Metric[0][0];
    this->s_d2r_Metric.Shift_function = 0.0;

    return OK;

}

Metric_type JNW_class::get_d2r_metric(double State_Vector[]) {

    this->update_d2r_metric(State_Vector);

    return this->s_d2r_Metric;

}

int JNW_class::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

    double& r_obs = p_Initial_Conditions->init_Pos[e_r];
    double& theta_obs = p_Initial_Conditions->init_Pos[e_theta];

    p_Initial_Conditions->init_Three_Momentum[e_phi] = -J_data[photon] * sin(theta_obs);
    p_Initial_Conditions->init_Three_Momentum[e_theta] = p_theta_data[photon];

    double& J = p_Initial_Conditions->init_Three_Momentum[e_phi];

    double rad_potential = 1 - pow(1 - JNW_R_SINGULARITY / r_obs, 2 * JNW_GAMMA - 1) * J * J / (r_obs * r_obs);

    double(*metric)[4] = p_Initial_Conditions->init_metric;

    p_Initial_Conditions->init_Three_Momentum[e_r] = sqrt(rad_potential) * metric[1][1];

    return OK;

}

int JNW_class::get_EOM(double State_vector[], double Derivatives[])
{

    double& r = State_vector[e_r];
    double& J = State_vector[e_p_phi];

    double sin1 = sin(State_vector[e_theta]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(State_vector[e_theta]);
    double cos2 = cos1 * cos1;

    double pow_gamma = pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA);
    double pow_gamma_minus_1 = pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA - 1);

    *(Derivatives + e_r      ) = pow_gamma * State_vector[e_p_r];
    *(Derivatives + e_theta  ) = pow_gamma_minus_1 / (r * r) * State_vector[e_p_theta];
    *(Derivatives + e_phi    ) = pow_gamma_minus_1 / (r * r * sin2) * J;
    *(Derivatives + e_p_phi  ) = 0.0;
    *(Derivatives + e_p_theta) = pow_gamma_minus_1 * cos1 / (r * r * sin1 * sin2) * J * J;

    double r_term_1 = -JNW_GAMMA * JNW_R_SINGULARITY / 2 / r / r * pow_gamma_minus_1 * (1.0 / pow_gamma / pow_gamma
                    + State_vector[e_p_r] * State_vector[e_p_r]);
    double r_term_2 = 1.0 / r / r / r * pow_gamma_minus_1 * (1 - JNW_R_SINGULARITY / 2 / r * (JNW_GAMMA - 1) / (1 - JNW_R_SINGULARITY / r))
                    * (State_vector[e_p_theta] * State_vector[e_p_theta] + J * J / sin2);

    *(Derivatives + e_p_r) = r_term_1 + r_term_2;

    return OK;

}

bool JNW_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool hit_singularity = State_vector[e_r] - JNW_R_SINGULARITY < 1e-4;

    bool scatter = State_vector[e_r] > 100 && Derivatives[e_r] < 0;

    if (JNW_GAMMA > 0.5) {

        return scatter || hit_singularity;

    }
    else {

        return scatter;

    }
};

Metric_Parameters_type JNW_class::get_parameters() {

    Metric_Parameters_type Parameters{};

    Parameters.JNW_Gamma_Parameter = this->Gamma;

    return Parameters;

}

bool JNW_class::load_parameters(Metric_Parameters_type Metric_Parameters) {

    if (!isnan(Metric_Parameters.JNW_Gamma_Parameter)) {

        this->Gamma = Metric_Parameters.JNW_Gamma_Parameter;

        return true;

    }

    return false;

}

void JNW_class::update_parameters(double Param_value, Metric_Parameter_Selector Parameter) {


    if (JNW_Gamma != Parameter) {

        std::cout << "Wrong Parameter enum for sim mode 2! -> Can only be 'JNW_Gamma'! Defaulting to it..." << "\n";

    }

    this->Gamma = Param_value;

}
