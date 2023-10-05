#pragma once

#include <cmath>
#include <iostream>
#include <vector>

#include "Constants.h"
#include "Enumerations.h"
#include "Spacetimes.h"
#include "General_GR_functions.h"

/***************************************
|                                      |
| Observer Class Functions Definitions |
|                                      |
***************************************/

extern Spacetime_Base_Class* Spacetimes[];

Observer_class::Observer_class(double r, double theta, double phi) {

    r_obs     = r;
    theta_obs = theta;
    phi_obs   = phi;

    double Obs_State_Vector[3] = {r, theta, phi};

    Metric_type obs_metric = Spacetimes[e_metric]->get_metric(Obs_State_Vector);

    /*

    Obs_velocity is given in contravatiant components

    */

    obs_velocity[0] = 1.0 / obs_metric.Lapse_function;
    obs_velocity[1] = 0;
    obs_velocity[2] = 0;
    obs_velocity[3] = obs_metric.Shift_function / obs_metric.Lapse_function;

}

double Observer_class::get_r_obs()     { return r_obs; };
double Observer_class::get_theta_obs() { return theta_obs; };
double Observer_class::get_phi_obs()   { return phi_obs; };

int Observer_class::get_obs_velocity(double Obs_velocity[4]) {

    for (int index = 0; index <= 3; index++) {

        Obs_velocity[index] = obs_velocity[index];

    }

    return OK;

}

/*******************************************
|                                          |
| Derived Kerr Class Functions Definitions |
|                                          |
*******************************************/

double* Kerr_class::get_ISCO() {

    double Z_1 = 1 + pow(1 - SPIN * SPIN, 1. / 3) * (pow(1 + SPIN, 1. / 3) + pow(1 - SPIN, 1. / 3));
    double Z_2 = sqrt(3 * SPIN * SPIN + Z_1 * Z_1);

    static double r_ISCO[2]{};

    r_ISCO[Inner] = MASS * (3 + Z_2 - sqrt((3 - Z_1) * (3 + Z_1 + 2 * Z_2)));
    r_ISCO[Outer] = MASS * (3 + Z_2 + sqrt((3 - Z_1) * (3 + Z_1 + 2 * Z_2)));

    return r_ISCO;

}

double* Kerr_class::get_Photon_Sphere() {

    static double photon_orbit[2]{};

    photon_orbit[Inner] = 2 * MASS * (1 + cos(2.0 / 3 * acos(SPIN)));
    photon_orbit[Outer] = 2 * MASS * (1 + cos(2.0 / 3 * acos(-SPIN)));

    return photon_orbit;

}

int Kerr_class::update_metric(double State_Vector[]) {

    double M = MASS;
    double a = SPIN;

    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r2 + a * a * cos_theta * cos_theta;
    double delta = r2 - 2 * M * r + a * a;

    this->s_Metric.Metric[0][0] = -(1 - 2 * M * r / rho2);
    this->s_Metric.Metric[0][3] = -2 * M * r * a * sin_theta * sin_theta / rho2;
    this->s_Metric.Metric[3][0] = this->s_Metric.Metric[0][3];
    this->s_Metric.Metric[1][1] = rho2 / delta;
    this->s_Metric.Metric[2][2] = rho2;
    this->s_Metric.Metric[3][3] = (r2 + a * a + 2 * M * r * a * a / rho2 * sin_theta * sin_theta) * sin_theta * sin_theta;

    double sigma2 = (r2 + a * a) * (r2 + a * a) - a * a * delta * sin_theta * sin_theta;

    this->s_Metric.Lapse_function = sqrt(rho2 * delta / sigma2);
    this->s_Metric.Shift_function = 2 * a * r / sigma2;

    return OK;

};

Metric_type Kerr_class::get_metric(double State_vector[]) {

    if (false == this->eval_bitmask[0] || true == this->ignore_flag) {

        this->update_metric(State_vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[0] = true;

        }

    }

    return this->s_Metric;
}

int Kerr_class::update_dr_metric(double State_Vector[]) {

    double M = MASS;
    double a = SPIN;

    Metric_type s_Metric = this->get_metric(State_Vector);

    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r2 + a * a * cos_theta * cos_theta;
    double delta = r2 - 2 * M * r + a * a;

    this->s_dr_Metric.Metric[0][0] = -2 * M / rho2 * (2 * r2 / rho2 - 1);
    this->s_dr_Metric.Metric[0][3] = 2 * M * a * sin_theta * sin_theta / rho2 * (2 * r2 / rho2 - 1);
    this->s_dr_Metric.Metric[3][0] = this->s_dr_Metric.Metric[0][3];
    this->s_dr_Metric.Metric[1][1] = 2 * r / delta * (1 - rho2 / delta * (1 - M / r));
    this->s_dr_Metric.Metric[2][2] = 2 * r;
    this->s_dr_Metric.Metric[3][3] = 2 * (r - M * a * a / rho2 * (2 * r2 / rho2 - 1) * sin_theta * sin_theta) * sin_theta * sin_theta;

    double sigma2 = rho2 * this->s_Metric.Metric[3][3] / sin_theta / sin_theta;
    double dr_sigma2 = 2 * r * this->s_Metric.Metric[3][3] + rho2 * this->s_dr_Metric.Metric[3][3];

    this->s_dr_Metric.Lapse_function = this->s_Metric.Lapse_function * (r / rho2 + (r - M) / delta - dr_sigma2 / 2 / sigma2);
    this->s_dr_Metric.Shift_function = this->s_Metric.Shift_function / r * (1 - r * dr_sigma2 / sigma2);

    return OK;

}

Metric_type Kerr_class::get_dr_metric(double State_Vector[]) {

    if (false == this->eval_bitmask[1] || true == this->ignore_flag) {

        this->update_dr_metric(State_Vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[1] = true;

        }

    }

    return this->s_dr_Metric;
}

int Kerr_class::update_d2r_metric(double State_Vector[]) {

    double M = MASS;
    double a = SPIN;

    Metric_type s_Metric    = this->get_metric(State_Vector);
    Metric_type s_dr_Metric = this->get_dr_metric(State_Vector);

    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r2 + a * a * cos_theta * cos_theta;
    double delta = r2 - 2 * M * r + a * a;

    this->s_d2r_Metric.Metric[0][0] = 4 * M * r / rho2 / rho2 * (4 * r2 / rho2 - 3);
    this->s_d2r_Metric.Metric[0][3] = -4 * M * a * r * sin_theta * sin_theta / rho2 / rho2 * (4 * r2 / rho2 - 3);
    this->s_d2r_Metric.Metric[3][0] = this->s_d2r_Metric.Metric[0][3];
    this->s_d2r_Metric.Metric[1][1] = 2 / delta * (1 - 4 * (r2 - r * M) / delta + rho2 / delta * (2 * (r - M) * (r - M) / delta - 1));
    this->s_d2r_Metric.Metric[2][2] = 2.0;
    this->s_d2r_Metric.Metric[3][3] = 2 * (1 + 2 * M * a * a * r / rho2 / rho2 * (4 * r2 / rho2 - 3) * sin_theta * sin_theta) * sin_theta * sin_theta;

    double sigma2     = rho2 * this->s_Metric.Metric[3][3] / sin_theta / sin_theta;
    double dr_sigma2  = (2 * r * this->s_Metric.Metric[3][3] + rho2 * this->s_dr_Metric.Metric[3][3]) / sin_theta / sin_theta;
    double d2r_sigma2 = (2 * this->s_Metric.Metric[3][3] + 4 * r * this->s_dr_Metric.Metric[3][3] + rho2 * this->s_d2r_Metric.Metric[3][3]) / sin_theta / sin_theta;

    double& N    = this->s_Metric.Lapse_function;
    double& dr_N = this->s_dr_Metric.Lapse_function;
    this->s_d2r_Metric.Lapse_function = dr_N * dr_N / N + N / rho2 * (1 - 2 * r2 / rho2 + rho2 / delta * (1 - (r - M) * (r - M) / delta) - rho2 / sigma2 / 2 * (d2r_sigma2 - dr_sigma2 * dr_sigma2 / sigma2));

    double& omega    = this->s_Metric.Shift_function;
    double& dr_omega = this->s_dr_Metric.Shift_function;
    this->s_d2r_Metric.Shift_function = -omega / r2 * (1 - r * dr_omega / omega + r * dr_sigma2 / sigma2) * (1 - r * dr_sigma2 / sigma2);

    return OK;

}

Metric_type Kerr_class::get_d2r_metric(double State_Vector[]) {

    if (false == this->eval_bitmask[2] || true == this->ignore_flag) {

        this->update_d2r_metric(State_Vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[2] = true;

        }

    }

    return this->s_d2r_Metric;
}


int Kerr_class::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

    double M = MASS;
    double a = SPIN;

    double& r_obs = p_Initial_Conditions->init_Pos[e_r];
    double& theta_obs = p_Initial_Conditions->init_Pos[e_theta];

    p_Initial_Conditions->init_Three_Momentum[e_phi] = -J_data[photon] * sin(theta_obs);
    p_Initial_Conditions->init_Three_Momentum[e_theta] = p_theta_data[photon];

    double& J = p_Initial_Conditions->init_Three_Momentum[e_phi];
    double& p_theta = p_Initial_Conditions->init_Three_Momentum[e_theta];

    double delta = pow(r_obs, 2) + pow(a, 2) - 2 * M * r_obs;
    double K = pow(p_theta, 2) + pow(cos(theta_obs), 2) * (pow(J / sin(theta_obs), 2) - pow(a, 2));

    double rad_potential = pow(r_obs * r_obs + a * a - a * J, 2) - delta * (pow(J - a, 2) + K);

    p_Initial_Conditions->init_Three_Momentum[e_r] = sqrt(rad_potential) / delta;

    return OK;
}

int Kerr_class::get_EOM(double inter_State_vector[], double Derivatives[], int iteration) {

    double& r = inter_State_vector[e_r + iteration * e_State_Number];
    double r2 = r * r;

    double& J = inter_State_vector[e_p_phi];

    double sin1 = sin(inter_State_vector[e_theta + iteration * e_State_Number]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(inter_State_vector[e_theta + iteration * e_State_Number]);
    double cos2 = cos1 * cos1;

    double rho2 = r2 + SPIN * SPIN * cos2;

    double P = r2 + SPIN * SPIN - SPIN * J;
    double delta = r2 - 2 * MASS * r + SPIN * SPIN;
    double F = P * P - delta * ((J - SPIN) * (J - SPIN) + cos2 * (J * J / sin2 - SPIN * SPIN));

    double& p_r     = inter_State_vector[e_p_r     + iteration * e_State_Number];
    double& p_theta = inter_State_vector[e_p_theta + iteration * e_State_Number];

    Derivatives[e_r      + iteration * e_State_Number] = delta / rho2 * p_r;
    Derivatives[e_theta  + iteration * e_State_Number] = 1.0 / rho2 * p_theta;
    Derivatives[e_phi    + iteration * e_State_Number] = 1.0 / (delta * rho2) * (P * SPIN + delta * (J / sin2 - SPIN));
    Derivatives[e_p_phi  + iteration * e_State_Number] = 0.0;

    double theta_term_1 = -(delta * p_r * p_r + p_theta * p_theta) * SPIN * SPIN * cos1 * sin1 / (rho2 * rho2);
    double theta_term_2 = F * SPIN * SPIN * cos1 * sin1 / (delta * rho2 * rho2) + (J * J * cos1 / (sin2 * sin1) - SPIN * SPIN * cos1 * sin1) / rho2;

    Derivatives[e_p_theta + iteration * e_State_Number] = theta_term_1 + theta_term_2;

    double r_term_1 = p_r * p_r / (rho2) * (MASS - r * (1 - delta / rho2)) + p_theta * p_theta * r / (rho2 * rho2);
    double r_term_2 = (2 * P * r - (r - MASS) * ((J - SPIN) * (J - SPIN) + cos2 * (J * J / (sin2) - SPIN * SPIN))) / (delta * rho2)
                    - F * (rho2 * (r - MASS) + r * delta) / (delta * delta * rho2 * rho2);

    Derivatives[e_p_r + iteration * e_State_Number] = r_term_1 + r_term_2;

    return OK;

}

bool Kerr_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool scatter = State_vector[e_r] > 30 && Derivatives[e_r] < 0;

    double r_horizon = MASS * (1 + sqrt(1 - SPIN * SPIN));

    bool hit_horizon = State_vector[e_r] - r_horizon < 1e-5;

    return scatter || hit_horizon;

};

void Kerr_class::reset_eval_bitmask() {

    for (int index = 0; index <= 2; index++) {

        this->eval_bitmask[index] = false;

    }

}

void Kerr_class::set_ignore_flag(bool flag) {

    this->ignore_flag = flag;

};

/*********************************************************
|                                                        |
| Derived Regular Black Hole Class Functions Definitions |
|                                                        |
*********************************************************/

double* RBH_class::get_ISCO() {

    static double r_ISCO[2]{};

    r_ISCO[Inner] = sqrt(36 * MASS * MASS - RBH_PARAM * RBH_PARAM);
    r_ISCO[Outer] = r_ISCO[Inner];

    return r_ISCO;

}

double* RBH_class::get_Photon_Sphere() {

    double M = MASS;

    static double photon_orbit[2]{};

    photon_orbit[Inner] = sqrt(9 * MASS * MASS - RBH_PARAM * RBH_PARAM);
    photon_orbit[Outer] = photon_orbit[Inner];

    return photon_orbit;

}

int RBH_class::update_metric(double State_Vector[]) {

    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double rho = sqrt(r2 + RBH_PARAM * RBH_PARAM);

    this->s_Metric.Metric[0][0] = -(1 - 2 * MASS / rho);
    this->s_Metric.Metric[0][3] = 0.0;
    this->s_Metric.Metric[1][1] = -1.0 / this->s_Metric.Metric[0][0];
    this->s_Metric.Metric[2][2] = rho * rho;
    this->s_Metric.Metric[3][3] = this->s_Metric.Metric[2][2] * sin_theta * sin_theta;

    this->s_Metric.Lapse_function = -this->s_Metric.Metric[0][0];
    this->s_Metric.Shift_function = 0.0;

    return OK;
}

Metric_type RBH_class::get_metric(double State_vector[]) {

    if (false == this->eval_bitmask[0] || true == this->ignore_flag) {

        this->update_metric(State_vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[0] = true;

        }

    }

    return this->s_Metric;
}

int RBH_class::update_dr_metric(double State_Vector[]) {

    Metric_type s_Metric = this->get_metric(State_Vector);

    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho = sqrt(r2 + RBH_PARAM * RBH_PARAM);
    double rho3 = rho * rho * rho;

    this->s_dr_Metric.Metric[0][0] = -2 * MASS * r / rho3;
    this->s_dr_Metric.Metric[0][3] = 0.0;
    this->s_dr_Metric.Metric[1][1] = 1.0 / (this->s_Metric.Metric[0][0] * this->s_Metric.Metric[0][0]) * this->s_dr_Metric.Metric[0][0];
    this->s_dr_Metric.Metric[2][2] = 2 * r;
    this->s_dr_Metric.Metric[3][3] = 2 * r * sin_theta * sin_theta;

    this->s_dr_Metric.Lapse_function = -this->s_dr_Metric.Metric[0][0];
    this->s_dr_Metric.Shift_function = 0.0;

    return OK;

}

Metric_type RBH_class::get_dr_metric(double State_vector[]) {

    if (false == this->eval_bitmask[1] || true == this->ignore_flag) {

        this->update_dr_metric(State_vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[1] = true;

        }

    }

    return this->s_dr_Metric;
}


int RBH_class::update_d2r_metric(double State_Vector[]) {

    Metric_type s_Metric    = this->get_metric(State_Vector);
    Metric_type s_dr_Metric = this->get_dr_metric(State_Vector);

    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho = sqrt(r2 + RBH_PARAM * RBH_PARAM);
    double rho3 = rho * rho * rho;
    double rho5 = rho * rho * rho * rho * rho;

    this->s_d2r_Metric.Metric[0][0] = -2 * MASS / rho3 + 6 * MASS * r2 / (rho5);
    this->s_d2r_Metric.Metric[0][3] = 0.0;
    this->s_d2r_Metric.Metric[1][1] = 1.0 / (this->s_Metric.Metric[0][0] * this->s_Metric.Metric[0][0]) * this->s_d2r_Metric.Metric[0][0] -
                       2.0 / (this->s_Metric.Metric[0][0] * this->s_Metric.Metric[0][0] * this->s_Metric.Metric[0][0]) * this->s_dr_Metric.Metric[0][0] * this->s_dr_Metric.Metric[0][0];
    this->s_d2r_Metric.Metric[2][2] = 2.0;
    this->s_d2r_Metric.Metric[3][3] = 2 * sin_theta * sin_theta;

    this->s_d2r_Metric.Lapse_function = -this->s_d2r_Metric.Metric[0][0];
    this->s_d2r_Metric.Shift_function = 0.0;

    return OK;

}

Metric_type RBH_class::get_d2r_metric(double State_vector[]) {

    if (false == this->eval_bitmask[2] || true == this->ignore_flag) {

        this->update_d2r_metric(State_vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[2] = true;

        }

    }

    return this->s_d2r_Metric;
}

int RBH_class::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

    double& r_obs = p_Initial_Conditions->init_Pos[e_r];
    double& theta_obs = p_Initial_Conditions->init_Pos[e_theta];

    p_Initial_Conditions->init_Three_Momentum[e_phi] = -J_data[photon] * sin(theta_obs);
    p_Initial_Conditions->init_Three_Momentum[e_theta] = p_theta_data[photon];

    double& J = p_Initial_Conditions->init_Three_Momentum[e_phi];

    double rho = sqrt(r_obs * r_obs + RBH_PARAM * RBH_PARAM);
    double rad_potential = 1 - (1 - 2 * MASS / rho) * J * J / (rho * rho);

    double (*metric)[4] = p_Initial_Conditions->init_metric;

    p_Initial_Conditions->init_Three_Momentum[e_r] = sqrt(rad_potential) * metric[1][1];

    return OK;
}

int RBH_class::get_EOM(double inter_State_vector[], double Derivatives[], int iteration) {

    double r = inter_State_vector[0 + iteration * e_State_Number];
    double rho = sqrt(r * r + RBH_PARAM * RBH_PARAM);

    double& J = inter_State_vector[e_p_phi];

    double sin1 = sin(inter_State_vector[1 + iteration * e_State_Number]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(inter_State_vector[1 + iteration * e_State_Number]);
    double cos2 = cos1 * cos1;

    Derivatives[e_r       + iteration * e_State_Number] = (1 - 2 * MASS / rho) * inter_State_vector[e_p_r + iteration * e_State_Number];
    Derivatives[e_theta   + iteration * e_State_Number] = 1.0 / (rho * rho) * inter_State_vector[e_p_theta + iteration * e_State_Number];
    Derivatives[e_phi     + iteration * e_State_Number] = J / (rho * rho * sin2);
    Derivatives[e_p_phi   + iteration * e_State_Number] = 0.0;
    Derivatives[e_p_theta + iteration * e_State_Number] = cos1 / (rho * rho * sin1 * sin2) * J * J;

    double r_term_1 = -MASS * r / (rho * rho * rho) * (1.0 / ((1 - 2 * MASS / rho) * (1 - 2 * MASS / rho)) + inter_State_vector[e_p_r + iteration * e_State_Number] * inter_State_vector[e_p_r + iteration * e_State_Number]);
    double r_term_2 = r / (rho * rho * rho * rho) * (inter_State_vector[e_p_theta + iteration * e_State_Number] * inter_State_vector[e_p_theta + iteration * e_State_Number] + J * J / sin2);

    Derivatives[e_p_r + iteration * e_State_Number] = r_term_1 + r_term_2;

    return 0;
}

bool RBH_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool scatter = State_vector[e_r] > 100 && Derivatives[e_r] < 0;

    double r_horizon = sqrt(4 * MASS * MASS - RBH_PARAM * RBH_PARAM);

    bool hit_horizon_RBH = State_vector[e_r] - r_horizon < 0.05;

    return scatter || hit_horizon_RBH;

}

void RBH_class::reset_eval_bitmask() {

    for (int index = 0; index <= 2; index++) {

        this->eval_bitmask[index] = false;

    }

}

void RBH_class::set_ignore_flag(bool flag) {

    this->ignore_flag = flag;

};

/***********************************************
|                                              |
| Derived Wormhole Class Functions Definitions |
|                                              |
***********************************************/

double* Wormhole_class::get_ISCO() {

    double M = MASS;

    static double r_ISCO[2]{};

    if (SPIN < 0.016) {

        r_ISCO[Inner] = 2 * M * (sqrt(4. / 9 * (6 * WH_REDSHIFT + 1)) * cosh(1. / 3 * acosh((1 + 9 * WH_REDSHIFT + 27. / 2 * WH_REDSHIFT * WH_REDSHIFT) / pow(6 * WH_REDSHIFT + 1, 3. / 2))) + 1. / 3);
        r_ISCO[Outer] = r_ISCO[Inner];
    }
    else {

        r_ISCO[Inner] = WH_R_THROAT;
        r_ISCO[Outer] = r_ISCO[Inner];

    }

    return r_ISCO;

}

double* Wormhole_class::get_Photon_Sphere() {

    double M = MASS;
    double a = SPIN;

    static double photon_orbit[2]{};

    photon_orbit[Inner] = M / 2 * (1 + sqrt(1 + 8 * WH_REDSHIFT));
    photon_orbit[Outer] = photon_orbit[Inner];

    return photon_orbit;

}

int Wormhole_class::update_metric(double State_Vector[]) {

    double& ell   = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r = sqrt(ell * ell + WH_R_THROAT * WH_R_THROAT);
    double r2 = r * r;
    double sin_theta = sin(theta);

    double exponent = -MASS / r - WH_REDSHIFT * MASS * MASS / r2;

    this->s_Metric.Lapse_function = exp(exponent);
    this->s_Metric.Shift_function = 2 * SPIN * MASS * MASS / r2 / r;

    this->s_Metric.Metric[0][0] = -this->s_Metric.Lapse_function * this->s_Metric.Lapse_function + 
                                   r2 * this->s_Metric.Shift_function * this->s_Metric.Shift_function * sin_theta * sin_theta;
    this->s_Metric.Metric[0][3] = -r2 * sin_theta * sin_theta * this->s_Metric.Shift_function;
    this->s_Metric.Metric[1][1] = 1. / (1 - WH_R_THROAT / r);
    this->s_Metric.Metric[2][2] = r2;
    this->s_Metric.Metric[3][3] = r2 * sin_theta * sin_theta;

    return OK;
}

Metric_type Wormhole_class::get_metric(double State_vector[]) {

    if (false == this->eval_bitmask[0] || true == this->ignore_flag) {

        this->update_metric(State_vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[0] = true;

        }

    }

    return this->s_Metric;
}

int Wormhole_class::update_dr_metric(double State_Vector[]) {

    Metric_type s_Metric = this->get_metric(State_Vector);

    double& ell   = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r = sqrt(ell * ell + WH_R_THROAT * WH_R_THROAT);
    double r2 = r * r;
    double sin_theta = sin(theta);

    this->s_dr_Metric.Lapse_function = this->s_Metric.Lapse_function * (1 / r2 + 2 * WH_REDSHIFT / (r2 * r));
    this->s_dr_Metric.Shift_function = -3 * this->s_Metric.Shift_function / r;

    double& N        = this->s_Metric.Lapse_function;
    double& dr_N     = this->s_dr_Metric.Lapse_function;
    double& omega    = this->s_Metric.Shift_function;
    double& dr_omega = this->s_dr_Metric.Shift_function;

    this->s_dr_Metric.Metric[0][0] = -2 * N * dr_N + 2 * r * omega * (omega + r * dr_omega) * sin_theta * sin_theta;
    this->s_dr_Metric.Metric[0][3] = -r * (2 * omega + r * dr_omega) * sin_theta * sin_theta;
    this->s_dr_Metric.Metric[1][1] = -1. / ((1 - WH_R_THROAT / r) * (1 - WH_R_THROAT / r)) * (WH_R_THROAT / r2);
    this->s_dr_Metric.Metric[2][2] = 2 * r;
    this->s_dr_Metric.Metric[3][3] = 2 * r * sin_theta * sin_theta;

    return OK;
}

Metric_type Wormhole_class::get_dr_metric(double State_vector[]) {

    if (false == this->eval_bitmask[1] || true == this->ignore_flag) {

        this->update_dr_metric(State_vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[1] = true;

        }

    }

    return this->s_dr_Metric;
}

int Wormhole_class::update_d2r_metric(double State_Vector[]) {

    Metric_type s_Metric    = this->get_metric(State_Vector);
    Metric_type s_dr_Metric = this->get_dr_metric(State_Vector);

    double& ell   = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r = sqrt(ell * ell + WH_R_THROAT * WH_R_THROAT);
    double r2 = r * r;
    double sin_theta = sin(theta);

    double& N        = this->s_Metric.Lapse_function;
    double& dr_N     = this->s_dr_Metric.Lapse_function;
    double& omega    = this->s_Metric.Shift_function;
    double& dr_omega = this->s_dr_Metric.Shift_function;

    this->s_d2r_Metric.Lapse_function = dr_N * (1 / r2 + 2 * WH_REDSHIFT / (r2 * r)) - N * (2. / (r2 * r) + 6 * WH_REDSHIFT / (r2 * r2));
    this->s_d2r_Metric.Shift_function = -3 * dr_omega / r + 3 * omega / r2;

    this->s_d2r_Metric.Metric[0][0] = -2 * dr_N * dr_N - 2 * N * this->s_d2r_Metric.Lapse_function + 2 * ((omega + r * dr_omega) * (omega + r * dr_omega) + 
                                      r * omega * (dr_omega + dr_omega + r * this->s_d2r_Metric.Shift_function)) * sin_theta * sin_theta;
     this->s_d2r_Metric.Metric[0][3] = -(2 * omega + r * dr_omega + r * (3 * dr_omega + r * this->s_d2r_Metric.Shift_function)) * sin_theta * sin_theta;
     this->s_d2r_Metric.Metric[1][1] = 2 / ((1 - WH_R_THROAT / r) * (1 - WH_R_THROAT / r)) * ((WH_R_THROAT / r2) * (WH_R_THROAT / r2) / (1 - WH_R_THROAT / r) + WH_R_THROAT / (r2 * r));
     this->s_d2r_Metric.Metric[2][2] = 2.0;
     this->s_d2r_Metric.Metric[3][3] = 2 * sin_theta * sin_theta;

    return OK;
}

Metric_type Wormhole_class::get_d2r_metric(double State_vector[]) {

    if (false == this->eval_bitmask[2] || true == this->ignore_flag) {

        this->update_d2r_metric(State_vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[2] = true;

        }

    }

    return this->s_d2r_Metric;
}

int Wormhole_class::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {


    double& r_obs = p_Initial_Conditions->init_Pos[e_r];
    double& theta_obs = p_Initial_Conditions->init_Pos[e_theta];

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

int Wormhole_class::get_EOM(double inter_State_vector[], double Derivatives[], int iteration) {

    double sqrt_r2 = sqrt(inter_State_vector[0 + iteration * e_State_Number] * inter_State_vector[0 + iteration * e_State_Number] + WH_R_THROAT * WH_R_THROAT);
    double d_ell_r = inter_State_vector[0 + iteration * e_State_Number] / sqrt_r2;

    double& J = inter_State_vector[e_p_phi];

    double omega = 2 * SPIN / (sqrt_r2 * sqrt_r2 * sqrt_r2);
    double d_ell_omega = -3 * omega / sqrt_r2 * d_ell_r;

    double exponent = -1 / sqrt_r2 - WH_REDSHIFT / (sqrt_r2 * sqrt_r2);
    double N = this->get_metric(inter_State_vector).Lapse_function;
    double d_ell_N = N * (1 / (sqrt_r2 * sqrt_r2) + 2 * WH_REDSHIFT / (sqrt_r2 * sqrt_r2 * sqrt_r2)) * d_ell_r;

    double N2 = N * N;

    double sin1 = sin(inter_State_vector[e_theta + iteration * e_State_Number]);
    double sin2 = sin1 * sin1;

    Derivatives[e_r       + iteration * e_State_Number] = 1.0 / (1 + WH_R_THROAT / sqrt_r2) * inter_State_vector[e_p_r + iteration * e_State_Number];
    Derivatives[e_theta   + iteration * e_State_Number] = 1.0 / (sqrt_r2 * sqrt_r2) * inter_State_vector[e_p_theta + iteration * e_State_Number];
    Derivatives[e_phi     + iteration * e_State_Number] = J / (sqrt_r2 * sqrt_r2 * sin2) + omega * (1 - omega * J) / N2;
    Derivatives[e_p_phi   + iteration * e_State_Number] = 0.0;
    Derivatives[e_p_theta + iteration * e_State_Number] = (cos(inter_State_vector[1 + iteration * e_State_Number]) / sin1) / (sqrt_r2 * sqrt_r2) * J * J / sin2;

    double term_1 = -1.0 / ((1 + WH_R_THROAT / sqrt_r2) * (1 + WH_R_THROAT / sqrt_r2)) * WH_R_THROAT * inter_State_vector[e_r + iteration * e_State_Number] / (sqrt_r2 * sqrt_r2 * sqrt_r2) * inter_State_vector[e_p_r + iteration * e_State_Number] * inter_State_vector[e_p_r + iteration * e_State_Number] / 2;
    double term_2 = 1.0 / (sqrt_r2 * sqrt_r2 * sqrt_r2) * (inter_State_vector[e_p_theta + iteration * e_State_Number] * inter_State_vector[e_p_theta + iteration * e_State_Number] + J * J / sin2) * d_ell_r;
    double term_3 = -(1.0 / (N2 * N) * d_ell_N * ((1 - omega * J) * (1 - omega * J)) - 1.0 / N2 * (-d_ell_omega * (1 - omega * J) * J));

    Derivatives[e_p_r + iteration * e_State_Number] = term_1 + term_2 + term_3;

    return OK;
}

bool Wormhole_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool scatter            = State_vector[e_r] >  sqrt(100 * 100 + WH_R_THROAT * WH_R_THROAT) && Derivatives[0] < 0;
    bool scatter_other_side = State_vector[e_r] < -sqrt(100 * 100 + WH_R_THROAT * WH_R_THROAT);
    bool stop_at_throat     = State_vector[e_r] < 1e-5;

    if (STOP_AT_THROAT) {

        return scatter || stop_at_throat;
    }
    else {

        return scatter || scatter_other_side;

    }
};

void Wormhole_class::reset_eval_bitmask() {

    for (int index = 0; index <= 2; index++) {

        this->eval_bitmask[index] = false;

    }

}

void Wormhole_class::set_ignore_flag(bool flag) {

    this->ignore_flag = flag;

};


/******************************************************************************
|                                                                             |
| Derived Janis-Newman-Winicour Naked Singularity Class Functions Definitions |
|                                                                             |
******************************************************************************/

double* JNW_class::get_ISCO() {

    static double r_ISCO[2]{};

    if (JNW_GAMMA > 1.0 / sqrt(5)) {

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

    if (false == this->eval_bitmask[0] || true == this->ignore_flag) {

        this->update_metric(State_Vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[0] = true;

        }
    }

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

    if (false == this->eval_bitmask[1] || true == this->ignore_flag) {

        this->update_dr_metric(State_Vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[1] = true;

        }

    }

    return this->s_dr_Metric;

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

    if (false == this->eval_bitmask[2] || true == this->ignore_flag) {

        this->update_d2r_metric(State_Vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[2] = true;

        }

    }

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

int JNW_class::get_EOM(double inter_State_vector[], double Derivatives[], int iteration)
{

    double r = inter_State_vector[0 + iteration * e_State_Number];

    double& J = inter_State_vector[e_p_phi];

    double sin1 = sin(inter_State_vector[1 + iteration * e_State_Number]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(inter_State_vector[1 + iteration * e_State_Number]);
    double cos2 = cos1 * cos1;

    double pow_gamma = pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA);
    double pow_gamma_minus_1 = pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA - 1);

    Derivatives[e_r       + iteration * e_State_Number] = pow_gamma * inter_State_vector[e_p_r + iteration * e_State_Number];
    Derivatives[e_theta   + iteration * e_State_Number] = pow_gamma_minus_1 / (r * r) * inter_State_vector[e_p_theta + iteration * e_State_Number];
    Derivatives[e_phi     + iteration * e_State_Number] = pow_gamma_minus_1 / (r * r * sin2) * J;
    Derivatives[e_p_phi   + iteration * e_State_Number] = 0.0;
    Derivatives[e_p_theta + iteration * e_State_Number] = pow_gamma_minus_1 * cos1 / (r * r * sin1 * sin2) * J * J;

    double r_term_1 = -JNW_GAMMA * JNW_R_SINGULARITY / 2 / r / r * pow_gamma_minus_1 * (1.0 / pow_gamma / pow_gamma
                    + inter_State_vector[e_p_r + iteration * e_State_Number] * inter_State_vector[e_p_r + iteration * e_State_Number]);
    double r_term_2 = 1.0 / r / r / r * pow_gamma_minus_1 * (1 - JNW_R_SINGULARITY / 2 / r * (JNW_GAMMA - 1) / (1 - JNW_R_SINGULARITY / r))
                    * (inter_State_vector[e_p_theta + iteration * e_State_Number] * inter_State_vector[e_p_theta + iteration * e_State_Number] + J * J / sin2);

    Derivatives[e_p_r + iteration * e_State_Number] = r_term_1 + r_term_2;

    return OK;

}

bool JNW_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool hit_singularity = State_vector[e_r] - JNW_R_SINGULARITY < 1e-9;

    bool scatter = State_vector[e_r] > 100 && Derivatives[e_r] < 0;

    return scatter ;

};

void JNW_class::reset_eval_bitmask() {

    for (int index = 0; index <= 2; index++) {

        this->eval_bitmask[index] = false;

    }

}

void JNW_class::set_ignore_flag(bool flag) {

    this->ignore_flag = flag;

};

/**********************************************************************************
|                                                                                 |
| Derived Gauss-Bonnet Naked Singularity / Black Hole Class Functions Definitions |
|                                                                                 |
**********************************************************************************/

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

    double Gamma2  = GAUSS_BONNET_GAMMA * GAUSS_BONNET_GAMMA;
    double Gamma4  = Gamma2 * Gamma2;
    double Gamma8  = Gamma4 * Gamma4;
    double Gamma10 = Gamma8 * Gamma2;

    r_ISCO[Outer] = fit_coeffs[0]  + 
                    fit_coeffs[1]  * GAUSS_BONNET_GAMMA +
                    fit_coeffs[2]  * Gamma2 +
                    fit_coeffs[3]  * Gamma2 * GAUSS_BONNET_GAMMA +
                    fit_coeffs[4]  * Gamma4 +
                    fit_coeffs[5]  * Gamma4 * GAUSS_BONNET_GAMMA +
                    fit_coeffs[6]  * Gamma4 * Gamma2 + 
                    fit_coeffs[7]  * Gamma8 / GAUSS_BONNET_GAMMA + 
                    fit_coeffs[8]  * Gamma8 + 
                    fit_coeffs[9]  * Gamma8 * GAUSS_BONNET_GAMMA + 
                    fit_coeffs[10] * Gamma10;

    r_ISCO[Inner] = pow(GAUSS_BONNET_GAMMA, 1.0 / 3);

    return r_ISCO;

};

double* Gauss_Bonnet_class::get_Photon_Sphere() {

    /* This expression is the root of a cubic equation */

    double q =  8 * MASS * GAUSS_BONNET_GAMMA;
    double p = -9 * MASS * MASS;

    static double photon_orbits[2]{};

    photon_orbits[Outer] = 2 * sqrt(-p / 3) * cos(1. / 3 * acos(3. / 2 * q / p * sqrt(-3. / p)));
    photon_orbits[Inner] = 2 * sqrt(-p / 3) * cos(1. / 3 * acos(3. / 2 * q / p * sqrt(-3. / p)) + 2. * M_PI / 3);

    return photon_orbits;

};

int Gauss_Bonnet_class::update_metric(double State_Vector[]) {

    double  M     = MASS;
    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double f = 1. + r2 / GAUSS_BONNET_GAMMA / 2. * (1. - sqrt(1. + 8. * GAUSS_BONNET_GAMMA * M / r2 / r));

    this->s_Metric.Metric[0][0] = -f;
    this->s_Metric.Metric[1][1] = 1. / f;
    this->s_Metric.Metric[2][2] = r2;
    this->s_Metric.Metric[3][3] = r2 * sin_theta * sin_theta;
    this->s_Metric.Metric[0][3] = 0.;
    this->s_Metric.Metric[3][0] = 0.;

    this->s_Metric.Lapse_function = -this->s_Metric.Metric[0][0];
    this->s_Metric.Shift_function = 0.;

    return OK;

}

Metric_type Gauss_Bonnet_class::get_metric(double State_Vector[]) {

    if (false == this->eval_bitmask[0] || true == this->ignore_flag) {

        this->update_metric(State_Vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[0] = true;

        }

    }

    return this->s_Metric;

}

int Gauss_Bonnet_class::update_dr_metric(double State_Vector[]) {

    double M      = MASS;
    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double f    = 1. + r2 / GAUSS_BONNET_GAMMA / 2. * (1. - sqrt(1. + 8. * GAUSS_BONNET_GAMMA * M / r2 / r));
    double dr_f = 2. / r * (f - 1.) + 6. * M / sqrt(r2 * r2 + 8. * GAUSS_BONNET_GAMMA * M * r);

    this->s_dr_Metric.Metric[0][0] = -dr_f;
    this->s_dr_Metric.Metric[1][1] = -1. / f / f * dr_f;
    this->s_dr_Metric.Metric[2][2] = 2. * r;
    this->s_dr_Metric.Metric[3][3] = 2. * r * sin_theta * sin_theta;
    this->s_dr_Metric.Metric[0][3] = 0.;
    this->s_dr_Metric.Metric[3][0] = 0.;

    this->s_dr_Metric.Lapse_function = -this->s_dr_Metric.Metric[0][0];
    this->s_dr_Metric.Shift_function = 0.;

    return OK;

}

Metric_type Gauss_Bonnet_class::get_dr_metric(double State_Vector[]) {

    if (false == this->eval_bitmask[1] || true == this->ignore_flag) {

        this->update_dr_metric(State_Vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[1] = true;

        }

    }

    return this->s_dr_Metric;

}

int Gauss_Bonnet_class::update_d2r_metric(double State_Vector[]) {

    double  M     = MASS;
    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double root = sqrt(r2 * r2 + 8 * GAUSS_BONNET_GAMMA * M * r);

    double f     =  1 + r2 / GAUSS_BONNET_GAMMA / 2. * (1 - sqrt(1. + 8. * GAUSS_BONNET_GAMMA * M / r2 / r));
    double dr_f  =  2. / r * (f - 1.) + 6 * M / root;
    double d2r_f = -2. / r2 * (f - 1.) + 2. / r * dr_f - 12. * M / root / root / root * (r2 * r + 2. * GAUSS_BONNET_GAMMA * M);

    this->s_d2r_Metric.Metric[0][0] = -d2r_f;
    this->s_d2r_Metric.Metric[1][1] = 2. / f / f / f * dr_f - 1. / f / f * d2r_f;
    this->s_d2r_Metric.Metric[2][2] = 2.;
    this->s_d2r_Metric.Metric[3][3] = 2. * sin_theta * sin_theta;
    this->s_d2r_Metric.Metric[0][3] = 0.;
    this->s_d2r_Metric.Metric[3][0] = 0.;

    this->s_d2r_Metric.Lapse_function = -this->s_d2r_Metric.Metric[0][0];
    this->s_d2r_Metric.Shift_function = 0.;

    return OK;

}

Metric_type Gauss_Bonnet_class::get_d2r_metric(double State_Vector[]) {

    if (false == this->eval_bitmask[2] || true == this->ignore_flag) {

        this->update_d2r_metric(State_Vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[2] = true;

        }

    }

    return this->s_d2r_Metric;

}

int Gauss_Bonnet_class::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

    double& r_obs     = p_Initial_Conditions->init_Pos[e_r];
    double& theta_obs = p_Initial_Conditions->init_Pos[e_theta];

    p_Initial_Conditions->init_Three_Momentum[e_phi] = -J_data[photon] * sin(theta_obs);
    p_Initial_Conditions->init_Three_Momentum[e_theta] = p_theta_data[photon];

    double& J = p_Initial_Conditions->init_Three_Momentum[e_phi];
    double  f = 1. + r_obs * r_obs / 2. / GAUSS_BONNET_GAMMA * (1. - sqrt(1. + 8. * GAUSS_BONNET_GAMMA * MASS / r_obs / r_obs / r_obs));

    double rad_potential = 1. - f * J * J / (r_obs * r_obs);

    double(*metric)[4] = p_Initial_Conditions->init_metric;

    p_Initial_Conditions->init_Three_Momentum[e_r] = sqrt(rad_potential) * metric[1][1];

    return OK;

}

int Gauss_Bonnet_class::get_EOM(double inter_State_vector[], double Derivatives[], int iteration){

    double& r = inter_State_vector[e_r + iteration * e_State_Number];

    double& J = inter_State_vector[e_p_phi];

    double sin1 = sin(inter_State_vector[e_theta + iteration * e_State_Number]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(inter_State_vector[e_theta + iteration * e_State_Number]);
    double cos2 = cos1 * cos1;

    double root = sqrt(1. + 8. * GAUSS_BONNET_GAMMA * MASS / r / r / r);

    double f    = 1. + r * r / GAUSS_BONNET_GAMMA / 2. * (1. - root);
    double dr_f = 2. / r * (f - 1.) + 6. * MASS / root / r / r;

    Derivatives[e_r       + iteration * e_State_Number] = f * inter_State_vector[e_p_r + iteration * e_State_Number];
    Derivatives[e_theta   + iteration * e_State_Number] = 1. / (r * r) * inter_State_vector[e_p_theta + iteration * e_State_Number];
    Derivatives[e_phi     + iteration * e_State_Number] = J / (r * r * sin2);
    Derivatives[e_p_phi   + iteration * e_State_Number] = 0.0;
    Derivatives[e_p_theta + iteration * e_State_Number] = cos1 / (r * r * sin1 * sin2) * J * J;

    double r_term_1 = -1. / 2 * (1.0 / f / f + inter_State_vector[e_p_r] * inter_State_vector[e_p_r]) * dr_f;
    double r_term_2 = 1.0 / r / r / r * (inter_State_vector[e_p_theta + iteration * e_State_Number] * inter_State_vector[e_p_theta + iteration * e_State_Number] + J * J / sin2);

    Derivatives[e_p_r + iteration * e_State_Number] = r_term_1 + r_term_2;

    return OK;

}

bool Gauss_Bonnet_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool scatter = State_vector[e_r] > 100 && Derivatives[e_r] < 0;

    bool too_high_order = State_vector[e_phi] * State_vector[e_phi] > 5 * M_PI * 5 * M_PI;

    double r_horizon{};
    bool hit_horizon = false;

    if (GAUSS_BONNET_GAMMA <= 1) {

        r_horizon = MASS + sqrt(MASS * MASS - GAUSS_BONNET_GAMMA * GAUSS_BONNET_GAMMA);
        hit_horizon = State_vector[e_r] - r_horizon < 1e-5;

    }

    return scatter || too_high_order || hit_horizon;

};

void Gauss_Bonnet_class::reset_eval_bitmask() {

    for (int index = 0; index <= 2; index++) {

        this->eval_bitmask[index] = false;

    }

}

void Gauss_Bonnet_class::set_ignore_flag(bool flag) {

    this->ignore_flag = flag;

};


/***********************************************************************
|                                                                      |
| Derived Black Hole with Dark Matter Halo Class Functions Definitions |
|                                                                      |
***********************************************************************/

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

    if (M_HALO > 1e2 + 1) {

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

    double Compactness2  = COMPACTNESS  * COMPACTNESS;
    double Compactness4  = Compactness2 * Compactness2;
    double Compactness8  = Compactness4 * Compactness4;
    double Compactness10 = Compactness8 * Compactness2;

    r_ISCO[Outer] = fit_coeffs[0]  +
                    fit_coeffs[1]  * COMPACTNESS +
                    fit_coeffs[2]  * Compactness2 +
                    fit_coeffs[3]  * Compactness2 * COMPACTNESS +
                    fit_coeffs[4]  * Compactness4 +
                    fit_coeffs[5]  * Compactness4 * COMPACTNESS +
                    fit_coeffs[6]  * Compactness4 * Compactness2 +
                    fit_coeffs[7]  * Compactness8 / COMPACTNESS +
                    fit_coeffs[8]  * Compactness8 +
                    fit_coeffs[9]  * Compactness8 * COMPACTNESS +
                    fit_coeffs[10] * Compactness10;

    r_ISCO[Inner] = r_ISCO[Outer];

    return r_ISCO;

};

int Black_Hole_w_Dark_Matter_Halo_class::update_metric(double State_Vector[]) {

    double  M     = MASS;
    double& r     = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double ksi = 2 * A_0 - M_HALO + 4 * M;
    double Y   = sqrt(M_HALO / ksi) * (2 * atan((r + A_0 + M_HALO) / sqrt(M_HALO * ksi)) - M_PI);
    double f   = (1 - 2 * M / r) * exp(Y);
    double m   = M + M_HALO * r2 / (A_0 + r) / (A_0 + r) * (1 - 2 * M / r) * (1 - 2 * M / r);

    this->s_Metric.Metric[0][0] = -f;
    this->s_Metric.Metric[1][1] = 1. / (1 - 2 * m / r);
    this->s_Metric.Metric[2][2] = r2;
    this->s_Metric.Metric[3][3] = r2 * sin_theta * sin_theta;
    this->s_Metric.Metric[0][3] = 0.;
    this->s_Metric.Metric[3][0] = 0.;

    this->s_Metric.Lapse_function = -this->s_Metric.Metric[0][0];
    this->s_Metric.Shift_function = 0.;

    return OK;
}

Metric_type Black_Hole_w_Dark_Matter_Halo_class::get_metric(double State_Vector[]) {

    if (false == this->eval_bitmask[0] || true == this->ignore_flag) {

        this->update_metric(State_Vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[0] = true;

        }

    }

    return this->s_Metric;

}

int Black_Hole_w_Dark_Matter_Halo_class::update_dr_metric(double State_Vector[]) {

    double M = MASS;
    double& r = State_Vector[e_r];
    double& theta = State_Vector[e_theta];

    double r2 = r * r;
    double sin_theta = sin(theta);

    double ksi  = 2 * A_0 - M_HALO + 4 * M;
    double Y    = sqrt(M_HALO / ksi) * (2 * atan((r + A_0 + M_HALO) / sqrt(M_HALO * ksi)) - M_PI);
    double dr_Y = 2 / ksi / (1 + (r + A_0 + M_HALO) * (r + A_0 + M_HALO) / M / ksi);

    double exp_y = exp(Y);

    double f    = (1 - 2 * M / r) * exp_y;
    double dr_f = 2 * M / r2 * exp_y + f * dr_Y;
    double m    = M + M_HALO * r2 / (A_0 + r) / (A_0 + r) * (1 - 2 * M / r) * (1 - 2 * M / r);
    double dr_m = 2 * (1 - 2 * M / r) * ((1 - 2 * M / r) * (1 - r / (r + A_0)) * r + 2 * M) * M_HALO / (r + A_0) / (r + A_0);

    this->s_Metric.Metric[0][0] = -dr_f;
    this->s_Metric.Metric[1][1] = -1. / (1 - 2 * m / r) / (1 - 2 * m / r) * (2 * m / r2 - 2 / r * dr_m);
    this->s_Metric.Metric[2][2] = r2;
    this->s_Metric.Metric[3][3] = r2 * sin_theta * sin_theta;
    this->s_Metric.Metric[0][3] = 0.;
    this->s_Metric.Metric[3][0] = 0.;

    this->s_Metric.Lapse_function = -this->s_Metric.Metric[0][0];
    this->s_Metric.Shift_function = 0.;

    return OK;
}

Metric_type Black_Hole_w_Dark_Matter_Halo_class::get_dr_metric(double State_Vector[]) {

    if (false == this->eval_bitmask[1] || true == this->ignore_flag) {

        this->update_dr_metric(State_Vector);

        if (false == this->ignore_flag) {

            this->eval_bitmask[1] = true;

        }

    }

    return this->s_dr_Metric;

}

int Black_Hole_w_Dark_Matter_Halo_class::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

    double& r_obs = p_Initial_Conditions->init_Pos[e_r];
    double& theta_obs = p_Initial_Conditions->init_Pos[e_theta];

    p_Initial_Conditions->init_Three_Momentum[e_phi] = -J_data[photon] * sin(theta_obs);
    p_Initial_Conditions->init_Three_Momentum[e_theta] = p_theta_data[photon];

    double& J = p_Initial_Conditions->init_Three_Momentum[e_phi];

    double ksi = 2 * A_0 - M_HALO + 4 * MASS;
    double Y = sqrt(M_HALO / ksi) * (2 * atan((r_obs + A_0 + M_HALO) / sqrt(M_HALO * ksi)) - M_PI);
    double dr_Y = 2 / ksi / (1 + (r_obs + A_0 + M_HALO) * (r_obs + A_0 + M_HALO) / MASS / ksi);

    double exp_y = exp(Y);

    double f = (1 - 2 * MASS / r_obs) * exp_y;

    double rad_potential = 1. - f * J * J / (r_obs * r_obs);

    double(*metric)[4] = p_Initial_Conditions->init_metric;

    p_Initial_Conditions->init_Three_Momentum[e_r] = sqrt(rad_potential) * metric[1][1];

    return OK;

}

int Black_Hole_w_Dark_Matter_Halo_class::get_EOM(double inter_State_vector[], double Derivatives[], int iteration) {

    double& r = inter_State_vector[e_r + iteration * e_State_Number];

    double& J = inter_State_vector[e_p_phi];

    double sin1 = sin(inter_State_vector[e_theta + iteration * e_State_Number]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(inter_State_vector[e_theta + iteration * e_State_Number]);
    double cos2 = cos1 * cos1;

    double M  = MASS;
    double r2 = r * r;

    double ksi   = 2 * A_0 - M_HALO + 4 * M;
    double Y     = sqrt(M_HALO / ksi) * (2 * atan((r + A_0 + M_HALO) / sqrt(M_HALO * ksi)) - M_PI);
    double dr_Y = 2 / ksi / (1 + (r + A_0 + M_HALO) * (r + A_0 + M_HALO) / M / ksi);

    double exp_y = exp(Y);

    double f    = (1 - 2 * M / r) * exp_y;
    double dr_f = 2 * M / r2 * exp_y + f * dr_Y;
    double m    = M + M_HALO * r2 / (A_0 + r) / (A_0 + r) * (1 - 2 * M / r) * (1 - 2 * M / r);
    double dr_m = 2 * (1 - 2 * M / r) * ((1 - 2 * M / r) * (1 - r / (r + A_0)) * r + 2 * M) * M_HALO / (r + A_0) / (r + A_0);

    Derivatives[e_r       + iteration * e_State_Number] = (1 - 2 * m / r) * inter_State_vector[e_p_r + iteration * e_State_Number];
    Derivatives[e_theta   + iteration * e_State_Number] = 1. / (r * r) * inter_State_vector[e_p_theta + iteration * e_State_Number];
    Derivatives[e_phi     + iteration * e_State_Number] = J / (r * r * sin2);
    Derivatives[e_p_phi   + iteration * e_State_Number] = 0.0;
    Derivatives[e_p_theta + iteration * e_State_Number] = cos1 / (r * r * sin1 * sin2) * J * J;

    double r_term_1 = -1. / 2 / f / f * dr_f + (dr_m / r - m / r2) * inter_State_vector[e_p_r + iteration * e_State_Number] * inter_State_vector[e_p_r + iteration * e_State_Number];
    double r_term_2 = 1.0 / r / r / r * (inter_State_vector[e_p_theta + iteration * e_State_Number] * inter_State_vector[e_p_theta + iteration * e_State_Number] + J * J / sin2);

    Derivatives[e_p_r + iteration * e_State_Number] = r_term_1 + r_term_2;

    return OK;

}

bool Black_Hole_w_Dark_Matter_Halo_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool scatter     = State_vector[e_r] > 100 && Derivatives[e_r] < 0;
    bool hit_horizon = State_vector[e_r] - 2 * MASS < 1e-5;

    return scatter || hit_horizon;

};

void Black_Hole_w_Dark_Matter_Halo_class::reset_eval_bitmask() {

    for (int index = 0; index <= 2; index++) {

        this->eval_bitmask[index] = false;

    }

}

void Black_Hole_w_Dark_Matter_Halo_class::set_ignore_flag(bool flag) {

    this->ignore_flag = flag;

};