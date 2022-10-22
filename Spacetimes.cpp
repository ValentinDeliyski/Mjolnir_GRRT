#pragma once

#include <cmath>
#include <iostream>
#include <vector>

#include "Constants.h"
#include "Enumerations.h"
#include "Spacetimes.h"
#include "General_functions.h"


/***************************************
|                                      |
| Observer Class Functions Definitions |
|                                      |
***************************************/

extern e_Spacetimes e_metric;

tag_observer::tag_observer(double r, double theta, double phi) {

    r_obs     = r;
    theta_obs = theta;
    phi_obs   = phi;

}

double tag_observer::get_r_obs()     { return r_obs; };
double tag_observer::get_theta_obs() { return theta_obs; };
double tag_observer::get_phi_obs()   { return phi_obs; };

int tag_observer::get_obs_velocity(double Obs_velocity[4], std::vector<c_Spacetime_Base*> Spacetimes) {


    double metric_obs[4][4], N_obs, omega_obs;

    Spacetimes[e_metric]->get_metric(metric_obs, &N_obs, &omega_obs, r_obs, theta_obs);

    /*
    Obs_velocity is given in contravatiant components 
    */

    Obs_velocity[0] = 1.0 / N_obs;
    Obs_velocity[1] = 0;
    Obs_velocity[2] = 0;
    Obs_velocity[3] = omega_obs / N_obs;

    return OK;

}

/*******************************************
|                                          |
| Derived Kerr Class Functions Definitions |
|                                          |
*******************************************/

double derived_Kerr_class::get_ISCO(Orbit_Orientation Orientation) {

    double Z_1 = 1 + pow(1 - SPIN * SPIN, 1. / 3) * (pow(1 + SPIN, 1. / 3) + pow(1 - SPIN, 1. / 3));
    double Z_2 = sqrt(3 * SPIN * SPIN + Z_1 * Z_1);

    switch (Orientation) {

        case Prograde:

            return MASS * (3 + Z_2 - sqrt((3 - Z_1) * (3 + Z_1 + 2 * Z_2)));

        case Retrograde:

            return MASS * (3 + Z_2 + sqrt((3 - Z_1) * (3 + Z_1 + 2 * Z_2)));

        default:

            std::cout << "Wrong Orbit Orientation! - Must be 'Prograde' or 'Retrograde'!" << '\n';

            return ERROR;

    }

}

double derived_Kerr_class::get_Photon_Sphere(Orbit_Orientation Orientation) {

    switch (Orientation) {

    case Prograde:

        return 2 * MASS * (1 + cos(2.0 / 3 * acos(SPIN)));

    case Retrograde:

        return 2 * MASS * (1 + cos(2.0 / 3 * acos(-SPIN)));

    default:

        std::cout << "Wrong Orbit Orientation! - Must be 'Prograde' or 'Retrograde'!" << '\n';

        return ERROR;

    }

}

int derived_Kerr_class::get_metric(double metric[4][4], double* N_metric, double* omega_metric, double r, double theta) {

    double M = MASS;
    double a = SPIN;

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r2 + a * a * cos_theta * cos_theta;
    double delta = r2 - 2 * M * r + a * a;

    metric[0][0] = -(1 - 2 * M * r / rho2);
    metric[0][3] = -2 * M * r * a * sin_theta * sin_theta / rho2;
    metric[3][0] = metric[0][3];
    metric[1][1] = rho2 / delta;
    metric[2][2] = rho2;
    metric[3][3] = (r2 + a * a + 2 * M * r * a * a / rho2 * sin_theta * sin_theta) * sin_theta * sin_theta;

    double sigma2 = metric[3][3] * rho2 / sin_theta / sin_theta;

    *N_metric = sqrt(rho2 * delta / sigma2);
    *omega_metric = 2 * a * r / sigma2;

    return OK;
}

int derived_Kerr_class::get_dr_metric(double dr_metric[4][4], double* dr_N, double* dr_omega, double r, double theta) {

    double M = MASS;
    double a = SPIN;

    double metric[4][4], N, omega;

    derived_Kerr_class::get_metric(metric, &N, &omega, r, theta);

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r2 + a * a * cos_theta * cos_theta;
    double delta = r2 - 2 * M * r + a * a;

    dr_metric[0][0] = -2 * M / rho2 * (2 * r2 / rho2 - 1);
    dr_metric[0][3] = 2 * M * a * sin_theta * sin_theta / rho2 * (2 * r2 / rho2 - 1);
    dr_metric[1][1] = 2 * r / delta * (1 - rho2 / delta * (1 - M / r));
    dr_metric[2][2] = 2 * r;
    dr_metric[3][3] = 2 * (r - M * a * a / rho2 * (2 * r2 / rho2 - 1) * sin_theta * sin_theta) * sin_theta * sin_theta;

    double sigma2 = rho2 * metric[3][3] / sin_theta / sin_theta;
    double dr_sigma2 = 2 * r * metric[3][3] + rho2 * dr_metric[3][3];

    *dr_N = N * (r / rho2 + (r - M) / delta - dr_sigma2 / 2 / sigma2);
    *dr_omega = omega / r * (1 - r * dr_sigma2 / sigma2);

    return OK;

}

int derived_Kerr_class::get_d2r_metric(double d2r_metric[4][4], double* d2r_N, double* d2r_omega, double r, double theta) {

    double M = MASS;
    double a = SPIN;

    double metric[4][4], N, omega;

    derived_Kerr_class::get_metric(metric, &N, &omega, r, theta);

    double dr_metric[4][4], dr_N, dr_omega;

    derived_Kerr_class::get_dr_metric(dr_metric, &dr_N, &dr_omega, r, theta);

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r2 + a * a * cos_theta * cos_theta;
    double delta = r2 - 2 * M * r + a * a;

    d2r_metric[0][0] = 4 * M * r / rho2 / rho2 * (4 * r2 / rho2 - 3);
    d2r_metric[0][3] = -4 * M * a * r * sin_theta * sin_theta / rho2 / rho2 * (4 * r2 / rho2 - 3);
    d2r_metric[1][1] = 2 / delta * (1 - 4 * (r2 - r * M) / delta + rho2 / delta * (2 * (r - M) * (r - M) / delta - 1));
    d2r_metric[2][2] = 2.0;
    d2r_metric[3][3] = 2 * (1 + 2 * M * a * a * r / rho2 / rho2 * (4 * r2 / rho2 - 3) * sin_theta * sin_theta) * sin_theta * sin_theta;

    double sigma2 = rho2 * metric[3][3] / sin_theta / sin_theta;
    double dr_sigma2 = (2 * r * metric[3][3] + rho2 * dr_metric[3][3]) / sin_theta / sin_theta;
    double d2r_sigma2 = (2 * metric[3][3] + 4 * r * dr_metric[3][3] + rho2 * d2r_metric[3][3]) / sin_theta / sin_theta;

    *d2r_N = dr_N * dr_N / N + N / rho2 * (1 - 2 * r2 / rho2 + rho2 / delta * (1 - (r - M) * (r - M) / delta) - rho2 / sigma2 / 2 * (d2r_sigma2 - dr_sigma2 * dr_sigma2 / sigma2));
    *d2r_omega = -omega / r2 * (1 - r * dr_omega / omega + r * dr_sigma2 / sigma2) * (1 - r * dr_sigma2 / sigma2);

    return OK;

}

int derived_Kerr_class::get_initial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                                                         int photon, double r_obs, double theta_obs, double metric[4][4], double N, double omega) {

    double M = MASS;
    double a = SPIN;

    *J = -J_data[photon] * sin(theta_obs);
    *p_theta = p_theta_data[photon];

    double delta = pow(r_obs, 2) + pow(a, 2) - 2 * M * r_obs;
    double K = pow(*p_theta, 2) + pow(cos(theta_obs), 2) * (pow(*J / sin(theta_obs), 2) - pow(a, 2));

    double rad_potential = pow(r_obs * r_obs + a * a - a * *J, 2) - delta * (pow(*J - a, 2) + K);

    *p_r = sqrt(rad_potential) / delta;

    return OK;
}

int derived_Kerr_class::get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration) {

    double& r = inter_State_vector[e_r + iteration * e_State_Number];
    double r2 = r * r;

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
    Derivatives[e_phi_FD + iteration * e_State_Number] = -1.0 / rho2; /* Mino Time */

    double theta_term_1 = -(delta * p_r * p_r + p_theta * p_theta) * SPIN * SPIN * cos1 * sin1 / (rho2 * rho2);
    double theta_term_2 = F * SPIN * SPIN * cos1 * sin1 / (delta * rho2 * rho2) + (J * J * cos1 / (sin2 * sin1) - SPIN * SPIN * cos1 * sin1) / rho2;

    Derivatives[e_p_theta + iteration * e_State_Number] = theta_term_1 + theta_term_2;

    double r_term_1 = p_r * p_r / (rho2) * (MASS - r * (1 - delta / rho2)) + p_theta * p_theta * r / (rho2 * rho2);
    double r_term_2 = (2 * P * r - (r - MASS) * ((J - SPIN) * (J - SPIN) + cos2 * (J * J / (sin2) - SPIN * SPIN))) / (delta * rho2)
                    - F * (rho2 * (r - MASS) + r * delta) / (delta * delta * rho2 * rho2);

    Derivatives[e_p_r + iteration * e_State_Number] = r_term_1 + r_term_2;

    return OK;

}

bool derived_Kerr_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool scatter = State_vector[e_r] > 300 && Derivatives[e_r] < 0;

    double r_horizon = MASS * (1 + sqrt(1 - SPIN * SPIN));

    bool hit_horizon = State_vector[e_r] - r_horizon < 0.05;


    return scatter || hit_horizon;

};

/*********************************************************
|                                                        |
| Derived Regular Black Hole Class Functions Definitions |
|                                                        |
*********************************************************/

double derived_RBH_class::get_ISCO(Orbit_Orientation Orientation) {

    double M = MASS;

    switch (Orientation) {

    case Prograde:

        return sqrt(36 * MASS * MASS - RBH_PARAM * RBH_PARAM);

    case Retrograde:

        return sqrt(36 * MASS * MASS - RBH_PARAM * RBH_PARAM);

    default:

        std::cout << "Wrong Orbit Orientation! - Must be 'Prograde' or 'Retrograde'!" << '\n';

        return ERROR;

    }

}

double derived_RBH_class::get_Photon_Sphere(Orbit_Orientation Orientation) {

    double M = MASS;

    switch (Orientation) {

    case Prograde:

        return sqrt(9 * MASS * MASS - RBH_PARAM * RBH_PARAM);

    case Retrograde:

        return sqrt(9 * MASS * MASS - RBH_PARAM * RBH_PARAM);

    default:

        std::cout << "Wrong Orbit Orientation! - Must be 'Prograde' or 'Retrograde'!" << '\n';

        return ERROR;

    }

}

int derived_RBH_class::get_metric(double metric[4][4], double* N_metric, double* omega_metric, double r, double theta) {

    double r2 = r * r;
    double sin_theta = sin(theta);
    double rho = sqrt(r2 + RBH_PARAM * RBH_PARAM);

    metric[0][0] = -(1 - 2 * MASS / rho);
    metric[0][3] = 0.0;
    metric[1][1] = -1.0 / metric[0][0];
    metric[2][2] = rho * rho;
    metric[3][3] = metric[2][2] * sin_theta * sin_theta;

    *N_metric = -metric[0][0];
    *omega_metric = 0.0;

    return OK;
}

int derived_RBH_class::get_dr_metric(double dr_metric[4][4], double* dr_N, double* dr_omega, double r, double theta) {

    double metric[4][4], N, omega;

    derived_RBH_class::get_metric(metric, &N, &omega, r, theta);

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho = sqrt(r2 + RBH_PARAM * RBH_PARAM);
    double rho3 = rho * rho * rho;

    dr_metric[0][0] = -2 * MASS * r / rho3;
    dr_metric[0][3] = 0.0;
    dr_metric[1][1] = 1.0 / (metric[0][0] * metric[0][0]) * dr_metric[0][0];
    dr_metric[2][2] = 2 * r;
    dr_metric[3][3] = 2 * r * sin_theta * sin_theta;

    *dr_N = -dr_metric[0][0];
    *dr_omega = 0.0;

    return OK;

}

int derived_RBH_class::get_d2r_metric(double d2r_metric[4][4], double* d2r_N, double* d2r_omega, double r, double theta) {

    double metric[4][4], N, omega;

    derived_RBH_class::get_metric(metric, &N, &omega, r, theta);

    double dr_metric[4][4], dr_N, dr_omega;

    derived_RBH_class::get_dr_metric(dr_metric, &dr_N, &dr_omega, r, theta);

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho = sqrt(r2 + RBH_PARAM * RBH_PARAM);
    double rho3 = rho * rho * rho;
    double rho5 = rho * rho * rho * rho * rho;

    d2r_metric[0][0] = -2 * MASS / rho3 + 6 * MASS * r2 / (rho5);
    d2r_metric[0][3] = 0.0;
    d2r_metric[1][1] = 1.0 / (metric[0][0] * metric[0][0]) * d2r_metric[0][0] - 
                       2.0 / (metric[0][0] * metric[0][0] * metric[0][0]) * dr_metric[0][0] * dr_metric[0][0];
    d2r_metric[2][2] = 2.0;
    d2r_metric[3][3] = 2 * sin_theta * sin_theta;

    *d2r_N = -d2r_metric[0][0];
    *d2r_omega = 0.0;

    return OK;

}

int derived_RBH_class::get_initial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                                                        int photon, double r_obs, double theta_obs, double metric[4][4], double N, double omega) {

    *J = -J_data[photon] * sin(theta_obs);
    *p_theta = p_theta_data[photon];

    double rho = sqrt(r_obs * r_obs + RBH_PARAM * RBH_PARAM);
    double rad_potential = 1 - (1 - 2 * MASS / rho) * *J * *J / (rho * rho);

    *p_r = sqrt(rad_potential) * metric[1][1];

    return OK;
}

int derived_RBH_class::get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration) {

    double r = inter_State_vector[0 + iteration * e_State_Number];
    double rho = sqrt(r * r + RBH_PARAM * RBH_PARAM);

    double sin1 = sin(inter_State_vector[1 + iteration * e_State_Number]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(inter_State_vector[1 + iteration * e_State_Number]);
    double cos2 = cos1 * cos1;

    Derivatives[e_r + iteration * e_State_Number] = (1 - 2 * MASS / rho) * inter_State_vector[e_p_r + iteration * e_State_Number];
    Derivatives[e_theta + iteration * e_State_Number] = 1.0 / (rho * rho) * inter_State_vector[e_p_theta + iteration * e_State_Number];
    Derivatives[e_phi + iteration * e_State_Number] = J / (rho * rho * sin2);
    Derivatives[e_phi_FD + iteration * e_State_Number] = 0.0;
    Derivatives[e_p_theta + iteration * e_State_Number] = cos1 / (rho * rho * sin1 * sin2) * J * J;

    double r_term_1 = -MASS * r / (rho * rho * rho) * (1.0 / ((1 - 2 * MASS / rho) * (1 - 2 * MASS / rho)) + inter_State_vector[e_p_r + iteration * e_State_Number] * inter_State_vector[e_p_r + iteration * e_State_Number]);
    double r_term_2 = r / (rho * rho * rho * rho) * (inter_State_vector[e_p_theta + iteration * e_State_Number] * inter_State_vector[e_p_theta + iteration * e_State_Number] + J * J / sin2);

    Derivatives[e_p_r + iteration * e_State_Number] = r_term_1 + r_term_2;

    return 0;
}

bool derived_RBH_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool scatter = State_vector[e_r] > 100 && Derivatives[e_r] < 0;

    double r_horizon = sqrt(4 * MASS * MASS - RBH_PARAM * RBH_PARAM);

    bool hit_horizon_RBH = State_vector[e_r] - r_horizon < 0.05;

    return scatter || hit_horizon_RBH;

}

/***********************************************
|                                              |
| Derived Wormhole Class Functions Definitions |
|                                              |
***********************************************/

double derived_Wormhole_class::get_ISCO(Orbit_Orientation Orientation) {

    double M = MASS;

    if (SPIN < 0.016) {

        return  2 * M * (sqrt(4. / 9 * (6 * WH_REDSHIFT + 1)) * cosh(1. / 3 * acosh((1 + 9 * WH_REDSHIFT + 27. / 2 * WH_REDSHIFT * WH_REDSHIFT) / pow(6 * WH_REDSHIFT + 1, 3. / 2))) + 1. / 3);

    }
    else {

        switch (Orientation) {

        case Prograde:

            return WH_R_THROAT;

        case Retrograde:

            return WH_R_THROAT;

        default:

            std::cout << "Wrong Orbit Orientation! - Must be 'Prograde' or 'Retrograde'!" << '\n';

            return ERROR;

        }

    }

}

double derived_Wormhole_class::get_Photon_Sphere(Orbit_Orientation Orientation) {

    double M = MASS;
    double a = SPIN;

    switch (Orientation) {

    case Prograde:

        return  M / 2 * (1 + sqrt(1 + 8 * WH_REDSHIFT));

    case Retrograde:

        return M / 2 * (1 + sqrt(1 + 8 * WH_REDSHIFT));;

    default:

        std::cout << "Wrong Orbit Orientation! - Must be 'Prograde' or 'Retrograde'!" << '\n';

        return ERROR;

    }

}

int derived_Wormhole_class::get_metric(double metric[4][4], double* N_metric, double* omega, double l, double theta) {

    double r = sqrt(l * l + WH_R_THROAT * WH_R_THROAT);
    double r2 = r * r;
    double sin_theta = sin(theta);

    double exponent = -MASS / r - WH_REDSHIFT * MASS * MASS / r2;

    *N_metric = exp(exponent);
    *omega = 2 * SPIN * MASS * MASS / r2 / r;

    metric[0][0] = -*N_metric * *N_metric + r2 * *omega * *omega * sin_theta * sin_theta;
    metric[0][3] = -r2 * sin_theta * sin_theta * *omega;
    metric[1][1] = 1. / (1 - WH_R_THROAT / r);
    metric[2][2] = r2;
    metric[3][3] = r2 * sin_theta * sin_theta;

    return 0;
}

int derived_Wormhole_class::get_dr_metric(double dr_metric[4][4], double* dr_N, double* dr_omega, double l, double theta) {

    double metric[4][4], N, omega;

    derived_Wormhole_class::get_metric(metric, &N, &omega, l, theta);

    double r = sqrt(l * l + WH_R_THROAT * WH_R_THROAT);
    double r2 = r * r;
    double sin_theta = sin(theta);

    *dr_N = N * (1 / r2 + 2 * WH_REDSHIFT / (r2 * r));
    *dr_omega = -3 * omega / r;

    dr_metric[0][0] = -2 * N * *dr_N + 2 * r * omega * (omega + r * *dr_omega) * sin_theta * sin_theta;
    dr_metric[0][3] = -r * (2 * omega + r * *dr_omega) * sin_theta * sin_theta;
    dr_metric[1][1] = -1. / ((1 - WH_R_THROAT / r) * (1 - WH_R_THROAT / r)) * (WH_R_THROAT / r2);
    dr_metric[2][2] = 2 * r;
    dr_metric[3][3] = 2 * r * sin_theta * sin_theta;

    return 0;
}

int derived_Wormhole_class::get_d2r_metric(double d2r_metric[4][4], double* d2r_N, double* d2r_omega, double l, double theta) {

    double metric[4][4], N, omega;

    derived_Wormhole_class::get_metric(metric, &N, &omega, l, theta);

    double dr_metric[4][4], dr_N, dr_omega;

    derived_Wormhole_class::get_dr_metric(dr_metric, &dr_N, &dr_omega, l, theta);

    double r = sqrt(l * l + WH_R_THROAT * WH_R_THROAT);
    double r2 = r * r;
    double sin_theta = sin(theta);

    *d2r_N = dr_N * (1 / r2 + 2 * WH_REDSHIFT / (r2 * r)) - N * (2. / (r2 * r) + 6 * WH_REDSHIFT / (r2 * r2));
    *d2r_omega = -3 * dr_omega / r + 3 * omega / r2;

    d2r_metric[0][0] = -2 * dr_N * dr_N - 2 * N * *d2r_N + 2 * ((omega + r * dr_omega) * (omega + r * dr_omega) + r * omega * (dr_omega + dr_omega + r * *d2r_omega)) * sin_theta * sin_theta;
    d2r_metric[0][3] = -(2 * omega + r * dr_omega + r * (3 * dr_omega + r * *d2r_omega)) * sin_theta * sin_theta;
    d2r_metric[1][1] = 2 / ((1 - WH_R_THROAT / r) * (1 - WH_R_THROAT / r)) * ((WH_R_THROAT / r2) * (WH_R_THROAT / r2) / (1 - WH_R_THROAT / r) + WH_R_THROAT / (r2 * r));
    d2r_metric[2][2] = 2.0;
    d2r_metric[3][3] = 2 * sin_theta * sin_theta;

    return 0;
}

int derived_Wormhole_class::get_initial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                                                            int photon, double r_obs, double theta_obs, double metric[4][4], double N, double omega) {

    *J = -J_data[photon] * sin(theta_obs);
    *p_theta = p_theta_data[photon];

    double rad_potential = -(*p_theta * *p_theta + *J * *J / sin(theta_obs) / sin(theta_obs)) * N * N / r_obs / r_obs + (1 - omega * *J) * (1 - omega * *J);

    *p_r = sqrt(rad_potential) / N / sqrt(metric[1][1]) * r_obs / sqrt(pow(r_obs, 2) - pow(1, 2)) * metric[1][1];

    return 0;
}

int derived_Wormhole_class::get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration) {

    double r_throat = WH_R_THROAT;

    double sqrt_r2 = sqrt(inter_State_vector[0 + iteration * e_State_Number] * inter_State_vector[0 + iteration * e_State_Number] + r_throat * r_throat);
    double d_ell_r = inter_State_vector[0 + iteration * e_State_Number] / sqrt_r2;

    double omega = 2 * SPIN / (sqrt_r2 * sqrt_r2 * sqrt_r2);
    double d_ell_omega = -3 * omega / sqrt_r2 * d_ell_r;

    double exponent = -1 / sqrt_r2 - WH_REDSHIFT / (sqrt_r2 * sqrt_r2);
    double N = exp(exponent);
    double d_ell_N = N * (1 / (sqrt_r2 * sqrt_r2) + 2 * WH_REDSHIFT / (sqrt_r2 * sqrt_r2 * sqrt_r2)) * d_ell_r;

    double N2 = N * N;

    double sin1 = sin(inter_State_vector[e_theta + iteration * e_State_Number]);
    double sin2 = sin1 * sin1;

    Derivatives[e_r       + iteration * e_State_Number] = 1.0 / (1 + r_throat / sqrt_r2) * inter_State_vector[e_p_r + iteration * e_State_Number];
    Derivatives[e_theta   + iteration * e_State_Number] = 1.0 / (sqrt_r2 * sqrt_r2) * inter_State_vector[e_p_theta + iteration * e_State_Number];
    Derivatives[e_phi     + iteration * e_State_Number] = J / (sqrt_r2 * sqrt_r2 * sin2);
    Derivatives[e_phi_FD  + iteration * e_State_Number] = omega * (1 - omega * J) / N2;
    Derivatives[e_p_theta + iteration * e_State_Number] = (cos(inter_State_vector[1 + iteration * e_State_Number]) / sin1) / (sqrt_r2 * sqrt_r2) * J * J / sin2;

    double term_1 = -1.0 / ((1 + r_throat / sqrt_r2) * (1 + r_throat / sqrt_r2)) * r_throat * inter_State_vector[e_r + iteration * e_State_Number] / (sqrt_r2 * sqrt_r2 * sqrt_r2) * inter_State_vector[e_p_r + iteration * e_State_Number] * inter_State_vector[e_p_r + iteration * e_State_Number] / 2;
    double term_2 = 1.0 / (sqrt_r2 * sqrt_r2 * sqrt_r2) * (inter_State_vector[e_p_theta + iteration * e_State_Number] * inter_State_vector[e_p_theta + iteration * e_State_Number] + J * J / sin2) * d_ell_r;
    double term_3 = -(1.0 / (N2 * N) * d_ell_N * ((1 - omega * J) * (1 - omega * J)) - 1.0 / N2 * (-d_ell_omega * (1 - omega * J) * J));

    Derivatives[e_p_r + iteration * e_State_Number] = term_1 + term_2 + term_3;

    return OK;
}

bool derived_Wormhole_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool scatter = scatter = State_vector[e_r] > sqrt(100 * 100 + WH_R_THROAT * WH_R_THROAT) && Derivatives[0] < 0;
    bool scatter_other_side = State_vector[e_r] < -sqrt(100 * 100 + WH_R_THROAT * WH_R_THROAT);


    return scatter || scatter_other_side;

};

/******************************************************************************
|                                                                             |
| Derived Janis-Newman-Winicour Naked Singularity Class Functions Definitions |
|                                                                             |
******************************************************************************/

double derived_JNW_class::get_ISCO(Orbit_Orientation Orientation) {

    if (JNW_GAMMA > 1.0 / sqrt(5)) {

        switch (Orientation) {

        case Prograde:

            return 1.0 / JNW_GAMMA * (3.0 * JNW_GAMMA + 1.0 + sqrt(5 * JNW_GAMMA * JNW_GAMMA - 1));

        case Retrograde:

            return 1.0 / JNW_GAMMA * (3.0 * JNW_GAMMA + 1.0 - sqrt(5 * JNW_GAMMA * JNW_GAMMA - 1));

        default:

            std::cout << "Wrong Orbit Orientation! - Must be 'Prograde' or 'Retrograde'!" << '\n';

            return ERROR;

        }
    }
    else {

        return JNW_R_SINGULARITY;

    }

 };

double derived_JNW_class::get_Photon_Sphere(Orbit_Orientation Orientation) {

    if (JNW_GAMMA > 0.5) { // Weak naked singularity

        return (2 * JNW_GAMMA + 1) * JNW_R_SINGULARITY / 2;

    }
    else {

        return JNW_R_SINGULARITY;

    }

};

int derived_JNW_class::get_metric(double metric[4][4], double* N_metric, double* omega_metric, double r, double theta) {

    double r2 = r * r;
    double sin_theta = sin(theta);

    metric[0][0] = -pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA);
    metric[1][1] = -1.0 / metric[0][0];
    metric[2][2] = pow(1 - JNW_R_SINGULARITY / r, 1 - JNW_GAMMA) * r2;
    metric[3][3] = metric[2][2] * sin_theta * sin_theta;
    metric[0][3] = 0;
    metric[3][0] = metric[0][3];

    *N_metric = -metric[0][0];
    *omega_metric = 0;

    return OK;

}

int derived_JNW_class::get_dr_metric( double dr_metric[4][4], double* dr_N, double* dr_omega, double r, double theta) {

    double metric[4][4], N, omega;

    derived_JNW_class::get_metric(metric, &N, &omega, r, theta);

    double r2 = r * r;
    double sin_theta = sin(theta);

    dr_metric[0][0] = -JNW_GAMMA * pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA - 1) * JNW_R_SINGULARITY / r2;
    dr_metric[1][1] = 1.0 / (metric[0][0] * metric[0][0]) * dr_metric[0][0];;
    dr_metric[2][2] = 2 * r * pow(1 - JNW_R_SINGULARITY / r, 1 - JNW_GAMMA) + (1 - JNW_GAMMA) * pow(1 - JNW_R_SINGULARITY / r, -JNW_GAMMA) * JNW_R_SINGULARITY;
    dr_metric[3][3] = dr_metric[2][2] * sin_theta * sin_theta;
    dr_metric[0][3] = 0;
    dr_metric[3][0] = dr_metric[0][3];

    *dr_N = -dr_metric[0][0];
    *dr_omega = 0.0;

    return OK;

};

int derived_JNW_class::get_d2r_metric(double d2r_metric[4][4], double* d2r_N, double* d2r_omega, double r, double theta) {

    double metric[4][4], N, omega;

    derived_JNW_class::get_metric(metric, &N, &omega, r, theta);

    double dr_metric[4][4], dr_N, dr_omega;

    derived_JNW_class::get_dr_metric(dr_metric, &dr_N, &dr_omega, r, theta);

    double r2 = r * r;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    d2r_metric[0][0] = -JNW_GAMMA * (JNW_GAMMA - 1) * pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA - 2) * JNW_R_SINGULARITY * JNW_R_SINGULARITY / r2 / r2
                     + 2 * JNW_GAMMA * pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA - 1) * JNW_R_SINGULARITY / r2 / r;

    d2r_metric[0][3] = 0.0;

    d2r_metric[1][1] = 1.0 / (metric[0][0] * metric[0][0]) * d2r_metric[0][0] 
                     - 2.0 / (metric[0][0] * metric[0][0] * metric[0][0]) * dr_metric[0][0] * dr_metric[0][0];

    d2r_metric[2][2] = 2 * pow(1 - JNW_R_SINGULARITY / r, 1 - JNW_GAMMA) + 2 * (1 - JNW_GAMMA) * pow(1 - JNW_R_SINGULARITY / r, -JNW_GAMMA) * JNW_R_SINGULARITY / r
                     - (1 - JNW_GAMMA) * JNW_GAMMA * pow(1 - JNW_R_SINGULARITY / r, -JNW_GAMMA - 1) * JNW_R_SINGULARITY * JNW_R_SINGULARITY / r2;

    d2r_metric[3][3] = d2r_metric[2][2] * sin_theta * sin_theta;

    *d2r_N = -d2r_metric[0][0];
    *d2r_omega = 0.0;

    return OK;

}

int derived_JNW_class::get_initial_conditions_from_file(double* J, double J_data[], double* p_theta, double p_theta_data[], double* p_r,
                                                        int photon, double r_obs, double theta_obs, double metric[4][4], double N, double omega) {

    *J = -J_data[photon] * sin(theta_obs);
    *p_theta = p_theta_data[photon];

    double rad_potential = 1 - pow(1 - JNW_R_SINGULARITY / r_obs, 2 * JNW_GAMMA - 1) * *J * *J / (r_obs * r_obs);

    *p_r = sqrt(rad_potential) * metric[1][1];

    return OK;

}

int derived_JNW_class::get_EOM(double inter_State_vector[], double J, double Derivatives[], int iteration)
{

    double r = inter_State_vector[0 + iteration * e_State_Number];

    double sin1 = sin(inter_State_vector[1 + iteration * e_State_Number]);
    double sin2 = sin1 * sin1;

    double cos1 = cos(inter_State_vector[1 + iteration * e_State_Number]);
    double cos2 = cos1 * cos1;

    double pow_gamma = pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA);
    double pow_gamma_minus_1 = pow(1 - JNW_R_SINGULARITY / r, JNW_GAMMA - 1);

    Derivatives[e_r + iteration * e_State_Number] = pow_gamma * inter_State_vector[e_p_r + iteration * e_State_Number];
    Derivatives[e_theta + iteration * e_State_Number] = pow_gamma_minus_1 / (r * r) * inter_State_vector[e_p_theta + iteration * e_State_Number];
    Derivatives[e_phi + iteration * e_State_Number] = pow_gamma_minus_1 / (r * r * sin2) * J;
    Derivatives[e_phi_FD + iteration * e_State_Number] = 0.0;
    Derivatives[e_p_theta + iteration * e_State_Number] = pow_gamma_minus_1 * cos1 / (r * r * sin1 * sin2) * J * J;

    double r_term_1 = -JNW_GAMMA * JNW_R_SINGULARITY / 2 / r / r * pow_gamma_minus_1 * (1.0 / pow_gamma / pow_gamma
                    + inter_State_vector[e_p_r + iteration * e_State_Number] * inter_State_vector[e_p_r + iteration * e_State_Number]);
    double r_term_2 = 1.0 / r / r / r * pow_gamma_minus_1 * (1 - JNW_R_SINGULARITY / 2 / r * (JNW_GAMMA - 1) / (1 - JNW_R_SINGULARITY / r))
                    * (inter_State_vector[e_p_theta + iteration * e_State_Number] * inter_State_vector[e_p_theta + iteration * e_State_Number] + J * J / sin2);

    Derivatives[e_p_r + iteration * e_State_Number] = r_term_1 + r_term_2;

    return OK;

}

bool derived_JNW_class::terminate_integration(double State_vector[], double Derivatives[]) {

    bool hit_singularity = State_vector[e_r] - JNW_R_SINGULARITY < 0.05;

    bool scatter = State_vector[e_r] > 100 && Derivatives[e_r] < 0;

    return scatter || hit_singularity;

};