#include "General_GR_functions.h"

void invert_metric(double Inv_metric[4][4], const double Metric[4][4]) {

    double g2 = Metric[0][3] * Metric[0][3] - Metric[3][3] * Metric[0][0];

    Inv_metric[0][0] = -Metric[3][3] / g2;
    Inv_metric[0][3] = Metric[0][3] / g2;
    Inv_metric[3][0] = Inv_metric[0][3];
    Inv_metric[1][1] = 1. / Metric[1][1];
    Inv_metric[2][2] = 1. / Metric[2][2];
    Inv_metric[3][3] = -Metric[0][0] / g2;

}

double get_metric_det(const double Metric[4][4]) {

    return Metric[0][0] * Metric[1][1] * Metric[2][2] * Metric[3][3] * (1 - Metric[0][3] * Metric[0][3] / (Metric[0][0] * Metric[3][3]));

}

double get_eq_induced_metric_det(const double metric[4][4]) {

    return metric[0][0] * metric[1][1] * metric[3][3] * (1 - metric[0][3] * metric[0][3] / (metric[0][0] * metric[3][3]));

}

double get_4vec_norm(const double* const Vector, const double Metric[4][4], Tensor_type_enums Vector_type) {

    double inv_metric[4][4]{}, vec_norm{};

    switch (Vector_type)
    {
    case Covariant:

        invert_metric(inv_metric, Metric);

        for (int left_idx = 0; left_idx <= 3; left_idx++) {

            for (int right_idx = 0; right_idx <= 3; right_idx++) {

                vec_norm += inv_metric[left_idx][right_idx] * Vector[left_idx] * Vector[right_idx];

            }

        }

        break;

    case Contravariant:

        for (int left_idx = 0; left_idx <= 3; left_idx++) {

            for (int right_idx = 0; right_idx <= 3; right_idx++) {

                vec_norm += Metric[left_idx][right_idx] * Vector[left_idx] * Vector[right_idx];

            }

        }

        break;

    default:

        std::cout << "Unsupported Tensor type - something broke in the get_4vec_norm function!" << "\n";

        exit(ERROR);

        break;
    }

    return vec_norm;

}

void Contravariant_coord_to_ZAMO(const double Metric[4][4], const double* const Contravariant_Vector, double* const ZAMO_Vector) {

    // TODO: Make this more obvious

    double alpha = sqrt(-Metric[0][0] + Metric[0][3] * Metric[0][3] / Metric[3][3]);
    double beta = Metric[0][3] / Metric[3][3];

    double p_t = Metric[0][0] * Contravariant_Vector[0] + Metric[0][3] * Contravariant_Vector[3];
    double p_phi = Metric[3][3] * Contravariant_Vector[3] + Metric[0][3] * Contravariant_Vector[0];

    ZAMO_Vector[0] = -(1 / alpha * p_t - beta / alpha * p_phi);
    ZAMO_Vector[1] = sqrt(Metric[1][1]) * Contravariant_Vector[1];
    ZAMO_Vector[2] = sqrt(Metric[2][2]) * Contravariant_Vector[2];
    ZAMO_Vector[3] = 1 / sqrt(Metric[3][3]) * p_phi;

}

void get_intitial_conditions_from_angles(Initial_conditions_type* p_Initial_Conditions, double V_angle_cam, double H_angle_cam) {

    /*
    
     n_cam is the direction vector of the light ray in the camera frame ( {r, theta, phi} components)
     n_FIDO is the "rotated" observer (a.e. his camera "y" axis is NOT aligned with the spin axis of the central object)
    
    */

    double n_cam[3] = { cos(V_angle_cam) * cos(H_angle_cam), sin(V_angle_cam), sin(H_angle_cam) * cos(V_angle_cam)}; 
    double n_FIDO[3]{};


    n_FIDO[e_phi - 1]   =  n_cam[e_phi - 1] * cos(p_Initial_Conditions->Observer_params.cam_rotation_angle) + n_cam[e_theta - 1] * sin(p_Initial_Conditions->Observer_params.cam_rotation_angle);
    n_FIDO[e_theta - 1] = -n_cam[e_phi - 1] * sin(p_Initial_Conditions->Observer_params.cam_rotation_angle) + n_cam[e_theta - 1] * cos(p_Initial_Conditions->Observer_params.cam_rotation_angle);
    n_FIDO[e_r - 1]     =  n_cam[e_r - 1];

    double V_angle = asin(n_FIDO[e_theta - 1]);
    double H_angle = atan2(n_FIDO[e_phi - 1], n_FIDO[e_r - 1]);

    double g2, gamma, ksi, L_z, E;

    double(*metric)[4] = p_Initial_Conditions->Init_metric.Metric;

    g2 = pow(metric[0][3], 2) - metric[0][0] * metric[3][3];
    ksi = sqrt(metric[3][3] / g2);
    gamma = -metric[0][3] / metric[3][3] * ksi;

    L_z = sqrt(metric[3][3]) * sin(H_angle + 2 * M_PI) * cos(V_angle);
    E = (1 + gamma * L_z) / ksi;

    p_Initial_Conditions->Init_Momentum[e_t]     = -1;
    p_Initial_Conditions->Init_Momentum[e_phi]   = L_z / E;
    p_Initial_Conditions->Init_Momentum[e_theta] = sqrt(metric[2][2]) * sin(V_angle) / E;
    p_Initial_Conditions->Init_Momentum[e_r]     = sqrt(metric[1][1]) * cos(H_angle + 2 * M_PI) * cos(V_angle) / E;

}

void get_image_coordinates(Initial_conditions_type* p_Initial_Conditions, double* const Image_coords) {

    double(*metric)[4] = p_Initial_Conditions->Init_metric.Metric;

    double g2    = pow(metric[0][3], 2) - metric[0][0] * metric[3][3];
    double ksi   = sqrt(metric[3][3] / g2);
    double gamma = -metric[0][3] / metric[3][3] * ksi;

    double& r_0  = p_Initial_Conditions->Observer_params.distance;
    double& J    = p_Initial_Conditions->Init_Momentum[e_phi];
    double& p_th = p_Initial_Conditions->Init_Momentum[e_theta];

    Image_coords[e_x] = -r_0 *  J   / (ksi - gamma * J) / sqrt(metric[3][3]);
    Image_coords[e_y] =  r_0 * p_th / (ksi - gamma * J) / sqrt(metric[2][2]);

}

double get_redshift(const double* const State_Vector, const double* const U_source, Observer_class* const Observer) {

    // Offsetting the State_Vector pointer in this function call, so dot_product indexes the ray momentum, and not the position
    return dot_product(Observer->get_obs_velocity(), State_Vector + e_p_t, 4) / dot_product(U_source, State_Vector + e_p_t, 4);

}

void get_Lorentz_boost_matrix(double Boost_matrix[4][4], const double* const U_source) {

    double V_r     = U_source[1] / U_source[0];
    double V_theta = U_source[2] / U_source[0];
    double V_phi   = U_source[3] / U_source[0];

    double V_squared = V_r * V_r + V_theta * V_theta + V_phi * V_phi;

    double gamma = 1.0 / sqrt(1 - V_squared);

    Boost_matrix[0][0] = gamma;
    Boost_matrix[0][1] = gamma * V_r;
    Boost_matrix[0][2] = gamma * V_theta;
    Boost_matrix[0][3] = gamma * V_phi;

    for (int index = 1; index <= 3; index += 1) {

        Boost_matrix[index][0] = Boost_matrix[0][index];

    }

    Boost_matrix[1][1] = 1 + (gamma - 1) * V_r * V_r / V_squared;
    Boost_matrix[1][2] = (gamma - 1) * V_r * V_theta / V_squared;
    Boost_matrix[1][3] = (gamma - 1) * V_r * V_phi / V_squared;

    Boost_matrix[2][1] = Boost_matrix[1][2];
    Boost_matrix[3][1] = Boost_matrix[1][3];

    Boost_matrix[2][2] = 1 + (gamma - 1) * V_theta * V_theta / V_squared;
    Boost_matrix[2][3] = (gamma - 1) * V_theta * V_phi / V_squared;

    Boost_matrix[3][2] = Boost_matrix[2][3];

    Boost_matrix[3][3] = 1 + (gamma - 1) * V_phi * V_phi / V_squared;

}

bool Check_for_theta_turning_point(const double* const State_Vector, const double* const Old_State) {

    return State_Vector[e_p_theta] * Old_State[e_p_theta] < 0;

}

int compute_image_order(const int N_theta_turning_points, Initial_conditions_type* const p_Initial_Conditions) {

    int order = N_theta_turning_points;

    if (p_Initial_Conditions->Observer_params.inclination > M_PI_2) {

        order -= bool(p_Initial_Conditions->Init_Momentum[e_theta] < 0);

    }
    else {

        order -= bool(p_Initial_Conditions->Init_Momentum[e_theta] > 0);

    }

    if (order > 3) {

        order = 3;

    }

    return order * bool(order > 0);

}

void get_connection_coefficients(const Metric_type s_Metric, const Metric_type s_dr_metric, const Metric_type s_dtheta_metric, double Connection_Coeffs[4][4][4]) {

    double inv_metric[4][4]{};

    invert_metric(inv_metric, s_Metric.Metric);

  /* ==================================================================== Г^t_{..} coefficients ================================================================ */

    Connection_Coeffs[e_t][e_t][e_t]         = 0.0;
    Connection_Coeffs[e_t][e_r][e_r]         = 0.0;
    Connection_Coeffs[e_t][e_theta][e_theta] = 0.0;
    Connection_Coeffs[e_t][e_phi][e_phi]     = 0.0;

    /* ------------------------------------------------------------------ Г^t_{t,r} coefficients --------------------------------------------------------------- */

    Connection_Coeffs[e_t][e_t][e_r] = inv_metric[e_t][e_t] * s_dr_metric.Metric[e_t][e_t] / 2 + 
                                       inv_metric[e_t][e_phi] * s_dr_metric.Metric[e_t][e_phi] / 2;
    Connection_Coeffs[e_t][e_r][e_t] = Connection_Coeffs[e_t][e_t][e_r];

    /* ------------------------------------------------------------------ Г^t_{t,phi} coefficients ------------------------------------------------------------- */

    Connection_Coeffs[e_t][e_t][e_phi] = 0.0;
    Connection_Coeffs[e_t][e_phi][e_t] = 0.0;

    /* ------------------------------------------------------------------ Г^t_{t,theta} coefficients ----------------------------------------------------------- */

    Connection_Coeffs[e_t][e_t][e_theta] = inv_metric[e_t][e_t]   * s_dtheta_metric.Metric[e_t][e_t] / 2 + 
                                           inv_metric[e_t][e_phi] * s_dtheta_metric.Metric[e_t][e_phi] / 2;
    Connection_Coeffs[e_t][e_theta][e_t] = Connection_Coeffs[e_t][e_t][e_theta];

   /* ------------------------------------------------------------------ Г^t_{phi,r} coefficients -------------------------------------------------------------- */

    Connection_Coeffs[e_t][e_phi][e_r] = inv_metric[e_t][e_t]   * s_dr_metric.Metric[e_t][e_phi] / 2 + 
                                         inv_metric[e_t][e_phi] * s_dr_metric.Metric[e_phi][e_phi] / 2;
    Connection_Coeffs[e_t][e_r][e_phi] = Connection_Coeffs[e_t][e_phi][e_r];

   /* ------------------------------------------------------------------ Г^t_{phi,theta} coefficients ---------------------------------------------------------- */

    Connection_Coeffs[e_t][e_phi][e_theta] = inv_metric[e_t][e_t]   * s_dtheta_metric.Metric[e_t][e_phi] / 2 + 
                                             inv_metric[e_t][e_phi] * s_dtheta_metric.Metric[e_phi][e_phi] / 2;
    Connection_Coeffs[e_t][e_theta][e_phi] = Connection_Coeffs[e_t][e_phi][e_theta];

   /* ------------------------------------------------------------------ Г^t_{r,theta} coefficients ------------------------------------------------------------- */

    Connection_Coeffs[e_t][e_theta][e_r] = 0.0;
    Connection_Coeffs[e_t][e_r][e_theta] = 0.0;


  /* ==================================================================== Г^r_{..} coefficients ==================================================================== */

    Connection_Coeffs[e_r][e_t][e_t]         = -inv_metric[e_r][e_r] * s_dr_metric.Metric[e_t][e_t] / 2;
    Connection_Coeffs[e_r][e_r][e_r]         =  inv_metric[e_r][e_r] * s_dr_metric.Metric[e_r][e_r] / 2;
    Connection_Coeffs[e_r][e_theta][e_theta] = -inv_metric[e_r][e_r] * s_dr_metric.Metric[e_theta][e_theta] / 2;
    Connection_Coeffs[e_r][e_phi][e_phi]     = -inv_metric[e_r][e_r] * s_dr_metric.Metric[e_phi][e_phi] / 2;

   /* ------------------------------------------------------------------ Г^r_{t,phi} coefficients ---------------------------------------------------------------- */

    Connection_Coeffs[e_r][e_t][e_phi] = -inv_metric[e_r][e_r] * s_dr_metric.Metric[e_t][e_phi] / 2;
    Connection_Coeffs[e_r][e_phi][e_t] = Connection_Coeffs[e_r][e_t][e_phi];

   /* ------------------------------------------------------------------ Г^r_{t,r} coefficients ------------------------------------------------------------------ */

    Connection_Coeffs[e_r][e_t][e_r] = 0.0;
    Connection_Coeffs[e_r][e_r][e_t] = Connection_Coeffs[e_r][e_t][e_r];

   /* ------------------------------------------------------------------ Г^r_{theta,r} coefficients -------------------------------------------------------------- */

    Connection_Coeffs[e_r][e_theta][e_r] = inv_metric[e_r][e_r] * s_dtheta_metric.Metric[e_r][e_r] / 2;
    Connection_Coeffs[e_r][e_r][e_theta] = Connection_Coeffs[e_r][e_theta][e_r];

   /* ------------------------------------------------------------------ Г^r_{r,phi} coefficients ---------------------------------------------------------------- */

    Connection_Coeffs[e_r][e_r][e_phi] = 0.0;
    Connection_Coeffs[e_r][e_phi][e_r] = Connection_Coeffs[e_r][e_r][e_phi];

   /* ------------------------------------------------------------------ Г^r_{t,theta} coefficients -------------------------------------------------------------- */

    Connection_Coeffs[e_r][e_t][e_theta] = 0.0;
    Connection_Coeffs[e_r][e_theta][e_t] = Connection_Coeffs[e_r][e_t][e_theta];
     
   /* ------------------------------------------------------------------ Г^r_{phi,theta} coefficients ------------------------------------------------------------ */

    Connection_Coeffs[e_r][e_phi][e_theta] = 0.0;
    Connection_Coeffs[e_r][e_theta][e_phi] = Connection_Coeffs[e_r][e_phi][e_theta];

    /* ==================================================================== Г^theta_{..} coefficients ================================================================ */

    Connection_Coeffs[e_theta][e_t][e_t]         = -inv_metric[e_theta][e_theta] * s_dtheta_metric.Metric[e_t][e_t] / 2;
    Connection_Coeffs[e_theta][e_r][e_r]         = -inv_metric[e_theta][e_theta] * s_dtheta_metric.Metric[e_r][e_r] / 2;
    Connection_Coeffs[e_theta][e_theta][e_theta] =  inv_metric[e_theta][e_theta] * s_dtheta_metric.Metric[e_theta][e_theta] / 2;
    Connection_Coeffs[e_theta][e_phi][e_phi]     = -inv_metric[e_theta][e_theta] * s_dtheta_metric.Metric[e_phi][e_phi] / 2;

    /* ------------------------------------------------------------------ Г^theta_{t,phi} coefficients ---------------------------------------------------------------- */

    Connection_Coeffs[e_theta][e_t][e_phi] = -inv_metric[e_theta][e_theta] * s_dtheta_metric.Metric[e_t][e_phi] / 2;
    Connection_Coeffs[e_theta][e_phi][e_t] = Connection_Coeffs[e_theta][e_t][e_phi];

    /* ------------------------------------------------------------------ Г^theta_{t,r} coefficients ------------------------------------------------------------------ */

    Connection_Coeffs[e_theta][e_t][e_r] = 0.0;
    Connection_Coeffs[e_theta][e_r][e_t] = 0.0;

    /* ------------------------------------------------------------------ Г^theta_{theta,r} coefficients -------------------------------------------------------------- */

    Connection_Coeffs[e_theta][e_theta][e_r] = inv_metric[e_theta][e_theta] * s_dr_metric.Metric[e_theta][e_theta] / 2;
    Connection_Coeffs[e_theta][e_r][e_theta] = Connection_Coeffs[e_theta][e_theta][e_r];

    /* ------------------------------------------------------------------ Г^theta_{r,phi} coefficients ---------------------------------------------------------------- */

    Connection_Coeffs[e_theta][e_r][e_phi] = 0.0;
    Connection_Coeffs[e_theta][e_phi][e_r] = 0.0;

    /* ------------------------------------------------------------------ Г^theta_{t,theta} coefficients -------------------------------------------------------------- */

    Connection_Coeffs[e_theta][e_t][e_theta] = 0.0;
    Connection_Coeffs[e_theta][e_theta][e_t] = 0.0;

    /* ------------------------------------------------------------------ Г^theta_{phi,theta} coefficients ------------------------------------------------------------ */

    Connection_Coeffs[e_theta][e_phi][e_theta] = 0.0;
    Connection_Coeffs[e_theta][e_theta][e_phi] = 0.0;

  /* ==================================================================== Г^phi_{..} coefficients ==================================================================== */

    Connection_Coeffs[e_phi][e_t][e_t] = 0.0;
    Connection_Coeffs[e_phi][e_r][e_r] = 0.0;
    Connection_Coeffs[e_phi][e_theta][e_theta] = 0.0;
    Connection_Coeffs[e_phi][e_phi][e_phi] = 0.0;

   /* ------------------------------------------------------------------ Г^phi_{t,phi} coefficients ---------------------------------------------------------------- */

    Connection_Coeffs[e_phi][e_t][e_phi] = 0.0;
    Connection_Coeffs[e_phi][e_phi][e_t] = 0.0;

   /* ------------------------------------------------------------------ Г^phi_{r,phi} coefficients ---------------------------------------------------------------- */

    Connection_Coeffs[e_phi][e_r][e_phi] = inv_metric[e_phi][e_phi] * s_dr_metric.Metric[e_phi][e_phi] / 2 + 
                                           inv_metric[e_phi][e_t] * s_dr_metric.Metric[e_phi][e_t] / 2;
    Connection_Coeffs[e_phi][e_phi][e_r] = Connection_Coeffs[e_phi][e_r][e_phi];

   /* ------------------------------------------------------------------ Г^phi_{theta,phi} coefficients ------------------------------------------------------------ */

    Connection_Coeffs[e_phi][e_theta][e_phi] = inv_metric[e_phi][e_phi] * s_dtheta_metric.Metric[e_phi][e_phi] / 2 +
                                               inv_metric[e_phi][e_t] * s_dtheta_metric.Metric[e_phi][e_t] / 2;
    Connection_Coeffs[e_phi][e_phi][e_theta] = Connection_Coeffs[e_phi][e_theta][e_phi];
    
   /* ------------------------------------------------------------------ Г^phi_{t,r} coefficients ------------------------------------------------------------------ */

    Connection_Coeffs[e_phi][e_t][e_r] = inv_metric[e_phi][e_phi] * s_dr_metric.Metric[e_phi][e_t] / 2 +
                                         inv_metric[e_phi][e_t] * s_dr_metric.Metric[e_t][e_t] / 2;
    Connection_Coeffs[e_phi][e_r][e_t] = Connection_Coeffs[e_phi][e_t][e_r];

   /* ------------------------------------------------------------------ Г^phi_{theta,r} coefficients -------------------------------------------------------------- */

    Connection_Coeffs[e_phi][e_theta][e_r] = 0.0;
    Connection_Coeffs[e_phi][e_r][e_theta] = Connection_Coeffs[e_phi][e_theta][e_r];
    
   /* ------------------------------------------------------------------ Г^phi_{t,theta} coefficients ------------------------------------------------------------------ */

    Connection_Coeffs[e_phi][e_t][e_theta] = inv_metric[e_phi][e_phi] * s_dtheta_metric.Metric[e_phi][e_t] / 2 +
                                             inv_metric[e_phi][e_t] * s_dtheta_metric.Metric[e_t][e_t] / 2;
    Connection_Coeffs[e_phi][e_theta][e_t] = Connection_Coeffs[e_phi][e_t][e_theta];

}