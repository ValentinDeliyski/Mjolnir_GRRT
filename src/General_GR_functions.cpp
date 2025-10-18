#include "General_GR_functions.h"

void get_derivative_of_inverse_metric(const double inv_Metric[4][4], const double deriv_Metric[4][4], double deriv_inv_Metric[4][4]) {

    memset(deriv_inv_Metric, 0, sizeof(double) * 16);

    double temp_matrix[4][4]{};

    matrix_matrix_multiply(deriv_Metric, inv_Metric, temp_matrix);
    matrix_matrix_multiply(inv_Metric, temp_matrix, deriv_inv_Metric);

    for (int left_idx = 0; left_idx < 4; left_idx++) {

        for (int right_idx = 0; right_idx < 4; right_idx++) {

            deriv_inv_Metric[left_idx][right_idx] *= -1;

        }

    }

}

void invert_metric(double Inv_metric[4][4], const double Metric[4][4]) {

    double g2 = Metric[e_t][e_phi] * Metric[e_t][e_phi] - Metric[e_phi][e_phi] * Metric[e_t][e_t];

    Inv_metric[e_t][e_t] = -Metric[e_phi][e_phi] / g2;
    Inv_metric[e_t][e_phi] = Metric[0][e_phi] / g2;
    Inv_metric[e_phi][e_t] = Inv_metric[0][e_phi];
    Inv_metric[e_r][e_r] = 1. / Metric[e_r][e_r];
    Inv_metric[e_theta][e_theta] = 1. / Metric[e_theta][e_theta];
    Inv_metric[e_phi][e_phi] = -Metric[e_t][e_t] / g2;

}

double get_metric_det(const double Metric[4][4]) {

    return Metric[e_t][e_t] * Metric[e_r][e_r] * Metric[e_theta][e_theta] * Metric[e_phi][e_phi] * (1 - Metric[e_t][e_phi] * Metric[e_t][e_phi] / (Metric[e_t][e_t] * Metric[e_phi][e_phi]));

}

double get_eq_induced_metric_det(const double Metric[4][4]) {

    return Metric[e_t][e_t] * Metric[e_r][e_r] * Metric[e_phi][e_phi] * (1 - Metric[e_t][e_phi] * Metric[e_t][e_phi] / (Metric[e_t][e_t] * Metric[e_phi][e_phi]));

}

double get_complex_4vec_norm(const std::complex<double>* const Vector, const double Metric[4][4], Tensor_type_enums Vector_type) {

    double inv_metric[4][4]{};
    std::complex<double> vec_norm_squrated{};

    switch (Vector_type)
    {
    case Covariant:

        invert_metric(inv_metric, Metric);

        for (int left_idx = 0; left_idx <= 3; left_idx++) {

            for (int right_idx = 0; right_idx <= 3; right_idx++) {

                vec_norm_squrated += inv_metric[left_idx][right_idx] * Vector[left_idx] * std::conj(Vector[right_idx]);

            }

        }

        break;

    case Contravariant:

        for (int left_idx = 0; left_idx <= 3; left_idx++) {

            for (int right_idx = 0; right_idx <= 3; right_idx++) {

                vec_norm_squrated += Metric[left_idx][right_idx] * Vector[left_idx] * std::conj(Vector[right_idx]);

            }

        }

        break;

    default:

        throw std::runtime_error("Unsupported Tensor type - something broke in the get_4vec_norm function!");

    }

    return std::sqrt(std::abs(vec_norm_squrated));

}

void Normalize_complex_vector(std::complex<double>* const Vector, const double Metric[4][4], Tensor_type_enums Vector_type) {

    double Vector_norm = get_complex_4vec_norm(Vector, Metric, Vector_type);

    if (isinf(1.0 / Vector_norm) || isnan(1.0 / Vector_norm)) { std::cout << "Invalid vector norm! \n"; exit(ERROR); }

    for (int idx = 0; idx < 4; idx++) {

        Vector[idx] *= 1.0 / Vector_norm;

    }

}

void Contravariant_coord_to_ZAMO(const Metric_type* const p_Metric, const double* const Contravariant_Vector, double* const ZAMO_Vector) {

    ZAMO_Vector[e_t] = p_Metric->Lapse_function * Contravariant_Vector[e_t];
    ZAMO_Vector[e_r] = sqrt(p_Metric->Metric[e_r][e_r])  * Contravariant_Vector[e_r];
    ZAMO_Vector[e_theta] = sqrt(p_Metric->Metric[e_theta][e_theta]) * Contravariant_Vector[e_theta];
    ZAMO_Vector[e_phi] = sqrt(p_Metric->Metric[e_phi][e_phi]) * (Contravariant_Vector[e_phi] - p_Metric->Shift_function * Contravariant_Vector[e_t]);

}

void ZAMO_to_Contravariant_coord(const Metric_type* const p_Metric, const double* const ZAMO_Vector, double* const Contravariant_Vector) {

    Contravariant_Vector[e_t] = ZAMO_Vector[e_t] / p_Metric->Lapse_function;
    Contravariant_Vector[e_r] = ZAMO_Vector[e_r] / sqrt(p_Metric->Metric[e_r][e_r]);
    Contravariant_Vector[e_theta] = ZAMO_Vector[e_theta] / sqrt(p_Metric->Metric[e_theta][e_theta]);
    Contravariant_Vector[e_phi] = p_Metric->Shift_function / p_Metric->Lapse_function * ZAMO_Vector[e_t] + ZAMO_Vector[e_phi] / sqrt(p_Metric->Metric[e_phi][e_phi]);

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

    g2 = pow(metric[e_t][e_phi], 2) - metric[e_t][e_t] * metric[e_phi][e_phi];
    ksi = sqrt(metric[e_phi][e_phi] / g2);
    gamma = -metric[e_t][e_phi] / metric[e_phi][e_phi] * ksi;

    L_z = sqrt(metric[3][3]) * sin(H_angle + 2 * M_PI) * cos(V_angle);
    E = (1 + gamma * L_z) / ksi;

    p_Initial_Conditions->Init_Momentum[e_t]     = -1;
    p_Initial_Conditions->Init_Momentum[e_phi]   = L_z / E;
    p_Initial_Conditions->Init_Momentum[e_theta] = sqrt(metric[e_theta][e_theta]) * sin(V_angle) / E;
    p_Initial_Conditions->Init_Momentum[e_r]     = sqrt(metric[e_r][e_r]) * cos(H_angle + 2 * M_PI) * cos(V_angle) / E;

}

void get_image_coordinates(Initial_conditions_type* p_Initial_Conditions, double* const Image_coords) {

    double(*metric)[4] = p_Initial_Conditions->Init_metric.Metric;

    double g2 = pow(metric[e_t][e_phi], 2) - metric[e_t][e_t] * metric[e_phi][e_phi];
    double ksi = sqrt(metric[e_phi][e_phi] / g2);
    double gamma = -metric[e_t][e_phi] / metric[e_phi][e_phi] * ksi;

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

void get_Lorentz_boost_matrix(double Boost_matrix[4][4], const double* const U_source_ZAMO, bool Inverse_boost) {

    double V_r     = U_source_ZAMO[e_r] / U_source_ZAMO[e_t];
    double V_theta = U_source_ZAMO[e_theta] / U_source_ZAMO[e_t];
    double V_phi   = U_source_ZAMO[e_phi] / U_source_ZAMO[e_t];

    if (Inverse_boost) {

        V_r *= -1;
        V_theta *= -1;
        V_phi *= -1;

    }

    double V_squared = V_r * V_r + V_theta * V_theta + V_phi * V_phi;

    double gamma = 1.0 / sqrt(1 - V_squared);

    Boost_matrix[e_t][e_t]     =  gamma;
    Boost_matrix[e_t][e_r]     = -gamma * V_r;
    Boost_matrix[e_t][e_theta] = -gamma * V_theta;
    Boost_matrix[e_t][e_phi]   = -gamma * V_phi;

    for (int index = e_r; index <= e_phi; index += 1) {

        Boost_matrix[index][e_t] = Boost_matrix[e_t][index];

    }

    Boost_matrix[e_r][e_r]     = 1 + (gamma - 1) * V_r * V_r / V_squared;
    Boost_matrix[e_r][e_theta] = (gamma - 1) * V_r * V_theta / V_squared;
    Boost_matrix[e_r][e_phi]   = (gamma - 1) * V_r * V_phi / V_squared;

    Boost_matrix[e_theta][e_r] = Boost_matrix[e_r][e_theta];
    Boost_matrix[e_phi][e_r]   = Boost_matrix[e_r][e_phi];

    Boost_matrix[e_theta][e_theta] = 1 + (gamma - 1) * V_theta * V_theta / V_squared;
    Boost_matrix[e_theta][e_phi]   = (gamma - 1) * V_theta * V_phi / V_squared;

    Boost_matrix[e_phi][e_theta] = Boost_matrix[e_theta][e_phi];

    Boost_matrix[e_phi][e_phi] = 1 + (gamma - 1) * V_phi * V_phi / V_squared;

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

void get_initial_conditions_from_image_coords(Initial_conditions_type* p_Initial_Conditions, double Image_X_coord, double Image_Y_coord) {

    double (*metric)[4] = p_Initial_Conditions->Init_metric.Metric;

    double inv_metric[4][4]{};
    invert_metric(inv_metric, p_Initial_Conditions->Init_metric.Metric);

    double g2 = pow(metric[e_t][e_phi], 2) - metric[e_t][e_t] * metric[e_phi][e_phi];
    double ksi = sqrt(metric[e_phi][e_phi] / g2);
    double gamma = -metric[e_t][e_phi] / metric[e_phi][e_phi] * ksi;

    double& r_0 = p_Initial_Conditions->Observer_params.distance;

    p_Initial_Conditions->Init_Momentum[e_t] = -1;
    p_Initial_Conditions->Init_Momentum[e_phi] = sqrt(metric[e_phi][e_phi]) * ksi * Image_X_coord / (sqrt(metric[e_phi][e_phi]) * gamma * Image_X_coord - r_0);
    p_Initial_Conditions->Init_Momentum[e_theta] = Image_Y_coord * (ksi - gamma * p_Initial_Conditions->Init_Momentum[e_phi]) * sqrt(metric[e_theta][e_theta]) / r_0;

    double effective_rad_potential{};

    for (int idx_1 = e_t; idx_1 <= e_phi; idx_1++) {

        for (int idx_2 = e_t; idx_2 <= e_phi; idx_2++) {

            if (e_r == idx_1 || e_r == idx_2) { continue; }

            effective_rad_potential += inv_metric[idx_1][idx_2] * p_Initial_Conditions->Init_Momentum[idx_1] * p_Initial_Conditions->Init_Momentum[idx_2];

        }
    }

    p_Initial_Conditions->Init_Momentum[e_r] = sqrt(-effective_rad_potential / inv_metric[e_r][e_r]);

}

std::complex<double> get_Penrose_Walker_constant(const double* const State_Vector, const Spacetime_Base_Class* const p_Spacetime, const std::complex<double>* const Polarization_Vector) {

    double Contravariant_momentum[4]{};
    double inv_metric[4][4]{};

    Metric_type s_Metric = p_Spacetime->get_metric(State_Vector);

    invert_metric(inv_metric, s_Metric.Metric);

    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            Contravariant_momentum[left_idx] += inv_metric[left_idx][right_idx] * State_Vector[right_idx + e_p_t];

        }

    }

    const double& p_t     = Contravariant_momentum[e_t];
    const double& p_r     = Contravariant_momentum[e_r];
    const double& p_theta = Contravariant_momentum[e_theta];
    const double& p_phi   = Contravariant_momentum[e_phi];

    std::complex<double> Kappa_1 = sqrt(-s_Metric.Metric[e_theta][e_theta] * s_Metric.Metric[e_r][e_r] * s_Metric.Metric[e_t][e_t]) * (p_t * Polarization_Vector[e_r] - p_r * Polarization_Vector[e_t]);
    std::complex<double> Kappa_2 = pow(s_Metric.Metric[e_theta][e_theta], 3.0 / 2) * sin(State_Vector[e_theta]) * (p_theta * Polarization_Vector[e_phi] - p_phi * Polarization_Vector[e_theta]);

    return Kappa_1 - complex_i * Kappa_2;

}