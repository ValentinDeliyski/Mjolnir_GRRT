#include "Spacetimes.h"
#include "General_math_functions.h"
#include "General_GR_functions.h"

Return_Values Numerical_metric::load_parameters(const Metric_parameters_type* const Metric_Parameters) {

    this->Parameters = Metric_Parameters->Numerical_metric_params;

    /* Computes the outer (mass normalized) event horizon radius is Boyer-Linguist coordinates, with the ADM mass of the numerical solution (the paper labels this capital R_H).
    NOTE: The ADM angular momentum parameter is normalized to the mass. */
    this->Parameters.Horizon_radius_BL = 1 + sqrt(1 - this->Parameters.a_ADM * this->Parameters.a_ADM);

    return OK;

}

double Numerical_metric::compactify_radial_coordiante(const double r) const {

    /* The reference for this implementation is http://gravitation.web.ua.pt/node/416. */

    /* Computes the shifted radial coordinate used in the paper (they label this with little "r").
       NOTE: The input to this function "r" is normaled to the mass. */
    double r_shifted_coordinate = this->Parameters.M_ADM * (r - this->Parameters.a_ADM * this->Parameters.a_ADM / this->Parameters.Horizon_radius_BL);

    /* Comute the other shifted coordinate that the paper uses in the numerical implementation (they label this little "x"). 
       NOTE: The horizon radius here is in the "r_shifted" coordinate system. */
    double x_uncompactified = sqrt(r_shifted_coordinate * r_shifted_coordinate - this->Parameters.Horizon_radius * this->Parameters.Horizon_radius);

    return x_uncompactified / (1 + x_uncompactified);
}

inline void Numerical_metric::get_control_point_matrix(const double* const Control_vector, const int r_idx, const int theta_idx, double Control_matrix[4][4]) const {

    /* The reference for this implementation is one of the python example scrips in this document: https://hal.science/hal-03017566/document. */

    for (int radial_offset = 0; radial_offset <= 3; radial_offset++) {

        for (int theta_offset = 0; theta_offset <= 3; theta_offset++) {

            Control_matrix[radial_offset][theta_offset] = Control_vector[r_idx + radial_offset + theta_offset * (this->Parameters.Radial_grid_size + 2) + theta_idx * (this->Parameters.Radial_grid_size + 2)];
            Control_matrix[radial_offset][theta_offset] = Control_vector[r_idx + radial_offset + theta_offset * (this->Parameters.Radial_grid_size + 2) + theta_idx * (this->Parameters.Radial_grid_size + 2)];
            Control_matrix[radial_offset][theta_offset] = Control_vector[r_idx + radial_offset + theta_offset * (this->Parameters.Radial_grid_size + 2) + theta_idx * (this->Parameters.Radial_grid_size + 2)];
            Control_matrix[radial_offset][theta_offset] = Control_vector[r_idx + radial_offset + theta_offset * (this->Parameters.Radial_grid_size + 2) + theta_idx * (this->Parameters.Radial_grid_size + 2)];

        }

    }


}

void Numerical_metric::get_polynomial_basis_vector(const double natural_parameter, double* const Polynomial_basis_vector) const {

    Polynomial_basis_vector[0] = (1 - natural_parameter) * (1 - natural_parameter) * (1 - natural_parameter);
    Polynomial_basis_vector[1] = 3 * natural_parameter * natural_parameter * natural_parameter - 6 * natural_parameter * natural_parameter + 4;
    Polynomial_basis_vector[2] = -3 * natural_parameter * natural_parameter * natural_parameter + 3 * natural_parameter * natural_parameter + 3 * natural_parameter + 1;
    Polynomial_basis_vector[3] = natural_parameter * natural_parameter * natural_parameter;

}

void Numerical_metric::get_derivative_polynomial_basis_vector(const double natural_parameter, double* const Polynomial_basis_vector) const {

    Polynomial_basis_vector[0] = -3 * (1 - natural_parameter) * (1 - natural_parameter);
    Polynomial_basis_vector[1] = 9 * natural_parameter * natural_parameter - 12 * natural_parameter;
    Polynomial_basis_vector[2] = -9 * natural_parameter * natural_parameter + 6 * natural_parameter + 3;
    Polynomial_basis_vector[3] = 3 * natural_parameter * natural_parameter;

}

void Numerical_metric::get_second_derivative_polynomial_basis_vector(const double natural_parameter, double* const Polynomial_basis_vector) const {

    Polynomial_basis_vector[0] = 6 * (1 - natural_parameter);
    Polynomial_basis_vector[1] = 18 * natural_parameter - 12;
    Polynomial_basis_vector[2] = -18 * natural_parameter + 6;
    Polynomial_basis_vector[3] = 6 * natural_parameter;

}

double Numerical_metric::evaluate_single_spline(const double Control_point_matrix[4][4], 
                                         const double Radial_natural_parameter, 
                                         const double Theta_natural_parameter,
                                         Derivative_selector_enums Derivative_selector) const {

    /* -------------- Compute the basis polynomials ------------ */

    double Radial_basis_polynomial[4]{};
    double Theta_basis_polynomial[4]{};

    switch (Derivative_selector) {

    case None:

        this->get_polynomial_basis_vector(Radial_natural_parameter, Radial_basis_polynomial);
        this->get_polynomial_basis_vector(Theta_natural_parameter, Theta_basis_polynomial);

        break;

    case First_radial_derivative:

        this->get_derivative_polynomial_basis_vector(Radial_natural_parameter, Radial_basis_polynomial);
        this->get_polynomial_basis_vector(Theta_natural_parameter, Theta_basis_polynomial);

        break;

    case First_theta_derivative:

        this->get_polynomial_basis_vector(Radial_natural_parameter, Radial_basis_polynomial);
        this->get_derivative_polynomial_basis_vector(Theta_natural_parameter, Theta_basis_polynomial);

        break;

    case Second_radial_derivative:

        this->get_second_derivative_polynomial_basis_vector(Radial_natural_parameter, Radial_basis_polynomial);
        this->get_polynomial_basis_vector(Theta_natural_parameter, Theta_basis_polynomial);

        break;

    default:

        std::cout << "Unsupported numerical metric derivative! \n";
        exit(ERROR);

    }

    /* -------------- Compute the double dot product between the two polynomial basis vectors and the control point matrix. --------------*/
    // NOTE: The control point matrix acts on the "theta" basis polynomial on the right, and on the "radial" one on the left.

    double Intermediate_result[4]{};
    mat_vec_multiply_4D(Control_point_matrix, Theta_basis_polynomial, Intermediate_result);

    return dot_product(Radial_basis_polynomial, Intermediate_result, 4) / 36;

}

Metric_type Numerical_metric::evaluate_all_splines(const double* const State_Vector, Derivative_selector_enums Derivative_selector) const {

    /* -------------- The spline works with the compactified radial coordinate, so here I have to convert to that. -------------- */
    double r_compactified = this->compactify_radial_coordiante(State_Vector[e_r]);

    /* -------------- Get the upper index of the grid interval where the current photon state vector is (in both compactified radial and there directions). --------------*/
    int Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size - 1, r_compactified) - this->Parameters.Compactified_radial_grid;
    int Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size - 1, State_Vector[e_theta]) - this->Parameters.Theta_grid;

    if (0 == Radial_grid_upper_idx) {

        /* If the point we are evaluating the metric at is below the first grid point, snap it to said first grid point. */
        Radial_grid_upper_idx = 1;

    }

    /* -------------- Compute the natural parameters along the coordinate directions for the given patch. These are the arguments for the polynomial basis vectors. --------------*/
    double radial_natural_parameter = (r_compactified - this->Parameters.Compactified_radial_grid[Radial_grid_upper_idx - 1]) / (this->Parameters.Compactified_radial_grid[Radial_grid_upper_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_upper_idx - 1]);
    double theta_narual_parameter = (State_Vector[e_theta] - this->Parameters.Theta_grid[Theta_grid_upper_idx - 1]) / (this->Parameters.Theta_grid[Theta_grid_upper_idx] - this->Parameters.Theta_grid[Theta_grid_upper_idx - 1]);

    Metric_type Metric{};

    /* -------------- Compute g_tt -------------- */

    double g_tt_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.g_tt_control_vector, Radial_grid_upper_idx - 1, Theta_grid_upper_idx - 1, g_tt_control_point_matrix);

    Metric.Metric[e_t][e_t] = this->evaluate_single_spline(g_tt_control_point_matrix, radial_natural_parameter, theta_narual_parameter, Derivative_selector);

    /* -------------- Compute g_rr -------------- */

    double g_rr_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.g_rr_control_vector, Radial_grid_upper_idx - 1, Theta_grid_upper_idx - 1, g_rr_control_point_matrix);

    Metric.Metric[e_r][e_r] = this->evaluate_single_spline(g_rr_control_point_matrix, radial_natural_parameter, theta_narual_parameter, Derivative_selector);

    /* -------------- Compute g_thth -------------- */

    double g_thth_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.g_thth_control_vector, Radial_grid_upper_idx - 1, Theta_grid_upper_idx - 1, g_thth_control_point_matrix);

    Metric.Metric[e_theta][e_theta] = this->evaluate_single_spline(g_thth_control_point_matrix, radial_natural_parameter, theta_narual_parameter, Derivative_selector);

    /* -------------- Compute g_phiphi -------------- */

    double g_phiphi_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.g_phiphi_control_vector, Radial_grid_upper_idx - 1, Theta_grid_upper_idx - 1, g_phiphi_control_point_matrix);

    Metric.Metric[e_phi][e_phi] = this->evaluate_single_spline(g_phiphi_control_point_matrix, radial_natural_parameter, theta_narual_parameter, Derivative_selector);

    /* -------------- Compute g_tphi -------------- */

    double g_tphi_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.g_tphi_control_vector, Radial_grid_upper_idx - 1, Theta_grid_upper_idx - 1, g_tphi_control_point_matrix);

    Metric.Metric[e_t][e_phi] = this->evaluate_single_spline(g_tphi_control_point_matrix, radial_natural_parameter, theta_narual_parameter, Derivative_selector);
    Metric.Metric[e_phi][e_t] = Metric.Metric[e_t][e_phi];

    return Metric;

}

Metric_type Numerical_metric::comute_Minkowski_metric(const double* const State_Vector) const {

    Metric_type s_Minkowski_metric{};

    const double& r = State_Vector[e_r];
    const double sin_theta = sin(State_Vector[e_theta]);

    s_Minkowski_metric.Metric[e_t][e_t] = -1;
    s_Minkowski_metric.Metric[e_r][e_r] = 1;
    s_Minkowski_metric.Metric[e_theta][e_theta] = r * r;
    s_Minkowski_metric.Metric[e_phi][e_phi] = r * r * sin_theta * sin_theta;

    s_Minkowski_metric.Lapse_function = 1;
    s_Minkowski_metric.Shift_function = 0;

    return s_Minkowski_metric;

}

Metric_type Numerical_metric::comute_dr_Minkowski_metric(const double* const State_Vector) const {

    Metric_type s_dr_Minkowski_metric{};

    const double& r = State_Vector[e_r];
    const double sin_theta = sin(State_Vector[e_theta]);

    s_dr_Minkowski_metric.Metric[e_t][e_t] = 0;
    s_dr_Minkowski_metric.Metric[e_r][e_r] = 0;
    s_dr_Minkowski_metric.Metric[e_theta][e_theta] = 2 *r;
    s_dr_Minkowski_metric.Metric[e_phi][e_phi] = 2 * r * sin_theta * sin_theta;

    s_dr_Minkowski_metric.Lapse_function = 0;
    s_dr_Minkowski_metric.Shift_function = 0;

    return s_dr_Minkowski_metric;

}

Metric_type Numerical_metric::comute_dtheta_Minkowski_metric(const double* const State_Vector) const {

    Metric_type s_dtheta_Minkowski_metric{};

    const double& r = State_Vector[e_r];
    const double sin_theta = sin(State_Vector[e_theta]);
    const double cos_theta = cos(State_Vector[e_theta]);

    s_dtheta_Minkowski_metric.Metric[e_t][e_t] = 0;
    s_dtheta_Minkowski_metric.Metric[e_r][e_r] = 0;
    s_dtheta_Minkowski_metric.Metric[e_theta][e_theta] = 0;
    s_dtheta_Minkowski_metric.Metric[e_phi][e_phi] = 2 * r * r * sin_theta * cos_theta;

    s_dtheta_Minkowski_metric.Lapse_function = 0;
    s_dtheta_Minkowski_metric.Shift_function = 0;

    return s_dtheta_Minkowski_metric;

}

Metric_type Numerical_metric::comute_d2r_Minkowski_metric(const double* const State_Vector) const {

    Metric_type s_d2r_Minkowski_metric{};

    const double& r = State_Vector[e_r];
    const double sin_theta = sin(State_Vector[e_theta]);
    const double cos_theta = cos(State_Vector[e_theta]);

    s_d2r_Minkowski_metric.Metric[e_t][e_t] = 0;
    s_d2r_Minkowski_metric.Metric[e_r][e_r] = 0;
    s_d2r_Minkowski_metric.Metric[e_theta][e_theta] = 2;
    s_d2r_Minkowski_metric.Metric[e_phi][e_phi] = 2 * sin_theta * sin_theta;

    s_d2r_Minkowski_metric.Lapse_function = 0;
    s_d2r_Minkowski_metric.Shift_function = 0;

    return s_d2r_Minkowski_metric;

}

Metric_type Numerical_metric::compute_metric_components_from_spline(const double* const State_Vector, Derivative_selector_enums Derivative_selector) const {

    /* -------------- The spline works with the compactified radial coordinate, so here I have to convert to that. -------------- */
    double r_compactified = this->compactify_radial_coordiante(State_Vector[e_r]);
    double r_shifted_coordinate = this->Parameters.M_ADM * (State_Vector[e_r] - this->Parameters.a_ADM * this->Parameters.a_ADM / this->Parameters.Horizon_radius_BL);

    /* -------------- Get the upper index of the grid interval where the current photon state vector is (in both compactified radial and there directions). --------------*/
    int Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size - 1, r_compactified) - this->Parameters.Compactified_radial_grid;
    int Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size - 1, State_Vector[e_theta]) - this->Parameters.Theta_grid;

    if (0 == Radial_grid_upper_idx) {

        /* If the point we are evaluating the metric at is below the first grid point, snap it to said first grid point. */
        Radial_grid_upper_idx = 1;

    }
    else if (Radial_grid_upper_idx >= this->Parameters.Radial_grid_size - 1) {

        switch (Derivative_selector) {

        case First_radial_derivative:

            return this->comute_dr_Minkowski_metric(State_Vector);

        case Second_radial_derivative:

            return this->comute_d2r_Minkowski_metric(State_Vector);

        case First_theta_derivative:

            return this->comute_dtheta_Minkowski_metric(State_Vector);

        default:

            return this->comute_Minkowski_metric(State_Vector);

        }

    }

    Metric_type temp_Metric_1{};
    Metric_type temp_Metric_2{};

    Metric_type Corrected_metric{};

    double derivative_correction_factor_1 = 1.0;
    double derivative_correction_factor_2 = 1.0;

    switch (Derivative_selector) {

    case First_radial_derivative:

        Corrected_metric = this->evaluate_all_splines(State_Vector, First_radial_derivative);

        /* ------- This is correcting by the factor d(natural_parameter)/d(x_compactified) ------- */
        derivative_correction_factor_1 = 1 / (this->Parameters.Compactified_radial_grid[Radial_grid_upper_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_upper_idx - 1]);

        /* ------- This is correcting by the factor d(x_compactified)/d(x_uncompactified) ------- */
        derivative_correction_factor_1 *= (1 - r_compactified) * (1 - r_compactified);

        /* ------- This is correcting by the factor d(x_uncompactified)/d(r_shifted_coordinate) ------- */
        derivative_correction_factor_1 *= r_shifted_coordinate * (1 - r_compactified) / r_compactified;

        /* ------- This is correcting by the factor d(r_shifted_coordinate)/d(r_BL) ------- */
        derivative_correction_factor_1 *= this->Parameters.M_ADM;

        for (int left_idx = 0; left_idx <= 3; left_idx++) {

            for (int right_idx = 0; right_idx <= 3; right_idx++) {

                Corrected_metric.Metric[left_idx][right_idx] *= derivative_correction_factor_1;

            }

        }

        return Corrected_metric;

    case First_theta_derivative:

        Corrected_metric = this->evaluate_all_splines(State_Vector, First_theta_derivative);

        /* ------- This is correcting by the factor d(natural_parameter)/d(theta) ------- */
        derivative_correction_factor_1 = 1 / (this->Parameters.Theta_grid[Theta_grid_upper_idx] - this->Parameters.Theta_grid[Theta_grid_upper_idx - 1]);

        for (int left_idx = 0; left_idx <= 3; left_idx++) {

            for (int right_idx = 0; right_idx <= 3; right_idx++) {

                Corrected_metric.Metric[left_idx][right_idx] *= derivative_correction_factor_1;

            }

        }

        return Corrected_metric;

    case Second_radial_derivative:

        temp_Metric_1 = this->evaluate_all_splines(State_Vector, Second_radial_derivative);
        temp_Metric_2 = this->evaluate_all_splines(State_Vector, First_radial_derivative);

        /* ------- This is correcting by the factor d(natural_parameter)/d(x_compactified) ------- */
        derivative_correction_factor_1 = 1 / (this->Parameters.Compactified_radial_grid[Radial_grid_upper_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_upper_idx - 1]);

        /* ------- This is correcting by the factor d(x_compactified)/d(x_uncompactified) ------- */
        derivative_correction_factor_1 *= (1 - r_compactified) * (1 - r_compactified);

        /* ------- This is correcting by the factor d(x_uncompactified)/d(r_shifted_coordinate) ------- */
        derivative_correction_factor_1 *= r_shifted_coordinate * (1 - r_compactified) / r_compactified;

        /* ------- This is correcting by the factor d(r_shifted_coordinate)/d(r_BL) ------- */
        derivative_correction_factor_1 *= this->Parameters.M_ADM;

        /* ------- This is the first term in d(derivative_correction_factor_1)/d(r_BL) ------- */
        derivative_correction_factor_2  = -this->Parameters.M_ADM * ((1 - r_compactified) * (1 - r_compactified) * (1 - r_compactified) + 3 * r_compactified * (1 - r_compactified) * (1 - r_compactified)) / r_compactified / r_compactified * r_shifted_coordinate * derivative_correction_factor_1;
        derivative_correction_factor_2 += this->Parameters.M_ADM * derivative_correction_factor_1 / r_shifted_coordinate;

        for (int left_idx = 0; left_idx <= 3; left_idx++) {

            for (int right_idx = 0; right_idx <= 3; right_idx++) {

                temp_Metric_1.Metric[left_idx][right_idx] *= derivative_correction_factor_1 * derivative_correction_factor_1;
                temp_Metric_2.Metric[left_idx][right_idx] *= derivative_correction_factor_2;

                Corrected_metric.Metric[left_idx][right_idx] = temp_Metric_1.Metric[left_idx][right_idx] + temp_Metric_2.Metric[left_idx][right_idx];

            }

        }

        return Corrected_metric;

    default:

        return this->evaluate_all_splines(State_Vector, None);

    }

}

Metric_type Numerical_metric::get_metric(const double* const State_Vector) const {

    return this->compute_metric_components_from_spline(State_Vector, None);

}

Metric_type Numerical_metric::get_dr_metric(const double* const State_Vector) const {

    return this->compute_metric_components_from_spline(State_Vector, First_radial_derivative);
}

Metric_type Numerical_metric::get_dtheta_metric(const double* const State_Vector) const {

    return this->compute_metric_components_from_spline(State_Vector, First_theta_derivative);
}

Metric_type Numerical_metric::get_d2r_metric(const double* const State_Vector) const {

    return this->compute_metric_components_from_spline(State_Vector, Second_radial_derivative);
}

int Numerical_metric::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

    return 0;
}

void Numerical_metric::get_EOM(double State_Vector[], double Derivatives[]) const {

    /* ----------- Temporary matrix, used to store intermediate calculations ----------- */
    double temp_matrix[4][4]{};

    Metric_type Metric = this->get_metric(State_Vector);
    Metric_type dr_Metric = this->get_dr_metric(State_Vector);
    Metric_type dtheta_Metric = this->get_dtheta_metric(State_Vector);

    double inv_metric[4][4]{};
    invert_metric(inv_metric, Metric.Metric);

    double dr_inv_metric[4][4]{};
    matrix_matrix_multiply(dr_Metric.Metric, inv_metric, temp_matrix);
    matrix_matrix_multiply(inv_metric, temp_matrix, dr_inv_metric);

    double dtheta_inv_metric[4][4]{};
    matrix_matrix_multiply(dtheta_Metric.Metric, inv_metric, temp_matrix);
    matrix_matrix_multiply(inv_metric, temp_matrix, dtheta_inv_metric);

    /* ----------- There is a minus sign infront of the whole expression for the derivative of an inverse of a matrix ----------- */
    for (int left_idx = 0; left_idx <= 3; left_idx++) {

        for (int right_idx = 0; right_idx <= 3; right_idx++) {

            dr_inv_metric[left_idx][right_idx] *= -1;
            dtheta_inv_metric[left_idx][right_idx] *= -1;

        }

    }

    for (int right_idx = 0; right_idx <= 3; right_idx++) {

        *(Derivatives + e_t) += inv_metric[e_t][right_idx] * State_Vector[right_idx + 4];
        *(Derivatives + e_r) += inv_metric[e_r][right_idx] * State_Vector[right_idx + 4];
        *(Derivatives + e_theta) += inv_metric[e_theta][right_idx] * State_Vector[right_idx + 4];
        *(Derivatives + e_phi) += inv_metric[e_phi][right_idx] * State_Vector[right_idx + 4];

        for (int left_idx = 0; left_idx <= 3; left_idx++) {

            *(Derivatives + e_p_r) += -1. / 2 * dr_inv_metric[left_idx][right_idx] * State_Vector[left_idx + 4] * State_Vector[right_idx + 4];
            *(Derivatives + e_p_theta) += -1. / 2 * dtheta_inv_metric[left_idx][right_idx] * State_Vector[left_idx + 4] * State_Vector[right_idx + 4];

        }

    }

    *(Derivatives + e_p_t) = 0.0;
    *(Derivatives + e_p_phi) = 0.0;
   
}

bool Numerical_metric::terminate_integration(double State_vector[], double Derivatives[]) {

    bool scatter = State_vector[e_r] > 30 && Derivatives[e_r] < 0;

    bool hit_horizon = State_vector[e_r] - this->Parameters.Horizon_radius_BL < 1e-1;

    return scatter || hit_horizon;
};
