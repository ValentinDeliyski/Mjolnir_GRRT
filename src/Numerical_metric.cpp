#include "Spacetimes.h"
#include "General_math_functions.h"
#include "General_GR_functions.h"

Return_Values Numerical_metric::load_parameters(const Metric_parameters_type* const p_Metric_Parameters) {

    if (isnan(p_Metric_Parameters->Scattering_radius) || isinf(p_Metric_Parameters->Scattering_radius) || p_Metric_Parameters->Scattering_radius < 0) {

        std::cout << "Invalid value for the scattering radius: " << p_Metric_Parameters->Scattering_radius << "\n";

        return ERROR;
    }

    if (isnan(p_Metric_Parameters->Min_distance_to_singular_point) || isinf(p_Metric_Parameters->Min_distance_to_singular_point)) {

        std::cout << "Invalid value for the distance to the throat: " << p_Metric_Parameters->Min_distance_to_singular_point << "\n";

        return ERROR;
    }

    this->Min_distance_to_singular_point = p_Metric_Parameters->Min_distance_to_singular_point;
    this->Scattering_radius = p_Metric_Parameters->Scattering_radius;

    this->Parameters = p_Metric_Parameters->Numerical_metric_params;

    this->GSL_workspace = gsl_eigen_nonsymmv_alloc(8);
    this->EOM_Eigenvalues = gsl_vector_complex_alloc(8);
    this->EOM_Eigenvectors = gsl_matrix_complex_alloc(8, 8);
    this->EOM_Jacobian_gsl_matrix = gsl_matrix_alloc(8, 8);

    return OK;

}

double Numerical_metric::compactify_radial_coordiante(const double r) const {

    /* Comute the shifted coordinate that the paper (http://gravitation.web.ua.pt/node/416) uses in the numerical implementation (they label this little "x"). */
    const double x_uncompactified = sqrt(r * r - this->Parameters.Horizon_radius * this->Parameters.Horizon_radius);

    return x_uncompactified / (1 + x_uncompactified);
}

inline void Numerical_metric::get_control_point_matrix(const double* const Control_vector, const int r_idx, const int theta_idx, double Control_matrix[4][4]) const {

    /* The reference for this implementation is one of the python example scrips in this document: https://hal.science/hal-03017566/document. */

    for (int radial_offset = 0; radial_offset < 4; radial_offset++) {

        for (int theta_offset = 0; theta_offset < 4; theta_offset++) {

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

    case Second_mixed_derivative:

        this->get_derivative_polynomial_basis_vector(Radial_natural_parameter, Radial_basis_polynomial);
        this->get_derivative_polynomial_basis_vector(Theta_natural_parameter, Theta_basis_polynomial);

        break;

    case Second_theta_derivative:

        this->get_polynomial_basis_vector(Radial_natural_parameter, Radial_basis_polynomial);
        this->get_second_derivative_polynomial_basis_vector(Theta_natural_parameter, Theta_basis_polynomial);

        break;

    default:

        std::cout << "Unsupported numerical metric derivative! \n";
        exit(ERROR);

    }

    /* -------------- Compute the double dot product between the two polynomial basis vectors and the control point matrix. --------------*/
    // NOTE: The control point matrix acts on the "theta" basis polynomial on the right, and on the "radial" one on the left.

    double Intermediate_result[4]{};
    mat_vec_multiply_4D(Control_point_matrix, Theta_basis_polynomial, Intermediate_result);

    return dot_product(Radial_basis_polynomial, Intermediate_result, 4) / 36.0;

}

Numerical_metric_potentials_type Numerical_metric::evaluate_all_splines(const double radial_natural_param, const double theta_narual_param, int Radial_grid_idx, int Theta_grid_idx, Derivative_selector_enums Derivative_selector) const {

    Numerical_metric_potentials_type s_Metric_potentials{};

    /* -------------- Compute F_0 -------------- */

    double F_0_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.F_0_control_vector, Radial_grid_idx - 1, Theta_grid_idx - 1, F_0_control_point_matrix);

    s_Metric_potentials.F_0 = this->evaluate_single_spline(F_0_control_point_matrix, radial_natural_param, theta_narual_param, Derivative_selector);

    /* -------------- Compute F_1 -------------- */

    double F_1_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.F_1_control_vector, Radial_grid_idx - 1, Theta_grid_idx - 1, F_1_control_point_matrix);

    s_Metric_potentials.F_1 = this->evaluate_single_spline(F_1_control_point_matrix, radial_natural_param, theta_narual_param, Derivative_selector);

    /* -------------- Compute F_2 -------------- */

    double F_2_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.F_2_control_vector, Radial_grid_idx - 1, Theta_grid_idx - 1, F_2_control_point_matrix);

    s_Metric_potentials.F_2 = this->evaluate_single_spline(F_2_control_point_matrix, radial_natural_param, theta_narual_param, Derivative_selector);

    /* -------------- Compute W -------------- */

    double W_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.W_control_vector, Radial_grid_idx - 1, Theta_grid_idx - 1, W_control_point_matrix);

    s_Metric_potentials.W = this->evaluate_single_spline(W_control_point_matrix, radial_natural_param, theta_narual_param, Derivative_selector);

    return s_Metric_potentials;

}

Numerical_metric_potentials_type Numerical_metric::compute_metric_components_from_spline(Spline_arguments_type s_Splnie_args, Derivative_selector_enums Derivative_selector) const {

    Numerical_metric_potentials_type temp_Potentials_1{};
    Numerical_metric_potentials_type temp_Potentials_2{};
    Numerical_metric_potentials_type temp_Potentials_3{};

    Numerical_metric_potentials_type Corrected_Potentials{};

    double derivative_correction_factor_1 = 1.0;
    double derivative_correction_factor_2 = 1.0;
    double derivative_correction_factor_3 = 1.0;
    double derivative_correction_factor_4 = 1.0;
    double derivative_correction_factor_5 = 1.0;

    switch (Derivative_selector) {

    case First_radial_derivative:

        Corrected_Potentials = this->evaluate_all_splines(s_Splnie_args.Radial_natural_param, s_Splnie_args.Theta_natural_param, s_Splnie_args.Radial_idx, s_Splnie_args.Theta_idx, Derivative_selector);

        /* ------- This is correcting by the factor d(natural_parameter)/d(x_compactified) ------- */
        derivative_correction_factor_1 = 1 / (this->Parameters.Compactified_radial_grid[s_Splnie_args.Radial_idx] - this->Parameters.Compactified_radial_grid[s_Splnie_args.Radial_idx - 1]);

        /* ------- This is correcting by the factor d(x_compactified)/d(x_uncompactified) ------- */
        derivative_correction_factor_1 *= (1 - s_Splnie_args.r_coord_compactified) * (1 - s_Splnie_args.r_coord_compactified);

        /* ------- This is correcting by the factor d(x_uncompactified)/d(r) ------- */
        derivative_correction_factor_1 *= s_Splnie_args.r_coord * (1 - s_Splnie_args.r_coord_compactified) / s_Splnie_args.r_coord_compactified;

        Corrected_Potentials.F_0 *= derivative_correction_factor_1;
        Corrected_Potentials.F_1 *= derivative_correction_factor_1;
        Corrected_Potentials.F_2 *= derivative_correction_factor_1;
        Corrected_Potentials.W   *= derivative_correction_factor_1;

        return Corrected_Potentials;

    case First_theta_derivative:

        Corrected_Potentials = this->evaluate_all_splines(s_Splnie_args.Radial_natural_param, s_Splnie_args.Theta_natural_param, s_Splnie_args.Radial_idx, s_Splnie_args.Theta_idx, Derivative_selector);

        /* ------- This is correcting by the factor d(natural_parameter)/d(theta) ------- */
        derivative_correction_factor_1 = 1 / (this->Parameters.Theta_grid[s_Splnie_args.Theta_idx] - this->Parameters.Theta_grid[s_Splnie_args.Theta_idx - 1]);

        Corrected_Potentials.F_0 *= derivative_correction_factor_1;
        Corrected_Potentials.F_1 *= derivative_correction_factor_1;
        Corrected_Potentials.F_2 *= derivative_correction_factor_1;
        Corrected_Potentials.W   *= derivative_correction_factor_1;

        return Corrected_Potentials;

    case Second_radial_derivative:

        temp_Potentials_1 = this->evaluate_all_splines(s_Splnie_args.Radial_natural_param, s_Splnie_args.Theta_natural_param, s_Splnie_args.Radial_idx, s_Splnie_args.Theta_idx, Second_radial_derivative);
        temp_Potentials_2 = this->evaluate_all_splines(s_Splnie_args.Radial_natural_param, s_Splnie_args.Theta_natural_param, s_Splnie_args.Radial_idx, s_Splnie_args.Theta_idx, First_radial_derivative);
        temp_Potentials_3 = temp_Potentials_2;

        /* ------- This is correcting by the factor d(natural_parameter)/d(x_compactified) ------- */
        derivative_correction_factor_1 = 1 / (this->Parameters.Compactified_radial_grid[s_Splnie_args.Radial_idx] - this->Parameters.Compactified_radial_grid[s_Splnie_args.Radial_idx - 1]);

        /* ------- This is correcting by the factor d(x_compactified)/d(x_uncompactified) ------- */
        derivative_correction_factor_2 = (1 - s_Splnie_args.r_coord_compactified) * (1 - s_Splnie_args.r_coord_compactified);

        /* ------- This is correcting by the factor d(x_uncompactified)/d(r) ------- */
        derivative_correction_factor_3 = s_Splnie_args.r_coord * (1 - s_Splnie_args.r_coord_compactified) / s_Splnie_args.r_coord_compactified;

        /* ------- This is correcting by the factor d^2(dx_compactified)/d(x_uncompactified)^2 ------- */
        derivative_correction_factor_4 = -2 * (1 - s_Splnie_args.r_coord_compactified) * (1 - s_Splnie_args.r_coord_compactified) * (1 - s_Splnie_args.r_coord_compactified);

        /* ------- This is correcting by the factor d^2(dx_uncompactified)/d(r)^2 ------- */
        derivative_correction_factor_5 = -this->Parameters.Horizon_radius * this->Parameters.Horizon_radius * (1 - s_Splnie_args.r_coord_compactified) * (1 - s_Splnie_args.r_coord_compactified) * (1 - s_Splnie_args.r_coord_compactified) / s_Splnie_args.r_coord_compactified / s_Splnie_args.r_coord_compactified / s_Splnie_args.r_coord_compactified;

        temp_Potentials_1.F_0 *= (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3) * (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3);
        temp_Potentials_1.F_1 *= (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3) * (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3);
        temp_Potentials_1.F_2 *= (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3) * (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3);
        temp_Potentials_1.W   *= (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3) * (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3);

        temp_Potentials_2.F_0 *= derivative_correction_factor_1 * derivative_correction_factor_4 * derivative_correction_factor_3 * derivative_correction_factor_3;
        temp_Potentials_2.F_1 *= derivative_correction_factor_1 * derivative_correction_factor_4 * derivative_correction_factor_3 * derivative_correction_factor_3;
        temp_Potentials_2.F_2 *= derivative_correction_factor_1 * derivative_correction_factor_4 * derivative_correction_factor_3 * derivative_correction_factor_3;
        temp_Potentials_2.W   *= derivative_correction_factor_1 * derivative_correction_factor_4 * derivative_correction_factor_3 * derivative_correction_factor_3;

        temp_Potentials_3.F_0 *= derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_5;
        temp_Potentials_3.F_1 *= derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_5;
        temp_Potentials_3.F_2 *= derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_5;
        temp_Potentials_3.W   *= derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_5;

        Corrected_Potentials.F_0 = temp_Potentials_1.F_0 + temp_Potentials_2.F_0 + temp_Potentials_3.F_0; 
        Corrected_Potentials.F_1 = temp_Potentials_1.F_1 + temp_Potentials_2.F_1 + temp_Potentials_3.F_1; 
        Corrected_Potentials.F_2 = temp_Potentials_1.F_2 + temp_Potentials_2.F_2 + temp_Potentials_3.F_2; 
        Corrected_Potentials.W   = temp_Potentials_1.W   + temp_Potentials_2.W   + temp_Potentials_3.W;  

        return Corrected_Potentials;

    case Second_mixed_derivative:

        Corrected_Potentials = this->evaluate_all_splines(s_Splnie_args.Radial_natural_param, s_Splnie_args.Theta_natural_param, s_Splnie_args.Radial_idx, s_Splnie_args.Theta_idx, Derivative_selector);

        /* ------- This is correcting by the factor d(natural_parameter)/d(x_compactified) ------- */
        derivative_correction_factor_1 = 1 / (this->Parameters.Compactified_radial_grid[s_Splnie_args.Radial_idx] - this->Parameters.Compactified_radial_grid[s_Splnie_args.Radial_idx - 1]);

        /* ------- This is correcting by the factor d(x_compactified)/d(x_uncompactified) ------- */
        derivative_correction_factor_1 *= (1 - s_Splnie_args.r_coord_compactified) * (1 - s_Splnie_args.r_coord_compactified);

        /* ------- This is correcting by the factor d(x_uncompactified)/d(r) ------- */
        derivative_correction_factor_1 *= s_Splnie_args.r_coord * (1 - s_Splnie_args.r_coord_compactified) / s_Splnie_args.r_coord_compactified;

        /* ------- This is correcting by the factor d(natural_parameter)/d(theta) ------- */
        derivative_correction_factor_1 *= 1 / (this->Parameters.Theta_grid[s_Splnie_args.Theta_idx] - this->Parameters.Theta_grid[s_Splnie_args.Theta_idx - 1]);

        Corrected_Potentials.F_0 *= derivative_correction_factor_1;
        Corrected_Potentials.F_1 *= derivative_correction_factor_1;
        Corrected_Potentials.F_2 *= derivative_correction_factor_1;
        Corrected_Potentials.W   *= derivative_correction_factor_1;

        return Corrected_Potentials;

    case Second_theta_derivative:

        Corrected_Potentials = this->evaluate_all_splines(s_Splnie_args.Radial_natural_param, s_Splnie_args.Theta_natural_param, s_Splnie_args.Radial_idx, s_Splnie_args.Theta_idx, Derivative_selector);

        /* ------- This is correcting by the factor d(natural_parameter)/d(theta) ------- */
        derivative_correction_factor_1 = 1 / (this->Parameters.Theta_grid[s_Splnie_args.Theta_idx] - this->Parameters.Theta_grid[s_Splnie_args.Theta_idx - 1]);

        Corrected_Potentials.F_0 *= derivative_correction_factor_1 * derivative_correction_factor_1;
        Corrected_Potentials.F_1 *= derivative_correction_factor_1 * derivative_correction_factor_1;
        Corrected_Potentials.F_2 *= derivative_correction_factor_1 * derivative_correction_factor_1;
        Corrected_Potentials.W   *= derivative_correction_factor_1 * derivative_correction_factor_1;

        return Corrected_Potentials;

    default:

        return this->evaluate_all_splines(s_Splnie_args.Radial_natural_param, s_Splnie_args.Theta_natural_param, s_Splnie_args.Radial_idx, s_Splnie_args.Theta_idx, Derivative_selector);

    }

}

Metric_type Numerical_metric::get_metric(const double* const State_Vector) const {

    /* ---------------- This is a wrapper function for compatability with the radiative transfer part of the code ---------------- */

    const double r_compactified = this->compactify_radial_coordiante(State_Vector[e_r]);

    const int Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size, r_compactified) - this->Parameters.Compactified_radial_grid;
    const int Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size, State_Vector[e_theta]) - this->Parameters.Theta_grid;

    return this->get_metric(State_Vector, Radial_grid_upper_idx, Theta_grid_upper_idx);

}

Metric_type Numerical_metric::get_metric(const double* const State_Vector, int Radial_grid_idx, int Theta_grid_idx) const {

    /* -------- The metric antatz is from https://arxiv.org/pdf/1501.04319. */

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    const double sin_theta = sin(theta);
    const double N = 1 - this->Parameters.Horizon_radius / r;
    
    /* -------------------------------------------------------------- Evaluate the metric potentials spline -------------------------------------------------------------- */

    Spline_arguments_type s_Spline_args{};

    s_Spline_args.r_coord = r;
    s_Spline_args.r_coord_compactified = compactify_radial_coordiante(s_Spline_args.r_coord);
    s_Spline_args.Radial_idx = Radial_grid_idx;
    s_Spline_args.Theta_idx = Theta_grid_idx;

    s_Spline_args.Radial_natural_param = (s_Spline_args.r_coord_compactified - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1]) 
                                      / (this->Parameters.Compactified_radial_grid[Radial_grid_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1]);

    s_Spline_args.Theta_natural_param = (State_Vector[e_theta] - this->Parameters.Theta_grid[Theta_grid_idx - 1]) 
                                      / (this->Parameters.Theta_grid[Theta_grid_idx] - this->Parameters.Theta_grid[Theta_grid_idx - 1]);

    const Numerical_metric_potentials_type s_Potentials = this->compute_metric_components_from_spline(s_Spline_args, None);

    /* ------------------------------------------------------------------------------------------------------------------------------------------------------------------- */

    /* ---- References for the sake of readability ---- */

    const double& W        = s_Potentials.W;
    const double exp_2F_0 = exp(2 * s_Potentials.F_0);
    const double exp_2F_1 = exp(2 * s_Potentials.F_1);
    const double exp_2F_2 = exp(2 * s_Potentials.F_2);

    /* ------------------------------------------------ */

    Metric_type s_Metric{};

    switch (this->Parameters.e_Anzatz) {

    case e_Anzatz_2:

        s_Metric.Metric[e_t][e_t] = -exp_2F_0 * N + exp_2F_2 * r * r * sin_theta * sin_theta * W * W;
        s_Metric.Metric[e_t][e_phi] = -exp_2F_2 * W * r * r * sin_theta * sin_theta;

        break;

    default:

        s_Metric.Metric[e_t][e_t] = -exp_2F_0 * N + exp_2F_2 * sin_theta * sin_theta * W * W;
        s_Metric.Metric[e_t][e_phi] = -exp_2F_2 * W * r * sin_theta * sin_theta;

        break;

    }

    s_Metric.Metric[e_phi][e_t] = s_Metric.Metric[e_t][e_phi];
    s_Metric.Metric[e_r][e_r] = exp_2F_1 / N;
    s_Metric.Metric[e_theta][e_theta] = exp_2F_1 * r * r;
    s_Metric.Metric[e_phi][e_phi] = exp_2F_2 * r * r * sin_theta * sin_theta;

    s_Metric.Lapse_function = sqrt(-s_Metric.Metric[e_t][e_t] + s_Metric.Metric[e_t][e_phi] * s_Metric.Metric[e_t][e_phi] / s_Metric.Metric[e_phi][e_phi]);
    s_Metric.Shift_function = -s_Metric.Metric[e_t][e_phi] / s_Metric.Metric[e_phi][e_phi];

    return s_Metric;

}

Metric_type Numerical_metric::get_dr_metric(const double* const State_Vector) const {

    /* ---------------- This is a wrapper function for compatability with the radiative transfer part of the code ---------------- */

    const double r_compactified = this->compactify_radial_coordiante(State_Vector[e_r]);

    const int Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size, r_compactified) - this->Parameters.Compactified_radial_grid;
    const int Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size, State_Vector[e_theta]) - this->Parameters.Theta_grid;

    return this->get_dr_metric(State_Vector, Radial_grid_upper_idx, Theta_grid_upper_idx);

}

Metric_type Numerical_metric::get_dr_metric(const double* const State_Vector, int Radial_grid_idx, int Theta_grid_idx) const {

    /* -------- The metric antatz is from https://arxiv.org/pdf/1501.04319. */

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    const double sin_theta = sin(theta);
    const double N    = 1 - this->Parameters.Horizon_radius / r;
    const double dr_N = this->Parameters.Horizon_radius / r / r;

    /* -------------------------------------------------------------- Evaluate the metric potentials spline -------------------------------------------------------------- */

    Spline_arguments_type s_Spline_args{};

    s_Spline_args.r_coord = r;
    s_Spline_args.r_coord_compactified = compactify_radial_coordiante(s_Spline_args.r_coord);
    s_Spline_args.Radial_idx = Radial_grid_idx;
    s_Spline_args.Theta_idx = Theta_grid_idx;

    s_Spline_args.Radial_natural_param = (s_Spline_args.r_coord_compactified - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1])
                                       / (this->Parameters.Compactified_radial_grid[Radial_grid_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1]);

    s_Spline_args.Theta_natural_param = (State_Vector[e_theta] - this->Parameters.Theta_grid[Theta_grid_idx - 1])
                                      / (this->Parameters.Theta_grid[Theta_grid_idx] - this->Parameters.Theta_grid[Theta_grid_idx - 1]);

    const Numerical_metric_potentials_type s_Potentials = this->compute_metric_components_from_spline(s_Spline_args, None);
    const Numerical_metric_potentials_type s_dr_Potentials = this->compute_metric_components_from_spline(s_Spline_args, First_radial_derivative);

    /* ------------------------------------------------------------------------------------------------------------------------------------------------------------------- */

    /* ---- References for the sake of readability ---- */

    const double& W = s_Potentials.W;
    const double& exp_2F_0 = exp(2 * s_Potentials.F_0);
    const double& exp_2F_1 = exp(2 * s_Potentials.F_1);
    const double& exp_2F_2 = exp(2 * s_Potentials.F_2);

    const double& dr_W   = s_dr_Potentials.W;
    const double& dr_F_0 = s_dr_Potentials.F_0;
    const double& dr_F_1 = s_dr_Potentials.F_1;
    const double& dr_F_2 = s_dr_Potentials.F_2;

    /* ------------------------------------------------ */

    Metric_type s_dr_Metric{};

    switch (this->Parameters.e_Anzatz) {

    case e_Anzatz_2:

        s_dr_Metric.Metric[e_t][e_t] = -exp_2F_0 * (2 * N * dr_F_0 + dr_N) + 2 * exp_2F_2 * r * sin_theta * sin_theta * W * (r * W * dr_F_2 + W + r * dr_W);
        s_dr_Metric.Metric[e_t][e_phi] = -exp_2F_2 * r * sin_theta * sin_theta * (2 * r * W * dr_F_2 + 2 * W + r * dr_W);

        break;

    default:

        s_dr_Metric.Metric[e_t][e_t] = -exp_2F_0 * (2 * N * dr_F_0 + dr_N) + 2 * exp_2F_2 * sin_theta * sin_theta * W * (W * dr_F_2 + dr_W);
        s_dr_Metric.Metric[e_t][e_phi] = -exp_2F_2 * sin_theta * sin_theta * (2 * r * W * dr_F_2 + W + r * dr_W);

        break;

    }

    s_dr_Metric.Metric[e_phi][e_t] = s_dr_Metric.Metric[e_t][e_phi];
    s_dr_Metric.Metric[e_r][e_r] = exp_2F_1 / N * (2 * dr_F_1 - dr_N / N);
    s_dr_Metric.Metric[e_theta][e_theta] = 2 * r * exp_2F_1 * (r * dr_F_1 + 1);
    s_dr_Metric.Metric[e_phi][e_phi] = 2 * r * exp_2F_2 * (r * dr_F_2 + 1) * sin_theta * sin_theta;

    return s_dr_Metric;

}

Metric_type Numerical_metric::get_dtheta_metric(const double* const State_Vector) const {

    /* ---------------- This is a wrapper function for compatability with the radiative transfer part of the code ---------------- */

    const double& r_compactified = this->compactify_radial_coordiante(State_Vector[e_r]);

    const int& Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size, r_compactified) - this->Parameters.Compactified_radial_grid;
    const int& Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size, State_Vector[e_theta]) - this->Parameters.Theta_grid;

    return this->get_dtheta_metric(State_Vector, Radial_grid_upper_idx, Theta_grid_upper_idx);

}

Metric_type Numerical_metric::get_dtheta_metric(const double* const State_Vector, int Radial_grid_idx, int Theta_grid_idx) const {

    /* -------- The metric antatz is from https://arxiv.org/pdf/1501.04319. */

    /* -------------------------------------------------------------- Evaluate the metric potentials spline -------------------------------------------------------------- */

    Spline_arguments_type s_Spline_args{};

    s_Spline_args.r_coord = State_Vector[e_r];
    s_Spline_args.r_coord_compactified = compactify_radial_coordiante(s_Spline_args.r_coord);
    s_Spline_args.Radial_idx = Radial_grid_idx;
    s_Spline_args.Theta_idx = Theta_grid_idx;

    s_Spline_args.Radial_natural_param = (s_Spline_args.r_coord_compactified - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1])
                                       / (this->Parameters.Compactified_radial_grid[Radial_grid_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1]);

    s_Spline_args.Theta_natural_param = (State_Vector[e_theta] - this->Parameters.Theta_grid[Theta_grid_idx - 1])
                                      / (this->Parameters.Theta_grid[Theta_grid_idx] - this->Parameters.Theta_grid[Theta_grid_idx - 1]);

    const Numerical_metric_potentials_type s_Potentials = this->compute_metric_components_from_spline(s_Spline_args, None);
    const Numerical_metric_potentials_type s_dtheta_Potentials = this->compute_metric_components_from_spline(s_Spline_args, First_theta_derivative);

    /* ------------------------------------------------------------------------------------------------------------------------------------------------------------------- */

    /* ---- References for the sake of readability ---- */

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    const double& W        = s_Potentials.W;
    const double exp_2F_0 = exp(2 * s_Potentials.F_0);
    const double exp_2F_1 = exp(2 * s_Potentials.F_1);
    const double exp_2F_2 = exp(2 * s_Potentials.F_2);

    const double& dtheta_W   = s_dtheta_Potentials.W;
    const double& dtheta_F_0 = s_dtheta_Potentials.F_0;
    const double& dtheta_F_1 = s_dtheta_Potentials.F_1;
    const double& dtheta_F_2 = s_dtheta_Potentials.F_2;

    /* ------------------------------------------------- */

    const double sin_theta = sin(theta);
    const double cos_theta = cos(theta);
    const double N = 1 - this->Parameters.Horizon_radius / r;

    Metric_type s_dtheta_Metric{};

    switch (this->Parameters.e_Anzatz) {

    case e_Anzatz_2:

        s_dtheta_Metric.Metric[e_t][e_t] = -2 * exp_2F_0 * N * dtheta_F_0 + 2 * exp_2F_2 * r * r * W * sin_theta * (W * sin_theta * dtheta_F_2 + sin_theta * dtheta_W + W * cos_theta);
        s_dtheta_Metric.Metric[e_t][e_phi] = -exp_2F_2 * r * r * sin_theta * (2 * sin_theta * W * dtheta_F_2 + sin_theta * dtheta_W + 2 * W * cos_theta);

        break;

    default:

        s_dtheta_Metric.Metric[e_t][e_t] = -2 * exp_2F_0 * N * dtheta_F_0 + 2 * exp_2F_2 * W * sin_theta * (W * sin_theta * dtheta_F_2 + sin_theta * dtheta_W + W * cos_theta);
        s_dtheta_Metric.Metric[e_t][e_phi] = -exp_2F_2 * r * sin_theta * (2 * sin_theta * W * dtheta_F_2 + sin_theta * dtheta_W + 2 * W * cos_theta);

        break;

    }

    s_dtheta_Metric.Metric[e_phi][e_t] = s_dtheta_Metric.Metric[e_t][e_phi];
    s_dtheta_Metric.Metric[e_r][e_r] = 2 * exp_2F_1 / N * dtheta_F_1;
    s_dtheta_Metric.Metric[e_theta][e_theta] = 2 * r * r * exp_2F_1 * dtheta_F_1;
    s_dtheta_Metric.Metric[e_phi][e_phi] = 2 * r * r * exp_2F_2 * sin_theta * (sin_theta * dtheta_F_2 + cos_theta);

    return s_dtheta_Metric;

}

Metric_type Numerical_metric::get_d2r_metric(const double* const State_Vector) const {

    /* ---------------- This is a wrapper function for compatability with the radiative transfer part of the code ---------------- */

    const double r_compactified = this->compactify_radial_coordiante(State_Vector[e_r]);

    const int Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size, r_compactified) - this->Parameters.Compactified_radial_grid;
    const int Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size, State_Vector[e_theta]) - this->Parameters.Theta_grid;

    return this->get_d2r_metric(State_Vector, Radial_grid_upper_idx, Theta_grid_upper_idx);

}

Metric_type Numerical_metric::get_d2r_metric(const double* const State_Vector, int Radial_grid_idx, int Theta_grid_idx) const {

    /* -------- The metric antatz is from https://arxiv.org/pdf/1501.04319. */

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    const double& sin_theta = sin(theta);
    const double& N = 1 - this->Parameters.Horizon_radius / r;
    const double& dr_N = this->Parameters.Horizon_radius / r / r;
    const double& d2r_N = -2 * this->Parameters.Horizon_radius / r / r / r;

    /* -------------------------------------------------------------- Evaluate the metric potentials spline -------------------------------------------------------------- */

    Spline_arguments_type s_Spline_args{};

    s_Spline_args.r_coord = r;
    s_Spline_args.r_coord_compactified = compactify_radial_coordiante(s_Spline_args.r_coord);
    s_Spline_args.Radial_idx = Radial_grid_idx;
    s_Spline_args.Theta_idx = Theta_grid_idx;

    s_Spline_args.Radial_natural_param = (s_Spline_args.r_coord_compactified - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1])
                                       / (this->Parameters.Compactified_radial_grid[Radial_grid_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1]);

    s_Spline_args.Theta_natural_param = (State_Vector[e_theta] - this->Parameters.Theta_grid[Theta_grid_idx - 1])
                                      / (this->Parameters.Theta_grid[Theta_grid_idx] - this->Parameters.Theta_grid[Theta_grid_idx - 1]);

    const Numerical_metric_potentials_type s_Potentials = this->compute_metric_components_from_spline(s_Spline_args, None);
    const Numerical_metric_potentials_type s_dr_Potentials = this->compute_metric_components_from_spline(s_Spline_args, First_radial_derivative);
    const Numerical_metric_potentials_type s_d2r_Potentials = this->compute_metric_components_from_spline(s_Spline_args, Second_radial_derivative);

    /* ------------------------------------------------------------------------------------------------------------------------------------------------------------------- */

    /* ---- References for the sake of readability ---- */

    const double& W = s_Potentials.W;
    const double exp_2F_0 = exp(2 * s_Potentials.F_0);
    const double exp_2F_1 = exp(2 * s_Potentials.F_1);
    const double exp_2F_2 = exp(2 * s_Potentials.F_2);

    const double& dr_W = s_dr_Potentials.W;
    const double& dr_F_0 = s_dr_Potentials.F_0;
    const double& dr_F_1 = s_dr_Potentials.F_1;
    const double& dr_F_2 = s_dr_Potentials.F_2;

    const double& d2r_W = s_d2r_Potentials.W;
    const double& d2r_F_0 = s_d2r_Potentials.F_0;
    const double& d2r_F_1 = s_d2r_Potentials.F_1;
    const double& d2r_F_2 = s_d2r_Potentials.F_2;

    /* ------------------------------------------------ */

    Metric_type s_d2r_Metric{};

    switch (this->Parameters.e_Anzatz) {

    case e_Anzatz_2:

        /* TODO: fix these */

        s_d2r_Metric.Metric[e_t][e_t] = -exp_2F_0 * (d2r_N + 4 * dr_N * dr_F_0 + 2 * N * d2r_F_0 + 4 * N * dr_F_0 * dr_F_0)
                                      + 2 * exp_2F_2 * r * r * sin_theta * sin_theta * (W * W * d2r_F_2 + 2 * W * W * dr_F_2 * dr_F_2 + 4 * W * dr_W * dr_F_2 + dr_W * dr_W + W * d2r_W)
                                      + 2 * exp_2F_2 * sin_theta * sin_theta * (4 * r * W * W *dr_F_2 + 4 * r * W * dr_W + W * W);
        s_d2r_Metric.Metric[e_t][e_phi] = -exp_2F_2 * sin_theta * sin_theta * (4 * r * dr_W + 8 * r * W * dr_F_2 + 4 * r * r * dr_W * dr_F_2 + 2 * r * r * W * d2r_F_2 + r * r * d2r_W + 4 * r * r * W * dr_F_2 * dr_F_2 + 2 * W);

        break;

    default:

        s_d2r_Metric.Metric[e_t][e_t] = -exp_2F_0 * (d2r_N + 4 * dr_N * dr_F_0 + 2 * N * d2r_F_0 + 4 * N * dr_F_0 * dr_F_0)
            + 2 * exp_2F_2 * sin_theta * sin_theta * (W * W * d2r_F_2 + 2 * W * W * dr_F_2 * dr_F_2 + 4 * W * dr_W * dr_F_2 + dr_W * dr_W + W * d2r_W);
        s_d2r_Metric.Metric[e_t][e_phi] = -exp_2F_2 * sin_theta * sin_theta * (2 * dr_W + 4 * W * dr_F_2 + 4 * r * dr_W * dr_F_2 + 2 * r * W * d2r_F_2 + 4 * r * W * dr_F_2 * dr_F_2 + r * d2r_W);

        break;

    }

    s_d2r_Metric.Metric[e_phi][e_t] = s_d2r_Metric.Metric[e_t][e_phi];
    s_d2r_Metric.Metric[e_r][e_r] = exp_2F_1 / N * (4 * dr_F_1 * dr_F_1 + 2 * d2r_F_1 - 4 * dr_F_1 * dr_N / N - d2r_N / N + 2 * dr_N * dr_N / N / N);
    s_d2r_Metric.Metric[e_theta][e_theta] = 2 * exp_2F_1 * (1 + 4 * r * dr_F_1 + 2 * r * r * dr_F_1 * dr_F_1 + r * r * d2r_F_1);
    s_d2r_Metric.Metric[e_phi][e_phi] = 2 * exp_2F_2 * (1 + 4 * r * dr_F_2 + 2 * r * r * dr_F_2 * dr_F_2 + r * r * d2r_F_2) * sin_theta * sin_theta;

    return s_d2r_Metric;
}

Metric_type Numerical_metric::get_d2theta_metric(const double* const State_Vector) const {

    /* ---------------- This is a wrapper function for compatability with the radiative transfer part of the code ---------------- */

    const double r_compactified = this->compactify_radial_coordiante(State_Vector[e_r]);

    const int Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size, r_compactified) - this->Parameters.Compactified_radial_grid;
    const int Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size, State_Vector[e_theta]) - this->Parameters.Theta_grid;

    return this->get_d2theta_metric(State_Vector, Radial_grid_upper_idx, Theta_grid_upper_idx);

}

Metric_type Numerical_metric::get_d2theta_metric(const double* const State_Vector, int Radial_grid_idx, int Theta_grid_idx) const {

    /* -------- The metric antatz is from https://arxiv.org/pdf/1501.04319. */

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    const double& sin_theta = sin(theta);
    const double& cos_theta = cos(theta);
    const double& N = 1 - this->Parameters.Horizon_radius / r;

    /* -------------------------------------------------------------- Evaluate the metric potentials spline -------------------------------------------------------------- */

    Spline_arguments_type s_Spline_args{};

    s_Spline_args.r_coord = r;
    s_Spline_args.r_coord_compactified = compactify_radial_coordiante(s_Spline_args.r_coord);
    s_Spline_args.Radial_idx = Radial_grid_idx;
    s_Spline_args.Theta_idx = Theta_grid_idx;

    s_Spline_args.Radial_natural_param = (s_Spline_args.r_coord_compactified - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1])
                                       / (this->Parameters.Compactified_radial_grid[Radial_grid_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1]);

    s_Spline_args.Theta_natural_param = (State_Vector[e_theta] - this->Parameters.Theta_grid[Theta_grid_idx - 1])
                                      / (this->Parameters.Theta_grid[Theta_grid_idx] - this->Parameters.Theta_grid[Theta_grid_idx - 1]);

    const Numerical_metric_potentials_type s_Potentials = this->compute_metric_components_from_spline(s_Spline_args, None);
    const Numerical_metric_potentials_type s_dtheta_Potentials = this->compute_metric_components_from_spline(s_Spline_args, First_theta_derivative);
    const Numerical_metric_potentials_type s_d2theta_Potentials = this->compute_metric_components_from_spline(s_Spline_args, Second_theta_derivative);

    /* ------------------------------------------------------------------------------------------------------------------------------------------------------------------- */

    /* ---- References for the sake of readability ---- */

    const double& W = s_Potentials.W;
    const double exp_2F_0 = exp(2 * s_Potentials.F_0);
    const double exp_2F_1 = exp(2 * s_Potentials.F_1);
    const double exp_2F_2 = exp(2 * s_Potentials.F_2);

    const double& dtheta_W   = s_dtheta_Potentials.W;
    const double& dtheta_F_0 = s_dtheta_Potentials.F_0;
    const double& dtheta_F_1 = s_dtheta_Potentials.F_1;
    const double& dtheta_F_2 = s_dtheta_Potentials.F_2;

    const double& d2theta_W   = s_d2theta_Potentials.W;
    const double& d2theta_F_0 = s_d2theta_Potentials.F_0;
    const double& d2theta_F_1 = s_d2theta_Potentials.F_1;
    const double& d2theta_F_2 = s_d2theta_Potentials.F_2;

    /* ------------------------------------------------ */

    Metric_type s_d2theta_Metric{};

    s_d2theta_Metric.Metric[e_t][e_t] = -2 * N * exp_2F_0 * (d2theta_F_0 + 2 * dtheta_F_0 * dtheta_F_0)
        + 2 * exp_2F_2 * sin_theta * sin_theta * (W * W * d2theta_F_2 + 2 * W * W * dtheta_F_2 * dtheta_F_2 + 4 * W * dtheta_W * dtheta_F_2 + dtheta_W * dtheta_W + W * d2theta_W - W * W)
        + 8 * exp_2F_2 * sin_theta * cos_theta * (dtheta_F_2 * W * W + W * dtheta_W) + 2 * exp_2F_2 * W * W * cos_theta * cos_theta;

    s_d2theta_Metric.Metric[e_t][e_phi] = -exp_2F_2 * r * sin_theta * sin_theta * (2 * d2theta_F_2 * W + 4 * dtheta_F_2 * dtheta_F_2 * W + 4 * dtheta_F_2 * dtheta_W + d2theta_W - 2 * W)
        - exp_2F_2 * r * sin_theta * cos_theta * (dtheta_W + 2 * dtheta_F_2 * W) - 2 * exp_2F_2 * W * r * cos_theta * cos_theta;

    switch (this->Parameters.e_Anzatz) {

    case e_Anzatz_2:

        s_d2theta_Metric.Metric[e_t][e_t] *= r * r;
        s_d2theta_Metric.Metric[e_t][e_phi] *= r;

        break;

    default:

        break;

    }

    s_d2theta_Metric.Metric[e_phi][e_t] = s_d2theta_Metric.Metric[e_t][e_phi];
    s_d2theta_Metric.Metric[e_r][e_r] = 2 * exp_2F_1 / N * (2 * dtheta_F_1 * dtheta_F_1 + d2theta_F_1);
    s_d2theta_Metric.Metric[e_theta][e_theta] = 2 * exp_2F_1 * r * r * (2 * dtheta_F_1 * dtheta_F_1 + d2theta_F_1);
    s_d2theta_Metric.Metric[e_phi][e_phi] = 2 * exp_2F_2 * r * r * (sin_theta * sin_theta * d2theta_F_2 + 2 * dtheta_F_2 * sin_theta * cos_theta + 2 * sin_theta * sin_theta * dtheta_F_2 * dtheta_F_2 + cos_theta * cos_theta - sin_theta * sin_theta);

    return s_d2theta_Metric;
}

Metric_type Numerical_metric::get_dr_dtheta_metric(const double* const State_Vector) const {

    /* ---------------- This is a wrapper function for compatability with the radiative transfer part of the code ---------------- */

    const double r_compactified = this->compactify_radial_coordiante(State_Vector[e_r]);

    const int Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size, r_compactified) - this->Parameters.Compactified_radial_grid;
    const int Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size, State_Vector[e_theta]) - this->Parameters.Theta_grid;

    return this->get_dr_dtheta_metric(State_Vector, Radial_grid_upper_idx, Theta_grid_upper_idx);

}

Metric_type Numerical_metric::get_dr_dtheta_metric(const double* const State_Vector, int Radial_grid_idx, int Theta_grid_idx) const {

    /* -------- The metric antatz is from https://arxiv.org/pdf/1501.04319. */

    const double& r = State_Vector[e_r];
    const double& theta = State_Vector[e_theta];

    const double& sin_theta = sin(theta);
    const double& cos_theta = cos(theta);
    const double& N = 1 - this->Parameters.Horizon_radius / r;
    const double& dr_N = this->Parameters.Horizon_radius / r / r;

    /* -------------------------------------------------------------- Evaluate the metric potentials spline -------------------------------------------------------------- */

    Spline_arguments_type s_Spline_args{};

    s_Spline_args.r_coord = r;
    s_Spline_args.r_coord_compactified = compactify_radial_coordiante(s_Spline_args.r_coord);
    s_Spline_args.Radial_idx = Radial_grid_idx;
    s_Spline_args.Theta_idx = Theta_grid_idx;

    s_Spline_args.Radial_natural_param = (s_Spline_args.r_coord_compactified - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1])
                                       / (this->Parameters.Compactified_radial_grid[Radial_grid_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1]);

    s_Spline_args.Theta_natural_param = (State_Vector[e_theta] - this->Parameters.Theta_grid[Theta_grid_idx - 1])
                                      / (this->Parameters.Theta_grid[Theta_grid_idx] - this->Parameters.Theta_grid[Theta_grid_idx - 1]);

    const Numerical_metric_potentials_type s_Potentials = this->compute_metric_components_from_spline(s_Spline_args, None);
    const Numerical_metric_potentials_type s_dr_Potentials = this->compute_metric_components_from_spline(s_Spline_args, First_radial_derivative);
    const Numerical_metric_potentials_type s_dtheta_Potentials = this->compute_metric_components_from_spline(s_Spline_args, First_theta_derivative);
    const Numerical_metric_potentials_type s_dr_dtheta_Potentials = this->compute_metric_components_from_spline(s_Spline_args, Second_mixed_derivative);

    /* ------------------------------------------------------------------------------------------------------------------------------------------------------------------- */

    /* ------ References for the sake of readability ------ */

    const double& W = s_Potentials.W;
    const double exp_2F_0 = exp(2 * s_Potentials.F_0);
    const double exp_2F_1 = exp(2 * s_Potentials.F_1);
    const double exp_2F_2 = exp(2 * s_Potentials.F_2);

    const double& dtheta_W   = s_dtheta_Potentials.W;
    const double& dtheta_F_0 = s_dtheta_Potentials.F_0;
    const double& dtheta_F_1 = s_dtheta_Potentials.F_1;
    const double& dtheta_F_2 = s_dtheta_Potentials.F_2;

    const double& dr_W   = s_dr_Potentials.W;
    const double& dr_F_0 = s_dr_Potentials.F_0;
    const double& dr_F_1 = s_dr_Potentials.F_1;
    const double& dr_F_2 = s_dr_Potentials.F_2;

    const double& dr_dtheta_W = s_dr_dtheta_Potentials.W;
    const double& dr_dtheta_F_0 = s_dr_dtheta_Potentials.F_0;
    const double& dr_dtheta_F_1 = s_dr_dtheta_Potentials.F_1;
    const double& dr_dtheta_F_2 = s_dr_dtheta_Potentials.F_2;

    /* ------------------------------------------------------ */

    Metric_type s_dr_dtheta_Metric{};

    switch (this->Parameters.e_Anzatz) {

    case e_Anzatz_2:

        /* TODO: fix these */

        s_dr_dtheta_Metric.Metric[e_t][e_t] = 0;

        s_dr_dtheta_Metric.Metric[e_t][e_phi] = 0;

        break;


    default:

        s_dr_dtheta_Metric.Metric[e_t][e_t] = -2 * exp_2F_0 * (dr_N * dtheta_F_0 + 2 * N * dr_F_0 * dtheta_F_0 + N * dr_dtheta_F_0)
                                            + (4 * dtheta_F_2 * dr_F_2 * W * W + 4 * W * dtheta_W * dr_F_2 + 2 * dr_dtheta_F_2 * W * W + 4 * dtheta_F_2 * dr_W * W + 2 * dtheta_W * dr_W + 2 * dr_dtheta_W * W) * exp_2F_2 * sin_theta * sin_theta
                                            + (4 * dr_F_2 * W * W + 4 * W * dr_W) * exp_2F_2 * sin_theta * cos_theta;

        s_dr_dtheta_Metric.Metric[e_t][e_phi] = -2 * (2 * dr_F_2 * W * r + dr_W * r + W) * (sin_theta * cos_theta + dtheta_F_2 * sin_theta * sin_theta) * exp_2F_2 
                                              - (2 * dr_dtheta_F_2 * W * r + 2 * dr_F_2 * dtheta_W * r + dr_dtheta_W * r + dtheta_W) * exp_2F_2 * sin_theta * sin_theta;

        break;

    }

    s_dr_dtheta_Metric.Metric[e_phi][e_t] = s_dr_dtheta_Metric.Metric[e_t][e_phi];
    s_dr_dtheta_Metric.Metric[e_r][e_r] = 2 * exp_2F_1 / N * dtheta_F_1 * (2 * dr_F_1 - dr_N / N) + 2 * exp_2F_1 / N * dr_dtheta_F_1;
    s_dr_dtheta_Metric.Metric[e_theta][e_theta] = 2 * exp_2F_1 * (2 * dr_F_1 * dtheta_F_1 * r * r + dr_dtheta_F_1 * r * r + 2 * dtheta_F_1 * r);
    s_dr_dtheta_Metric.Metric[e_phi][e_phi] = 2 * exp_2F_2 * (2 * dr_F_2 * dtheta_F_2 * r * r + dr_dtheta_F_2 * r * r + 2 * dtheta_F_2 * r) * sin_theta * sin_theta
                                            + 4 * exp_2F_2 * (r * r * dr_F_2 + r) * sin_theta * cos_theta;



    return s_dr_dtheta_Metric;

}

void Numerical_metric::get_EOM_Jacobian(const double* const State_Vector, double Jacobian[e_Dynamic_state_size][e_Dynamic_state_size]) {

    /* The Jacobian for the EOM in an axially symmetric spacetime in vacuum has the following form:
        |---------------------------------------------------------------------------------------------------------------------------------------------|
        |0,          d_r g^{t a}p_a,             d_\theta g^{t a}p_a,       0,  |   g^{tt},         0,                      0,             g^{t\phi}  |
        |0,          d_r g^{r a}p_a,             d_\theta g^{r a}p_a,       0,  |     0,         g^{rr},                    0,                 0      |
        |0,        d_r g^{\theta a}p_a,       d_\theta g^{\theta a}p_a,     0,  |     0,            0,               g^{\theta\theta},         0      |
        |0,         d_r g^{\phi a}p_a,         d_\theta g^{\phi a}p_a,      0,  | g^{t\phi},        0,                      0,            g^{\phi\phi}|
        |-----------------------------------------------------------------------|---------------------------------------------------------------------|
        |0,                0,                              0,               0,  |     0,     -d_r g^{t a}p_a,      -d_\theta g^{t a}p_a,       0,     |
        |0,     -1/2 d^2_r g^{ab}p_ap_b,   -1/2 d^2_{r\theta} g^{ab}p_ap_b, 0,  |     0,     -d_r g^{r a}p_a,      -d_\theta g^{r a}p_a,       0,     |
        |0, -1/2 d^2_{r\theta} g^{ab}p_ap_b, -1/2 d^2_\theta g^{ab}p_ap_b,  0,  |     0,   -d_r g^{\theta a}p_a, -d_\theta g^{\theta a}p_a,    0,     |
        |0,                0,                              0,               0,  |     0,    -d_r g^{\phi a}p_a,   -d_\theta g^{\phi a}p_a,     0,     |
        |---------------------------------------------------------------------------------------------------------------------------------------------| */

    const double r_compactified = this->compactify_radial_coordiante(State_Vector[e_r]);

    const int Radial_grid_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size, r_compactified) - this->Parameters.Compactified_radial_grid;
    const int Theta_grid_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size, State_Vector[e_theta]) - this->Parameters.Theta_grid;

    memset(Jacobian, 0, sizeof(double) * 64);

    const Metric_type s_Metric = this->get_metric(State_Vector, Radial_grid_idx, Theta_grid_idx);
    const Metric_type s_dr_Metric = this->get_dr_metric(State_Vector, Radial_grid_idx, Theta_grid_idx);
    const Metric_type s_dtheta_Metric = this->get_dtheta_metric(State_Vector, Radial_grid_idx, Theta_grid_idx);
    const Metric_type s_d2r_Metric = this->get_d2r_metric(State_Vector, Radial_grid_idx, Theta_grid_idx);
    const Metric_type s_d2theta_Metric = this->get_d2theta_metric(State_Vector, Radial_grid_idx, Theta_grid_idx);
    const Metric_type s_dr_dtheta_Metric = this->get_dr_dtheta_metric(State_Vector, Radial_grid_idx, Theta_grid_idx);

    /* ------- Compute the derivatives in the inverse metric ------- */

    double inv_Metric[4][4]{};
    invert_metric(inv_Metric, s_Metric.Metric);

    double dr_inv_metric[4][4]{};
    get_derivative_of_inverse_metric(inv_Metric, s_dr_Metric.Metric, dr_inv_metric);

    double dtheta_inv_metric[4][4]{};
    get_derivative_of_inverse_metric(inv_Metric, s_dtheta_Metric.Metric, dtheta_inv_metric);

    double d2r_inv_metric[4][4]{};
    get_second_derivative_of_inverse_metric(inv_Metric, dr_inv_metric, s_dr_Metric.Metric, s_d2r_Metric.Metric, d2r_inv_metric);

    double d2theta_inv_metric[4][4]{};
    get_second_derivative_of_inverse_metric(inv_Metric, dtheta_inv_metric, s_dtheta_Metric.Metric, s_d2theta_Metric.Metric, d2theta_inv_metric);

    double dr_dtheta_inv_metric[4][4]{};
    get_second_mixed_derivative_of_inverse_metric(inv_Metric, dtheta_inv_metric, s_dr_Metric.Metric, s_dtheta_Metric.Metric, s_dr_dtheta_Metric.Metric, dr_dtheta_inv_metric);

    /* -------------------------- Compute the contractions with the photon momentum --------------------------*/

    for (int idx = 0; idx < 4; idx++) {

        Jacobian[0][1] += dr_inv_metric[e_t][idx] * State_Vector[idx + e_p_t];
        Jacobian[1][1] += dr_inv_metric[e_r][idx] * State_Vector[idx + e_p_t];
        Jacobian[2][1] += dr_inv_metric[e_theta][idx] * State_Vector[idx + e_p_t];
        Jacobian[3][1] += dr_inv_metric[e_phi][idx] * State_Vector[idx + e_p_t];

        Jacobian[0][2] += dtheta_inv_metric[e_t][idx] * State_Vector[idx + e_p_t];
        Jacobian[1][2] += dtheta_inv_metric[e_r][idx] * State_Vector[idx + e_p_t];
        Jacobian[2][2] += dtheta_inv_metric[e_theta][idx] * State_Vector[idx + e_p_t];
        Jacobian[3][2] += dtheta_inv_metric[e_phi][idx] * State_Vector[idx + e_p_t];

        Jacobian[4][5] -= dr_inv_metric[e_t][idx] * State_Vector[idx + e_p_t];
        Jacobian[5][5] -= dr_inv_metric[e_r][idx] * State_Vector[idx + e_p_t];
        Jacobian[6][5] -= dr_inv_metric[e_theta][idx] * State_Vector[idx + e_p_t];
        Jacobian[7][5] -= dr_inv_metric[e_phi][idx] * State_Vector[idx + e_p_t];

        Jacobian[4][6] -= dtheta_inv_metric[e_t][idx] * State_Vector[idx + e_p_t];
        Jacobian[5][6] -= dtheta_inv_metric[e_r][idx] * State_Vector[idx + e_p_t];
        Jacobian[6][6] -= dtheta_inv_metric[e_theta][idx] * State_Vector[idx + e_p_t];
        Jacobian[7][6] -= dtheta_inv_metric[e_phi][idx] * State_Vector[idx + e_p_t];

        for (int idx_2 = 0; idx_2 < 4; idx_2++) {

            Jacobian[5][1] -= 0.5 * d2r_inv_metric[idx][idx_2] * State_Vector[idx + e_p_t] * State_Vector[idx_2 + e_p_t];
            Jacobian[6][1] -= 0.5 * dr_dtheta_inv_metric[idx][idx_2] * State_Vector[idx + e_p_t] * State_Vector[idx_2 + e_p_t];
            Jacobian[5][2] -= 0.5 * dr_dtheta_inv_metric[idx][idx_2] * State_Vector[idx + e_p_t] * State_Vector[idx_2 + e_p_t];
            Jacobian[6][2] -= 0.5 * d2theta_inv_metric[idx][idx_2] * State_Vector[idx + e_p_t] * State_Vector[idx_2 + e_p_t];

        }

    }

    Jacobian[0][4] = inv_Metric[e_t][e_t];
    Jacobian[0][7] = inv_Metric[e_t][e_phi];
    Jacobian[1][5] = inv_Metric[e_r][e_r];
    Jacobian[2][6] = inv_Metric[e_theta][e_theta];
    Jacobian[3][7] = inv_Metric[e_phi][e_phi];
    Jacobian[3][4] = inv_Metric[e_t][e_phi];

}

std::complex<double> Numerical_metric::get_largest_EOM_eigenvalue(const double* const State_Vector) {

    double EOM_Jacobian[e_Dynamic_state_size][e_Dynamic_state_size]{};

    this->get_EOM_Jacobian(State_Vector, EOM_Jacobian);

    for (int idx_1 = 0; idx_1 < e_Dynamic_state_size; idx_1++) {

        for (int idx_2 = 0; idx_2 < e_Dynamic_state_size; idx_2++) {

            gsl_matrix_set(this->EOM_Jacobian_gsl_matrix, idx_1, idx_2, EOM_Jacobian[idx_1][idx_2]);

        }
    }

    gsl_eigen_nonsymmv(this->EOM_Jacobian_gsl_matrix, this->EOM_Eigenvalues, this->EOM_Eigenvectors, this->GSL_workspace);

    std::complex<double> Largest_eigenvalue{};

    for (int idx = 0; idx < 8; idx++) {
    
        double Current_eigenvalue_norm = gsl_complex_abs(gsl_vector_complex_get(this->EOM_Eigenvalues, idx));

        if (Current_eigenvalue_norm > std::abs(Largest_eigenvalue)) {

            Largest_eigenvalue.real(GSL_REAL(gsl_vector_complex_get(this->EOM_Eigenvalues, idx)));
            Largest_eigenvalue.imag(GSL_IMAG(gsl_vector_complex_get(this->EOM_Eigenvalues, idx)));

        }

    }

    return Largest_eigenvalue;

}

int Numerical_metric::get_initial_conditions_from_file(Initial_conditions_type* p_Initial_Conditions, double J_data[], double p_theta_data[], int photon) {

    return 0;
}

void Numerical_metric::get_EOM(const double* const State_Vector, double* const Derivatives) const {

    const double r_compactified = this->compactify_radial_coordiante(State_Vector[e_r]);
    const int Radial_grid_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size, r_compactified) - this->Parameters.Compactified_radial_grid;
    const int Theta_grid_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size, State_Vector[e_theta]) - this->Parameters.Theta_grid;

    const Metric_type Metric = this->get_metric(State_Vector, Radial_grid_idx, Theta_grid_idx);
    const Metric_type dr_Metric = this->get_dr_metric(State_Vector, Radial_grid_idx, Theta_grid_idx);
    const Metric_type dtheta_Metric = this->get_dtheta_metric(State_Vector, Radial_grid_idx, Theta_grid_idx);

    /* ----------- Temporary matrix, used to store intermediate calculations ----------- */
    double temp_matrix[4][4]{};

    double inv_metric[4][4]{};
    invert_metric(inv_metric, Metric.Metric);

    double dr_inv_metric[4][4]{};
    matrix_matrix_multiply(dr_Metric.Metric, inv_metric, temp_matrix);
    matrix_matrix_multiply(inv_metric, temp_matrix, dr_inv_metric);
    
    double dtheta_inv_metric[4][4]{};
    matrix_matrix_multiply(dtheta_Metric.Metric, inv_metric, temp_matrix);
    matrix_matrix_multiply(inv_metric, temp_matrix, dtheta_inv_metric);

    /* ----------- There is a minus sign infront of the whole expression for the derivative of an inverse of a matrix ----------- */
    for (int left_idx = 0; left_idx < 4; left_idx++) {

        for (int right_idx = 0; right_idx < 4; right_idx++) {

            dr_inv_metric[left_idx][right_idx] *= -1;
            dtheta_inv_metric[left_idx][right_idx] *= -1;

        }

    }

    for (int right_idx = 0; right_idx < 4; right_idx++) {

        Derivatives[e_t] += inv_metric[e_t][right_idx] * State_Vector[right_idx + 4];
        Derivatives[e_r] += inv_metric[e_r][right_idx] * State_Vector[right_idx + 4];
        Derivatives[e_theta] += inv_metric[e_theta][right_idx] * State_Vector[right_idx + 4];
        Derivatives[e_phi] += inv_metric[e_phi][right_idx] * State_Vector[right_idx + 4];

        for (int left_idx = 0; left_idx < 4; left_idx++) {

            Derivatives[e_p_r] += -1. / 2 * dr_inv_metric[left_idx][right_idx] * State_Vector[left_idx + 4] * State_Vector[right_idx + 4];
            Derivatives[e_p_theta] += -1. / 2 * dtheta_inv_metric[left_idx][right_idx] * State_Vector[left_idx + 4] * State_Vector[right_idx + 4];

        }

    }

    Derivatives[e_p_t] = 0.0;
    Derivatives[e_p_phi] = 0.0;

}

bool Numerical_metric::terminate_integration(const double* const State_vector) {

    const bool scatter = State_vector[e_r] > this->Scattering_radius && State_vector[e_p_r] < 0;

    const bool hit_horizon = State_vector[e_r] - this->Parameters.Horizon_radius < this->Min_distance_to_singular_point;

    return scatter || hit_horizon;
};
