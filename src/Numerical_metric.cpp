#include "Spacetimes.h"

Numerical_metric::Numerical_metric(const Metric_parameters_type* const p_Metric_Parameters) {

    if (isnan(p_Metric_Parameters->Scattering_radius) or isinf(p_Metric_Parameters->Scattering_radius) or p_Metric_Parameters->Scattering_radius < 0) {

        throw std::runtime_error(std::format("Invalid value for the scattering radius: {}", p_Metric_Parameters->Scattering_radius));

    }

    if (isnan(p_Metric_Parameters->Min_distance_to_singular_point) or isinf(p_Metric_Parameters->Min_distance_to_singular_point)) {

        throw std::runtime_error(std::format("Invalid value for the distance to the throat: {}", p_Metric_Parameters->Min_distance_to_singular_point));

    }

    this->Min_distance_to_singular_point = p_Metric_Parameters->Min_distance_to_singular_point;
    this->Scattering_radius = p_Metric_Parameters->Scattering_radius;

    this->Parameters = p_Metric_Parameters->Numerical_metric_params;

    if (this->Parameters.e_Spline_type != Spline_selection_enums::Custom_cubic) {

        if (this->Parameters.e_Spline_type == Spline_selection_enums::GSL_linear) {

            this->Spline_instance_exp_2F_0 = gsl_spline2d_alloc(gsl_interp2d_bilinear, this->Parameters.Radial_grid_size, this->Parameters.Theta_grid_size);
            this->Spline_instance_exp_2F_1 = gsl_spline2d_alloc(gsl_interp2d_bilinear, this->Parameters.Radial_grid_size, this->Parameters.Theta_grid_size);
            this->Spline_instance_exp_2F_2 = gsl_spline2d_alloc(gsl_interp2d_bilinear, this->Parameters.Radial_grid_size, this->Parameters.Theta_grid_size);
            this->Spline_instance_W = gsl_spline2d_alloc(gsl_interp2d_bilinear, this->Parameters.Radial_grid_size, this->Parameters.Theta_grid_size);

        }
        else {

            this->Spline_instance_exp_2F_0 = gsl_spline2d_alloc(gsl_interp2d_bicubic, this->Parameters.Radial_grid_size, this->Parameters.Theta_grid_size);
            this->Spline_instance_exp_2F_1 = gsl_spline2d_alloc(gsl_interp2d_bicubic, this->Parameters.Radial_grid_size, this->Parameters.Theta_grid_size);
            this->Spline_instance_exp_2F_2 = gsl_spline2d_alloc(gsl_interp2d_bicubic, this->Parameters.Radial_grid_size, this->Parameters.Theta_grid_size);
            this->Spline_instance_W = gsl_spline2d_alloc(gsl_interp2d_bicubic, this->Parameters.Radial_grid_size, this->Parameters.Theta_grid_size);

        }

        this->Radial_interp_accelerator = gsl_interp_accel_alloc();
        this->Theta_interp_accelerator = gsl_interp_accel_alloc();

        gsl_spline2d_init(this->Spline_instance_exp_2F_0,
                          this->Parameters.Compactified_radial_grid, 
                          this->Parameters.Theta_grid, 
                          this->Parameters.Raw_F_0_data, 
                          this->Parameters.Radial_grid_size, 
                          this->Parameters.Theta_grid_size);

        gsl_spline2d_init(this->Spline_instance_exp_2F_1,
                          this->Parameters.Compactified_radial_grid,
                          this->Parameters.Theta_grid,
                          this->Parameters.Raw_F_1_data,
                          this->Parameters.Radial_grid_size,
                          this->Parameters.Theta_grid_size);

        gsl_spline2d_init(this->Spline_instance_exp_2F_2,
                          this->Parameters.Compactified_radial_grid,
                          this->Parameters.Theta_grid,
                          this->Parameters.Raw_F_2_data,
                          this->Parameters.Radial_grid_size,
                          this->Parameters.Theta_grid_size);

        gsl_spline2d_init(this->Spline_instance_W,
                          this->Parameters.Compactified_radial_grid,
                          this->Parameters.Theta_grid,
                          this->Parameters.Raw_W_data,
                          this->Parameters.Radial_grid_size,
                          this->Parameters.Theta_grid_size);

    }
    else {

        this->Spline_instance_exp_2F_0 = nullptr;
        this->Spline_instance_exp_2F_1 = nullptr;
        this->Spline_instance_exp_2F_2 = nullptr;
        this->Spline_instance_W = nullptr;

        this->Radial_interp_accelerator = nullptr;
        this->Theta_interp_accelerator = nullptr;

    }

}

Numerical_metric::~Numerical_metric() {

    gsl_spline2d_free(this->Spline_instance_exp_2F_0);
    gsl_spline2d_free(this->Spline_instance_exp_2F_1);
    gsl_spline2d_free(this->Spline_instance_exp_2F_2);
    gsl_spline2d_free(this->Spline_instance_W);

    gsl_interp_accel_free(this->Radial_interp_accelerator);
    gsl_interp_accel_free(this->Theta_interp_accelerator);

    delete this->Parameters.Raw_F_0_data;
    delete this->Parameters.Raw_F_1_data;
    delete this->Parameters.Raw_F_2_data;
    delete this->Parameters.Raw_W_data;

    delete this->Parameters.Compactified_radial_grid;
    delete this->Parameters.Theta_grid;

}

double Numerical_metric::compactify_radial_coordiante(const double r) const {

    /* Comute the shifted coordinate that the paper (http://gravitation.web.ua.pt/node/416) uses in the numerical implementation (they label this little "x"). */
    const double x_uncompactified = sqrt(r * r - this->Parameters.Horizon_radius * this->Parameters.Horizon_radius);

    return x_uncompactified / (1 + x_uncompactified);
}

inline Delta_coeffs_type Numerical_metric::get_delta_matrix(const double* const Grid_step_array, const long long idx) const {

    Delta_coeffs_type Coeff_struct{};

    Coeff_struct.a_coeff = Grid_step_array[idx] * Grid_step_array[idx] / ((Grid_step_array[idx - 2] + Grid_step_array[idx - 1] + Grid_step_array[idx]) * (Grid_step_array[idx - 1] + Grid_step_array[idx]));

    Coeff_struct.b_coeff = Grid_step_array[idx] * Grid_step_array[idx] / ((Grid_step_array[idx - 1] + Grid_step_array[idx] + Grid_step_array[idx + 1]) * (Grid_step_array[idx - 1] + Grid_step_array[idx]));

    Coeff_struct.c_coeff = Grid_step_array[idx] * Grid_step_array[idx] / ((Grid_step_array[idx - 1] + Grid_step_array[idx] + Grid_step_array[idx + 1]) * (Grid_step_array[idx] + Grid_step_array[idx + 1]));

    Coeff_struct.d_coeff = Grid_step_array[idx] * Grid_step_array[idx] / ((Grid_step_array[idx] + Grid_step_array[idx + 1] + Grid_step_array[idx + 2]) * (Grid_step_array[idx] + Grid_step_array[idx + 1]));

    Coeff_struct.e_coeff = Grid_step_array[idx] * Grid_step_array[idx - 1] / ((Grid_step_array[idx - 1] + Grid_step_array[idx] + Grid_step_array[idx + 1]) * (Grid_step_array[idx - 1] + Grid_step_array[idx]));

    Coeff_struct.f_coeff = Grid_step_array[idx - 1] * Grid_step_array[idx - 1] / ((Grid_step_array[idx - 1] + Grid_step_array[idx] + Grid_step_array[idx + 1]) * (Grid_step_array[idx - 1] + Grid_step_array[idx]));

    return Coeff_struct;

}

inline void Numerical_metric::get_control_point_matrix(const double* const Control_vector, const long long r_idx, const long long theta_idx, double Control_matrix[4][4]) const {

    /* The reference for this implementation is one of the python example scrips in this document: https://hal.science/hal-03017566/document. */

    for (long long radial_offset = 0; radial_offset < 4; radial_offset++) {

        for (long long theta_offset = 0; theta_offset < 4; theta_offset++) {

            Control_matrix[radial_offset][theta_offset] = Control_vector[r_idx + radial_offset + theta_offset * (this->Parameters.Radial_grid_size + 2) + theta_idx * (this->Parameters.Radial_grid_size + 2)];
            Control_matrix[radial_offset][theta_offset] = Control_vector[r_idx + radial_offset + theta_offset * (this->Parameters.Radial_grid_size + 2) + theta_idx * (this->Parameters.Radial_grid_size + 2)];
            Control_matrix[radial_offset][theta_offset] = Control_vector[r_idx + radial_offset + theta_offset * (this->Parameters.Radial_grid_size + 2) + theta_idx * (this->Parameters.Radial_grid_size + 2)];
            Control_matrix[radial_offset][theta_offset] = Control_vector[r_idx + radial_offset + theta_offset * (this->Parameters.Radial_grid_size + 2) + theta_idx * (this->Parameters.Radial_grid_size + 2)];

        }

    }

}

void Numerical_metric::get_polynomial_basis_vector(const double natural_parameter, const Delta_coeffs_type* const Delta_coeffs, double* const Polynomial_basis_vector) const {

    Polynomial_basis_vector[0] = -Delta_coeffs->a_coeff * natural_parameter * natural_parameter * natural_parameter
                               + 3 * Delta_coeffs->a_coeff * natural_parameter * natural_parameter
                               - 3 * Delta_coeffs->a_coeff * natural_parameter
                               + Delta_coeffs->a_coeff;

    Polynomial_basis_vector[1] = (Delta_coeffs->a_coeff + Delta_coeffs->b_coeff + Delta_coeffs->c_coeff) * natural_parameter * natural_parameter * natural_parameter
                               + (-3 * Delta_coeffs->a_coeff - 3 * Delta_coeffs->b_coeff) * natural_parameter * natural_parameter
                               + ( 3 * Delta_coeffs->a_coeff - 3 * Delta_coeffs->e_coeff) * natural_parameter
                               + 1 - Delta_coeffs->a_coeff - Delta_coeffs->f_coeff;

    Polynomial_basis_vector[2] = (-Delta_coeffs->b_coeff - Delta_coeffs->c_coeff - Delta_coeffs->d_coeff) * natural_parameter * natural_parameter * natural_parameter
                               + 3 * Delta_coeffs->b_coeff * natural_parameter * natural_parameter
                               + 3 * Delta_coeffs->e_coeff * natural_parameter 
                               + Delta_coeffs->f_coeff;

    Polynomial_basis_vector[3] = Delta_coeffs->d_coeff * natural_parameter * natural_parameter * natural_parameter;

}

void Numerical_metric::get_derivative_polynomial_basis_vector(const double natural_parameter, const Delta_coeffs_type* const Delta_coeffs, double* const Polynomial_basis_vector) const {

    Polynomial_basis_vector[0] = -3 * Delta_coeffs->a_coeff * natural_parameter * natural_parameter
                               + 6 * Delta_coeffs->a_coeff * natural_parameter
                               - 3 * Delta_coeffs->a_coeff;

    Polynomial_basis_vector[1] = 3 * (Delta_coeffs->a_coeff + Delta_coeffs->b_coeff + Delta_coeffs->c_coeff) * natural_parameter * natural_parameter
                               + 2 * (-3 * Delta_coeffs->a_coeff - 3 * Delta_coeffs->b_coeff) * natural_parameter
                               + (3 * Delta_coeffs->a_coeff - 3 * Delta_coeffs->e_coeff);

    Polynomial_basis_vector[2] = 3 * (-Delta_coeffs->b_coeff - Delta_coeffs->c_coeff - Delta_coeffs->d_coeff) * natural_parameter * natural_parameter
                                + 6 * Delta_coeffs->b_coeff * natural_parameter
                                + 3 * Delta_coeffs->e_coeff;

    Polynomial_basis_vector[3] = 3 * Delta_coeffs->d_coeff * natural_parameter * natural_parameter;

}

void Numerical_metric::get_second_derivative_polynomial_basis_vector(const double natural_parameter, const Delta_coeffs_type* const Delta_coeffs, double* const Polynomial_basis_vector) const {

    Polynomial_basis_vector[0] = -6 * Delta_coeffs->a_coeff * natural_parameter
                                + 6 * Delta_coeffs->a_coeff;

    Polynomial_basis_vector[1] = 6 * (Delta_coeffs->a_coeff + Delta_coeffs->b_coeff + Delta_coeffs->c_coeff) * natural_parameter
                               + 2 * (-3 * Delta_coeffs->a_coeff - 3 * Delta_coeffs->b_coeff);

    Polynomial_basis_vector[2] = 6 * (-Delta_coeffs->b_coeff - Delta_coeffs->c_coeff - Delta_coeffs->d_coeff) * natural_parameter
                               + 6 * Delta_coeffs->b_coeff;

    Polynomial_basis_vector[3] = 6 * Delta_coeffs->d_coeff * natural_parameter;

}

double Numerical_metric::evaluate_single_spline(const double Control_point_matrix[4][4], 
                                                const double Radial_natural_parameter, 
                                                const double Theta_natural_parameter,
                                                const Delta_coeffs_type* const Radial_coeffs,
                                                const Delta_coeffs_type* const Theta_coeffs,
                                                Derivative_selector_enums Derivative_selector) const {

    /* -------------- Compute the basis polynomials ------------ */

    double Radial_basis_polynomial[4]{};
    double Theta_basis_polynomial[4]{};

    switch (Derivative_selector) {

    case None:

        this->get_polynomial_basis_vector(Radial_natural_parameter, Radial_coeffs, Radial_basis_polynomial);
        this->get_polynomial_basis_vector(Theta_natural_parameter, Theta_coeffs, Theta_basis_polynomial);

        break;

    case First_radial_derivative:

        this->get_derivative_polynomial_basis_vector(Radial_natural_parameter, Radial_coeffs, Radial_basis_polynomial);
        this->get_polynomial_basis_vector(Theta_natural_parameter, Theta_coeffs, Theta_basis_polynomial);

        break;

    case First_theta_derivative:

        this->get_polynomial_basis_vector(Radial_natural_parameter, Radial_coeffs, Radial_basis_polynomial);
        this->get_derivative_polynomial_basis_vector(Theta_natural_parameter, Theta_coeffs, Theta_basis_polynomial);

        break;

    case Second_radial_derivative:

        this->get_second_derivative_polynomial_basis_vector(Radial_natural_parameter, Radial_coeffs, Radial_basis_polynomial);
        this->get_polynomial_basis_vector(Theta_natural_parameter, Theta_coeffs, Theta_basis_polynomial);

        break;

    case Second_mixed_derivative:

        this->get_derivative_polynomial_basis_vector(Radial_natural_parameter, Radial_coeffs, Radial_basis_polynomial);
        this->get_derivative_polynomial_basis_vector(Theta_natural_parameter, Theta_coeffs, Theta_basis_polynomial);

        break;

    case Second_theta_derivative:

        this->get_polynomial_basis_vector(Radial_natural_parameter, Radial_coeffs, Radial_basis_polynomial);
        this->get_second_derivative_polynomial_basis_vector(Theta_natural_parameter, Theta_coeffs, Theta_basis_polynomial);

        break;

    default: throw std::runtime_error("Unsupported numerical metric derivative! \n");

    }

    /* -------------- Compute the double dot product between the two polynomial basis vectors and the control point matrix. --------------*/
    // NOTE: The control point matrix acts on the "theta" basis polynomial on the right, and on the "radial" one on the left.

    double Intermediate_result[4]{};
    mat_vec_multiply_4D(Control_point_matrix, Theta_basis_polynomial, Intermediate_result);

    return dot_product(Radial_basis_polynomial, Intermediate_result, 4);

}

Numerical_metric_potentials_type Numerical_metric::evaluate_all_custom_splines(Spline_arguments_type s_Spline_args, Derivative_selector_enums Derivative_selector) const {

    Numerical_metric_potentials_type s_Metric_potentials{};

    // The +3 offset comes from the padding of the steps array at the begining
    const Delta_coeffs_type Radial_coeffs = this->get_delta_matrix(this->Parameters.Compactified_radial_grid_steps, s_Spline_args.Radial_idx - 1 + 3);
    const Delta_coeffs_type Theta_coeffs = this->get_delta_matrix(this->Parameters.Theta_grid_steps, s_Spline_args.Theta_idx - 1 + 3);

    /* -------------- Compute F_0 -------------- */

    double F_0_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.F_0_control_vector, s_Spline_args.Radial_idx - 1, s_Spline_args.Theta_idx - 1, F_0_control_point_matrix);

    s_Metric_potentials.exp_2F_0 = this->evaluate_single_spline(F_0_control_point_matrix, s_Spline_args.Radial_natural_param, s_Spline_args.Theta_natural_param, &Radial_coeffs, &Theta_coeffs, Derivative_selector);

    /* -------------- Compute F_1 -------------- */

    double F_1_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.F_1_control_vector, s_Spline_args.Radial_idx - 1, s_Spline_args.Theta_idx - 1, F_1_control_point_matrix);

    s_Metric_potentials.exp_2F_1 = this->evaluate_single_spline(F_1_control_point_matrix, s_Spline_args.Radial_natural_param, s_Spline_args.Theta_natural_param, &Radial_coeffs, &Theta_coeffs, Derivative_selector);

    /* -------------- Compute F_2 -------------- */

    double F_2_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.F_2_control_vector, s_Spline_args.Radial_idx - 1, s_Spline_args.Theta_idx - 1, F_2_control_point_matrix);

    s_Metric_potentials.exp_2F_2 = this->evaluate_single_spline(F_2_control_point_matrix, s_Spline_args.Radial_natural_param, s_Spline_args.Theta_natural_param, &Radial_coeffs, &Theta_coeffs, Derivative_selector);

    /* -------------- Compute W -------------- */

    double W_control_point_matrix[4][4]{};
    this->get_control_point_matrix(this->Parameters.W_control_vector, s_Spline_args.Radial_idx - 1, s_Spline_args.Theta_idx - 1, W_control_point_matrix);

    s_Metric_potentials.W = this->evaluate_single_spline(W_control_point_matrix, s_Spline_args.Radial_natural_param, s_Spline_args.Theta_natural_param, &Radial_coeffs, &Theta_coeffs, Derivative_selector);

    return s_Metric_potentials;

}

Numerical_metric_potentials_type Numerical_metric::evaluate_all_gsl_splines(Spline_arguments_type s_Spline_args, Derivative_selector_enums Derivative_selector) const {

    // Hack-y error handling so I dont get a "y value out of range" GSL error
    // TODO: Make this less hacky

    if (s_Spline_args.Theta_coord > M_PI) { s_Spline_args.Theta_coord = M_PI; }

    if (s_Spline_args.Theta_coord < 0) { s_Spline_args.Theta_coord = 0; }

    Numerical_metric_potentials_type s_Metric_potentials{};

    switch (Derivative_selector) {

    case Derivative_selector_enums::None:

        s_Metric_potentials.W = gsl_spline2d_eval(this->Spline_instance_W, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_0 = gsl_spline2d_eval(this->Spline_instance_exp_2F_0, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_1 = gsl_spline2d_eval(this->Spline_instance_exp_2F_1, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_2 = gsl_spline2d_eval(this->Spline_instance_exp_2F_2, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);

        break;

    case Derivative_selector_enums::First_radial_derivative:

        s_Metric_potentials.W = gsl_spline2d_eval_deriv_x(this->Spline_instance_W, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_0 = gsl_spline2d_eval_deriv_x(this->Spline_instance_exp_2F_0, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_1 = gsl_spline2d_eval_deriv_x(this->Spline_instance_exp_2F_1, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_2 = gsl_spline2d_eval_deriv_x(this->Spline_instance_exp_2F_2, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);

        break;

    case Derivative_selector_enums::First_theta_derivative:

        s_Metric_potentials.W = gsl_spline2d_eval_deriv_y(this->Spline_instance_W, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_0 = gsl_spline2d_eval_deriv_y(this->Spline_instance_exp_2F_0, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_1 = gsl_spline2d_eval_deriv_y(this->Spline_instance_exp_2F_1, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_2 = gsl_spline2d_eval_deriv_y(this->Spline_instance_exp_2F_2, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);

        break;

    case Derivative_selector_enums::Second_radial_derivative:

        s_Metric_potentials.W = gsl_spline2d_eval_deriv_xx(this->Spline_instance_W, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_0 = gsl_spline2d_eval_deriv_xx(this->Spline_instance_exp_2F_0, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_1 = gsl_spline2d_eval_deriv_xx(this->Spline_instance_exp_2F_1, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_2 = gsl_spline2d_eval_deriv_xx(this->Spline_instance_exp_2F_2, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);

        break;

    case Derivative_selector_enums::Second_mixed_derivative:

        s_Metric_potentials.W = gsl_spline2d_eval_deriv_xy(this->Spline_instance_W, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_0 = gsl_spline2d_eval_deriv_xy(this->Spline_instance_exp_2F_0, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_1 = gsl_spline2d_eval_deriv_xy(this->Spline_instance_exp_2F_1, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_2 = gsl_spline2d_eval_deriv_xy(this->Spline_instance_exp_2F_2, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);

        break;

    case Derivative_selector_enums::Second_theta_derivative:

        s_Metric_potentials.W = gsl_spline2d_eval_deriv_yy(this->Spline_instance_W, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_0 = gsl_spline2d_eval_deriv_yy(this->Spline_instance_exp_2F_0, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_1 = gsl_spline2d_eval_deriv_yy(this->Spline_instance_exp_2F_1, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);
        s_Metric_potentials.exp_2F_2 = gsl_spline2d_eval_deriv_yy(this->Spline_instance_exp_2F_2, s_Spline_args.r_coord_compactified, s_Spline_args.Theta_coord, this->Radial_interp_accelerator, this->Theta_interp_accelerator);

        break;

    }

    return s_Metric_potentials;

}

Numerical_metric_potentials_type Numerical_metric::evaluate_all_splines(Spline_arguments_type s_Spline_args, Derivative_selector_enums Derivative_selector) const {

    if (this->Parameters.e_Spline_type == Spline_selection_enums::Custom_cubic) {

        return this->evaluate_all_custom_splines(s_Spline_args, Derivative_selector);

    }
    else {

        return this->evaluate_all_gsl_splines(s_Spline_args, Derivative_selector);

    }

}

Numerical_metric_potentials_type Numerical_metric::compute_metric_components_from_spline(Spline_arguments_type s_Spline_args, Derivative_selector_enums Derivative_selector) const {

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

        Corrected_Potentials = this->evaluate_all_splines(s_Spline_args, Derivative_selector);

        if (this->Parameters.e_Spline_type == Spline_selection_enums::Custom_cubic) {

            /* ------- This is correcting by the factor d(natural_parameter)/d(x_compactified) ------- */
            derivative_correction_factor_1 = 1 / (this->Parameters.Compactified_radial_grid[s_Spline_args.Radial_idx] - this->Parameters.Compactified_radial_grid[s_Spline_args.Radial_idx - 1]);

        }

        /* ------- This is correcting by the factor d(x_compactified)/d(x_uncompactified) ------- */
        derivative_correction_factor_1 *= (1 - s_Spline_args.r_coord_compactified) * (1 - s_Spline_args.r_coord_compactified);

        /* ------- This is correcting by the factor d(x_uncompactified)/d(r) ------- */
        derivative_correction_factor_1 *= s_Spline_args.r_coord * (1 - s_Spline_args.r_coord_compactified) / s_Spline_args.r_coord_compactified;

        Corrected_Potentials.exp_2F_0 *= derivative_correction_factor_1;
        Corrected_Potentials.exp_2F_1 *= derivative_correction_factor_1;
        Corrected_Potentials.exp_2F_2 *= derivative_correction_factor_1;
        Corrected_Potentials.W   *= derivative_correction_factor_1;

        return Corrected_Potentials;

    case First_theta_derivative:

        Corrected_Potentials = this->evaluate_all_splines(s_Spline_args, Derivative_selector);

        if (this->Parameters.e_Spline_type != Spline_selection_enums::Custom_cubic) {  return Corrected_Potentials; }

        /* ------- This is correcting by the factor d(natural_parameter)/d(x_compactified) ------- */
        derivative_correction_factor_1 = 1 / (this->Parameters.Theta_grid[s_Spline_args.Theta_idx] - this->Parameters.Theta_grid[s_Spline_args.Theta_idx - 1]);

        Corrected_Potentials.exp_2F_0 *= derivative_correction_factor_1;
        Corrected_Potentials.exp_2F_1 *= derivative_correction_factor_1;
        Corrected_Potentials.exp_2F_2 *= derivative_correction_factor_1;
        Corrected_Potentials.W   *= derivative_correction_factor_1;

        return Corrected_Potentials;

    case Second_radial_derivative:

        temp_Potentials_1 = this->evaluate_all_splines(s_Spline_args, Second_radial_derivative);
        temp_Potentials_2 = this->evaluate_all_splines(s_Spline_args, First_radial_derivative);
        temp_Potentials_3 = temp_Potentials_2;

        if (this->Parameters.e_Spline_type == Spline_selection_enums::Custom_cubic) {

            /* ------- This is correcting by the factor d(natural_parameter)/d(x_compactified) ------- */
            derivative_correction_factor_1 = 1 / (this->Parameters.Compactified_radial_grid[s_Spline_args.Radial_idx] - this->Parameters.Compactified_radial_grid[s_Spline_args.Radial_idx - 1]);

        }

        /* ------- This is correcting by the factor d(x_compactified)/d(x_uncompactified) ------- */
        derivative_correction_factor_2 = (1 - s_Spline_args.r_coord_compactified) * (1 - s_Spline_args.r_coord_compactified);

        /* ------- This is correcting by the factor d(x_uncompactified)/d(r) ------- */
        derivative_correction_factor_3 = s_Spline_args.r_coord * (1 - s_Spline_args.r_coord_compactified) / s_Spline_args.r_coord_compactified;

        /* ------- This is correcting by the factor d^2(dx_compactified)/d(x_uncompactified)^2 ------- */
        derivative_correction_factor_4 = -2 * (1 - s_Spline_args.r_coord_compactified) * (1 - s_Spline_args.r_coord_compactified) * (1 - s_Spline_args.r_coord_compactified);

        /* ------- This is correcting by the factor d^2(dx_uncompactified)/d(r)^2 ------- */
        derivative_correction_factor_5 = -s_Spline_args.r_coord * s_Spline_args.r_coord * (1 - s_Spline_args.r_coord_compactified) * (1 - s_Spline_args.r_coord_compactified) * (1 - s_Spline_args.r_coord_compactified) / s_Spline_args.r_coord_compactified / s_Spline_args.r_coord_compactified / s_Spline_args.r_coord_compactified;
        derivative_correction_factor_5 += (1 - s_Spline_args.r_coord_compactified) / s_Spline_args.r_coord_compactified;

        temp_Potentials_1.exp_2F_0 *= (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3) * (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3);
        temp_Potentials_1.exp_2F_1 *= (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3) * (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3);
        temp_Potentials_1.exp_2F_2 *= (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3) * (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3);
        temp_Potentials_1.W   *= (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3) * (derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_3);

        temp_Potentials_2.exp_2F_0 *= derivative_correction_factor_1 * derivative_correction_factor_4 * derivative_correction_factor_3 * derivative_correction_factor_3;
        temp_Potentials_2.exp_2F_1 *= derivative_correction_factor_1 * derivative_correction_factor_4 * derivative_correction_factor_3 * derivative_correction_factor_3;
        temp_Potentials_2.exp_2F_2 *= derivative_correction_factor_1 * derivative_correction_factor_4 * derivative_correction_factor_3 * derivative_correction_factor_3;
        temp_Potentials_2.W   *= derivative_correction_factor_1 * derivative_correction_factor_4 * derivative_correction_factor_3 * derivative_correction_factor_3;

        temp_Potentials_3.exp_2F_0 *= derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_5;
        temp_Potentials_3.exp_2F_1 *= derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_5;
        temp_Potentials_3.exp_2F_2 *= derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_5;
        temp_Potentials_3.W   *= derivative_correction_factor_1 * derivative_correction_factor_2 * derivative_correction_factor_5;

        Corrected_Potentials.exp_2F_0 = temp_Potentials_1.exp_2F_0 + temp_Potentials_2.exp_2F_0 + temp_Potentials_3.exp_2F_0;
        Corrected_Potentials.exp_2F_1 = temp_Potentials_1.exp_2F_1 + temp_Potentials_2.exp_2F_1 + temp_Potentials_3.exp_2F_1;
        Corrected_Potentials.exp_2F_2 = temp_Potentials_1.exp_2F_2 + temp_Potentials_2.exp_2F_2 + temp_Potentials_3.exp_2F_2;
        Corrected_Potentials.W   = temp_Potentials_1.W  + temp_Potentials_2.W  + temp_Potentials_3.W;  

        return Corrected_Potentials;

    default:

        return this->evaluate_all_splines(s_Spline_args, Derivative_selector);

    }

}

Metric_type Numerical_metric::get_local_metric(const double* const Local_State_Vector) const {

    /* ---------------- This is a wrapper function for compatability with the radiative transfer part of the code ---------------- */

    double r_compactified{};
    long long int Radial_grid_upper_idx{}, Theta_grid_upper_idx{};

    if (this->Parameters.e_Spline_type == Spline_selection_enums::Custom_cubic) {

        r_compactified = this->compactify_radial_coordiante(Local_State_Vector[e_r]);

        Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size - 1, r_compactified) - this->Parameters.Compactified_radial_grid;
        Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size - 1, Local_State_Vector[e_theta]) - this->Parameters.Theta_grid;
    }

    return this->get_metric(Local_State_Vector, std::as_const(Radial_grid_upper_idx), std::as_const(Theta_grid_upper_idx));

}

Metric_type Numerical_metric::get_global_metric(const double* const Global_State_Vector) const {

    return this->get_local_metric(Global_State_Vector);

}

Metric_type Numerical_metric::get_metric(const double* const Local_State_Vector, long long Radial_grid_idx, long long Theta_grid_idx) const {

    /* -------- The metric antatz is from https://arxiv.org/pdf/1501.04319. -------- */

    /* -------------------------------------------------------------- Evaluate the metric potentials spline -------------------------------------------------------------- */

    Spline_arguments_type s_Spline_args{};

    s_Spline_args.r_coord = Local_State_Vector[e_r];
    s_Spline_args.Theta_coord = Local_State_Vector[e_theta];

    s_Spline_args.r_coord_compactified = compactify_radial_coordiante(Local_State_Vector[e_r]);

    s_Spline_args.Radial_idx = Radial_grid_idx;
    s_Spline_args.Theta_idx = Theta_grid_idx;

    s_Spline_args.Radial_natural_param = (s_Spline_args.r_coord_compactified - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1]) 
                                      / (this->Parameters.Compactified_radial_grid[Radial_grid_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1]);

    s_Spline_args.Theta_natural_param = (Local_State_Vector[e_theta] - this->Parameters.Theta_grid[Theta_grid_idx - 1])
                                      / (this->Parameters.Theta_grid[Theta_grid_idx] - this->Parameters.Theta_grid[Theta_grid_idx - 1]);

    const Numerical_metric_potentials_type s_Potentials = this->compute_metric_components_from_spline(s_Spline_args, None);

    /* ------------------------------------------------------------------------------------------------------------------------------------------------------------------- */

    /* ---- References for the sake of readability ---- */

    const double& W = s_Potentials.W;
    const double& exp_2F_0 = s_Potentials.exp_2F_0;
    const double& exp_2F_1 = s_Potentials.exp_2F_1;
    const double& exp_2F_2 = s_Potentials.exp_2F_2;

    const double sin_theta = sin(Local_State_Vector[e_theta]);
    const double N = 1 - this->Parameters.Horizon_radius / Local_State_Vector[e_r];

    const double& r = Local_State_Vector[e_r];

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

Metric_type Numerical_metric::get_dr_local_metric(const double* const Local_State_Vector) const {

    double r_compactified{};
    long long int Radial_grid_upper_idx{}, Theta_grid_upper_idx{};

    if (this->Parameters.e_Spline_type == Spline_selection_enums::Custom_cubic) {

        r_compactified = this->compactify_radial_coordiante(Local_State_Vector[e_r]);

        Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size - 1, r_compactified) - this->Parameters.Compactified_radial_grid;
        Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size - 1, Local_State_Vector[e_theta]) - this->Parameters.Theta_grid;

    }

    return this->get_dr_metric(Local_State_Vector, Radial_grid_upper_idx, Theta_grid_upper_idx);

}

Metric_type Numerical_metric::get_dr_global_metric(const double* const Global_State_Vector) const {

    return this->get_dr_local_metric(Global_State_Vector);

}

Metric_type Numerical_metric::get_dr_metric(const double* const Local_State_Vector, long long Radial_grid_idx, long long Theta_grid_idx) const {

    /* -------- The metric antatz is from https://arxiv.org/pdf/1501.04319. */

    /* -------------------------------------------------------------- Evaluate the metric potentials spline -------------------------------------------------------------- */

    Spline_arguments_type s_Spline_args{};

    s_Spline_args.r_coord = Local_State_Vector[e_r];
    s_Spline_args.Theta_coord = Local_State_Vector[e_theta];

    s_Spline_args.r_coord_compactified = compactify_radial_coordiante(Local_State_Vector[e_r]);

    s_Spline_args.Radial_idx = Radial_grid_idx;
    s_Spline_args.Theta_idx = Theta_grid_idx;

    s_Spline_args.Radial_natural_param = (s_Spline_args.r_coord_compactified - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1])
                                       / (this->Parameters.Compactified_radial_grid[Radial_grid_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1]);

    s_Spline_args.Theta_natural_param = (Local_State_Vector[e_theta] - this->Parameters.Theta_grid[Theta_grid_idx - 1])
                                      / (this->Parameters.Theta_grid[Theta_grid_idx] - this->Parameters.Theta_grid[Theta_grid_idx - 1]);

    const Numerical_metric_potentials_type s_Potentials = this->compute_metric_components_from_spline(s_Spline_args, None);
    const Numerical_metric_potentials_type s_dr_Potentials = this->compute_metric_components_from_spline(s_Spline_args, First_radial_derivative);

    /* ------------------------------------------------------------------------------------------------------------------------------------------------------------------- */

    /* ---- References for the sake of readability ---- */

    const double& W = s_Potentials.W;
    const double& exp_2F_0 = s_Potentials.exp_2F_0;
    const double& exp_2F_1 = s_Potentials.exp_2F_1;
    const double& exp_2F_2 = s_Potentials.exp_2F_2;

    const double& dr_W = s_dr_Potentials.W;
    const double& dr_exp_2F_0 = s_dr_Potentials.exp_2F_0;
    const double& dr_exp_2F_1 = s_dr_Potentials.exp_2F_1;
    const double& dr_exp_2F_2 = s_dr_Potentials.exp_2F_2;

    const double& r = Local_State_Vector[e_r];

    const double sin_theta = sin(Local_State_Vector[e_theta]);
    const double N = 1 - this->Parameters.Horizon_radius / r;
    const double dr_N = this->Parameters.Horizon_radius / r / r;

    /* ------------------------------------------------ */

    Metric_type s_dr_Metric{};

    switch (this->Parameters.e_Anzatz) {

    case e_Anzatz_2:

        s_dr_Metric.Metric[e_t][e_t] = -dr_N * exp_2F_0 - N * dr_exp_2F_0 + (dr_exp_2F_2 * W * W * r * r + 2 * exp_2F_2 * W * dr_W * r * r + 2 * exp_2F_2 * W * W * r) * sin_theta * sin_theta;
        s_dr_Metric.Metric[e_t][e_phi] = -(dr_exp_2F_2 * W * r * r + exp_2F_2 * dr_W * r * r + 2 * exp_2F_2 * W * r) * sin_theta * sin_theta;

        break;

    default:

        s_dr_Metric.Metric[e_t][e_t] = -dr_N * exp_2F_0 - N * dr_exp_2F_0 + dr_exp_2F_2 * W * W * sin_theta * sin_theta + 2 * exp_2F_2 * W * dr_W * sin_theta * sin_theta;
        s_dr_Metric.Metric[e_t][e_phi] = -(dr_exp_2F_2 * W * r + exp_2F_2 * W + exp_2F_2 * r * dr_W) * sin_theta * sin_theta;

        break;

    }

    s_dr_Metric.Metric[e_phi][e_t] = s_dr_Metric.Metric[e_t][e_phi];
    s_dr_Metric.Metric[e_r][e_r] = dr_exp_2F_1 / N  - dr_N / (N * N) * exp_2F_1;
    s_dr_Metric.Metric[e_theta][e_theta] = 2 * r * exp_2F_1 + r * r * dr_exp_2F_1;
    s_dr_Metric.Metric[e_phi][e_phi] = (2 * r * exp_2F_2 + r * r * dr_exp_2F_2) * sin_theta * sin_theta;

    return s_dr_Metric;

}

Metric_type Numerical_metric::get_dtheta_local_metric(const double* const Local_State_Vector) const {

    /* ---------------- This is a wrapper function for compatability with the radiative transfer part of the code ---------------- */

    double r_compactified{};
    long long int Radial_grid_upper_idx{}, Theta_grid_upper_idx{};

    if (this->Parameters.e_Spline_type == Spline_selection_enums::Custom_cubic) {

        r_compactified = this->compactify_radial_coordiante(Local_State_Vector[e_r]);

        Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size - 1, r_compactified) - this->Parameters.Compactified_radial_grid;
        Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size - 1, Local_State_Vector[e_theta]) - this->Parameters.Theta_grid;

    }

    return this->get_dtheta_metric(Local_State_Vector, Radial_grid_upper_idx, Theta_grid_upper_idx);

}

Metric_type Numerical_metric::get_dtheta_global_metric(const double* const Global_State_Vector) const {

    return this->get_dtheta_local_metric(Global_State_Vector);

}

Metric_type Numerical_metric::get_dtheta_metric(const double* const Local_State_Vector, long long Radial_grid_idx, long long Theta_grid_idx) const {

    /* -------- The metric antatz is from https://arxiv.org/pdf/1501.04319. */

    /* -------------------------------------------------------------- Evaluate the metric potentials spline -------------------------------------------------------------- */

    Spline_arguments_type s_Spline_args{};

    s_Spline_args.r_coord = Local_State_Vector[e_r];
    s_Spline_args.Theta_coord = Local_State_Vector[e_theta];

    s_Spline_args.r_coord_compactified = compactify_radial_coordiante(Local_State_Vector[e_r]);

    s_Spline_args.Radial_idx = Radial_grid_idx;
    s_Spline_args.Theta_idx = Theta_grid_idx;

    s_Spline_args.Radial_natural_param = (s_Spline_args.r_coord_compactified - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1])
                                       / (this->Parameters.Compactified_radial_grid[Radial_grid_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1]);

    s_Spline_args.Theta_natural_param = (Local_State_Vector[e_theta] - this->Parameters.Theta_grid[Theta_grid_idx - 1])
                                      / (this->Parameters.Theta_grid[Theta_grid_idx] - this->Parameters.Theta_grid[Theta_grid_idx - 1]);

    const Numerical_metric_potentials_type s_Potentials = this->compute_metric_components_from_spline(s_Spline_args, None);
    const Numerical_metric_potentials_type s_dtheta_Potentials = this->compute_metric_components_from_spline(s_Spline_args, First_theta_derivative);

    /* ------------------------------------------------------------------------------------------------------------------------------------------------------------------- */

    /* ---- References for the sake of readability ---- */

    const double& W = s_Potentials.W;
    const double& exp_2F_0 = s_Potentials.exp_2F_0;
    const double& exp_2F_1 = s_Potentials.exp_2F_1;
    const double& exp_2F_2 = s_Potentials.exp_2F_2;

    const double& dtheta_W = s_dtheta_Potentials.W;
    const double& dtheta_exp_2F_0 = s_dtheta_Potentials.exp_2F_0;
    const double& dtheta_exp_2F_1 = s_dtheta_Potentials.exp_2F_1;
    const double& dtheta_exp_2F_2 = s_dtheta_Potentials.exp_2F_2;

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    const double sin_theta = sin(theta);
    const double cos_theta = cos(theta);
    const double N = 1 - this->Parameters.Horizon_radius / r;

    /* ------------------------------------------------- */

    Metric_type s_dtheta_Metric{};

    switch (this->Parameters.e_Anzatz) {

    case e_Anzatz_2:

        s_dtheta_Metric.Metric[e_t][e_t] = -N * dtheta_exp_2F_0 + (dtheta_exp_2F_2 * W * W * sin_theta * sin_theta + 2 * exp_2F_2 * W * dtheta_W * sin_theta * sin_theta + 2 * exp_2F_2 * W * W * sin_theta * cos_theta) * r * r;
        s_dtheta_Metric.Metric[e_t][e_phi] = -dtheta_exp_2F_2 * W * r * r * sin_theta * sin_theta - exp_2F_2 * dtheta_W * r * r * sin_theta * sin_theta - 2 * exp_2F_2 * W * r * r * sin_theta * cos_theta;

        break;

    default:

        s_dtheta_Metric.Metric[e_t][e_t] = -N * dtheta_exp_2F_0 + dtheta_exp_2F_2 * W * W * sin_theta * sin_theta + 2 * exp_2F_2 * W * dtheta_W * sin_theta * sin_theta + 2 * exp_2F_2 * W * W * sin_theta * cos_theta;
        s_dtheta_Metric.Metric[e_t][e_phi] = -dtheta_exp_2F_2 * W * r * sin_theta * sin_theta - exp_2F_2 * dtheta_W * r * sin_theta * sin_theta - 2 * exp_2F_2 * W * r * sin_theta * cos_theta;

        break;

    }

    s_dtheta_Metric.Metric[e_phi][e_t] = s_dtheta_Metric.Metric[e_t][e_phi];
    s_dtheta_Metric.Metric[e_r][e_r] = dtheta_exp_2F_1 / N;
    s_dtheta_Metric.Metric[e_theta][e_theta] = r * r * dtheta_exp_2F_1;
    s_dtheta_Metric.Metric[e_phi][e_phi] = r * r * dtheta_exp_2F_2 * sin_theta * sin_theta + 2 * r * r * exp_2F_2 * sin_theta * cos_theta;

    return s_dtheta_Metric;

}

Metric_type Numerical_metric::get_d2r_local_metric(const double* const Local_State_Vector) const {

    /* ---------------- This is a wrapper function for compatability with the radiative transfer part of the code ---------------- */

    const double r_compactified = this->compactify_radial_coordiante(Local_State_Vector[e_r]);

    const auto Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size - 1, r_compactified) - this->Parameters.Compactified_radial_grid;
    const auto Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size - 1, Local_State_Vector[e_theta]) - this->Parameters.Theta_grid;

    return this->get_d2r_metric(Local_State_Vector, Radial_grid_upper_idx, Theta_grid_upper_idx);

}

Metric_type Numerical_metric::get_d2r_metric(const double* const Local_State_Vector, long long Radial_grid_idx, long long Theta_grid_idx) const {

    /* -------- The metric antatz is from https://arxiv.org/pdf/1501.04319. */

    /* -------------------------------------------------------------- Evaluate the metric potentials spline -------------------------------------------------------------- */

    Spline_arguments_type s_Spline_args{};

    s_Spline_args.r_coord = Local_State_Vector[e_r];
    s_Spline_args.Theta_coord = Local_State_Vector[e_theta];

    s_Spline_args.r_coord_compactified = compactify_radial_coordiante(Local_State_Vector[e_r]);

    s_Spline_args.Radial_idx = Radial_grid_idx;
    s_Spline_args.Theta_idx = Theta_grid_idx;

    s_Spline_args.Radial_natural_param = (s_Spline_args.r_coord_compactified - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1])
                                       / (this->Parameters.Compactified_radial_grid[Radial_grid_idx] - this->Parameters.Compactified_radial_grid[Radial_grid_idx - 1]);

    s_Spline_args.Theta_natural_param = (Local_State_Vector[e_theta] - this->Parameters.Theta_grid[Theta_grid_idx - 1])
                                      / (this->Parameters.Theta_grid[Theta_grid_idx] - this->Parameters.Theta_grid[Theta_grid_idx - 1]);

    const Numerical_metric_potentials_type s_Potentials = this->compute_metric_components_from_spline(s_Spline_args, None);
    const Numerical_metric_potentials_type s_dr_Potentials = this->compute_metric_components_from_spline(s_Spline_args, First_radial_derivative);
    const Numerical_metric_potentials_type s_d2r_Potentials = this->compute_metric_components_from_spline(s_Spline_args, Second_radial_derivative);

    /* ------------------------------------------------------------------------------------------------------------------------------------------------------------------- */

    /* ---- References for the sake of readability ---- */

    const double& W = s_Potentials.W;
    const double& exp_2F_0 = s_Potentials.exp_2F_0;
    const double& exp_2F_1 = s_Potentials.exp_2F_1;
    const double& exp_2F_2 = s_Potentials.exp_2F_2;

    const double& dr_W = s_dr_Potentials.W;
    const double& dr_exp_2F_0 = s_dr_Potentials.exp_2F_0;
    const double& dr_exp_2F_1 = s_dr_Potentials.exp_2F_1;
    const double& dr_exp_2F_2 = s_dr_Potentials.exp_2F_2;

    const double& d2r_W = s_d2r_Potentials.W;
    const double& d2r_exp_2F_0 = s_d2r_Potentials.exp_2F_0;
    const double& d2r_exp_2F_1 = s_d2r_Potentials.exp_2F_1;
    const double& d2r_exp_2F_2 = s_d2r_Potentials.exp_2F_2;

    const double& r = Local_State_Vector[e_r];
    const double& theta = Local_State_Vector[e_theta];

    const double sin_theta = sin(theta);
    const double N = 1 - this->Parameters.Horizon_radius / r;
    const double dr_N = this->Parameters.Horizon_radius / r / r;
    const double d2r_N = -2 * this->Parameters.Horizon_radius / r / r / r;

    /* ------------------------------------------------ */

    Metric_type s_d2r_Metric{};

    switch (this->Parameters.e_Anzatz) {

    case e_Anzatz_2:

        s_d2r_Metric.Metric[e_t][e_t] = -d2r_N * exp_2F_0 - 2 * dr_N * dr_exp_2F_0 - N * d2r_exp_2F_0
                                      + (d2r_exp_2F_2 * W * W * r * r + 4 * dr_exp_2F_2 * W * dr_W * r * r + 4 * dr_exp_2F_2 * W * W * r 
                                          + 2 * exp_2F_2 * dr_W * dr_W * r * r + 2 * exp_2F_2 * W * d2r_W * r * r + 8 * exp_2F_2 * W * dr_W * r + 2 * exp_2F_2 * W * W) * sin_theta * sin_theta;

        s_d2r_Metric.Metric[e_t][e_phi] = -(d2r_exp_2F_2 * W * r * r + 2 * dr_exp_2F_2 * dr_W * r * r + 4 * dr_exp_2F_2 * W * r
                                            + exp_2F_2 * d2r_W * r * r + 4 * exp_2F_2 * dr_W * r + 2 * exp_2F_2 * W) * sin_theta * sin_theta;

        break;

    default:

        s_d2r_Metric.Metric[e_t][e_t] = -d2r_N * exp_2F_0 - 2 * dr_N * dr_exp_2F_0 - N * d2r_exp_2F_0 
                                      + (d2r_exp_2F_2 * W * W + 4 * dr_exp_2F_2 * W * dr_W + 2 * exp_2F_2 * dr_W * dr_W + 2 * exp_2F_2 * W * d2r_W) * sin_theta * sin_theta;
        s_d2r_Metric.Metric[e_t][e_phi] = -(d2r_exp_2F_2 * W * r + 2 * dr_exp_2F_2 * dr_W * r + 2 * dr_exp_2F_2 * W + 2 * exp_2F_2 * dr_W + exp_2F_2 * r * d2r_W) * sin_theta * sin_theta;

        break;

    }

    s_d2r_Metric.Metric[e_phi][e_t] = s_d2r_Metric.Metric[e_t][e_phi];
    s_d2r_Metric.Metric[e_r][e_r] = d2r_exp_2F_1 / N - 2 * dr_exp_2F_1 * dr_N / (N * N) + 2 * exp_2F_1 * dr_N * dr_N / (N * N * N) - exp_2F_1 * d2r_N / (N * N);
    s_d2r_Metric.Metric[e_theta][e_theta] = 2 * exp_2F_1 + 4 * r * dr_exp_2F_1 + r * r * d2r_exp_2F_1;
    s_d2r_Metric.Metric[e_phi][e_phi] = (2 * exp_2F_2 + 4 * r * dr_exp_2F_2 + r * r * d2r_exp_2F_2) * sin_theta * sin_theta;

    return s_d2r_Metric;
}

void Numerical_metric::get_EOM(const double* const State_Vector, double* const Derivatives) {

    memset(Derivatives, 0, e_Dynamic_state_size * sizeof(double));

    double r_compactified{};
    long long int Radial_grid_upper_idx{}, Theta_grid_upper_idx{};

    if (this->Parameters.e_Spline_type == Spline_selection_enums::Custom_cubic) {

        r_compactified = this->compactify_radial_coordiante(State_Vector[e_r]);

        Radial_grid_upper_idx = std::upper_bound(this->Parameters.Compactified_radial_grid, this->Parameters.Compactified_radial_grid + this->Parameters.Radial_grid_size - 1, r_compactified) - this->Parameters.Compactified_radial_grid;
        Theta_grid_upper_idx = std::upper_bound(this->Parameters.Theta_grid, this->Parameters.Theta_grid + this->Parameters.Theta_grid_size - 1, State_Vector[e_theta]) - this->Parameters.Theta_grid;

    }

    const Metric_type Metric = this->get_metric(State_Vector, Radial_grid_upper_idx, Theta_grid_upper_idx);
    const Metric_type dr_Metric = this->get_dr_metric(State_Vector, Radial_grid_upper_idx, Theta_grid_upper_idx);
    const Metric_type dtheta_Metric = this->get_dtheta_metric(State_Vector, Radial_grid_upper_idx, Theta_grid_upper_idx);

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

    const bool scatter = State_vector[e_r] > this->Scattering_radius and State_vector[e_p_r] < 0;

    const bool hit_horizon = State_vector[e_r] - this->Parameters.Horizon_radius < this->Min_distance_to_singular_point;

    return scatter or hit_horizon;
};

void Numerical_metric::Convert_global_to_local_coords(const double* const State_Vector_Global, const double* const Global_Vec_to_Convert, double* Local_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

    switch (Entry_to_convert) {

    case e_Full_State_Vector:

        memcpy(Local_Vec_to_Convert, Global_Vec_to_Convert, e_Full_state_size * sizeof(double));
        break;

    case e_Contravariant_vector:

        memcpy(Local_Vec_to_Convert, Global_Vec_to_Convert, 4 * sizeof(double));
        break;

    case e_Covariant_vector:
        memcpy(Local_Vec_to_Convert, Global_Vec_to_Convert, 4 * sizeof(double));
        break;

    case e_Coordinates:
        memcpy(Local_Vec_to_Convert, Global_Vec_to_Convert, 4 * sizeof(double));
        break;

    default:

        throw std::runtime_error("Unsupported coordinate conversion type. Something Broke in Convert_global_to_local_coords()!");

    }

}

void Numerical_metric::Convert_local_to_global_coords(const double* const State_Vector_Local, const double* const Local_Vec_to_Convert, double* Global_Vec_to_Convert, Coord_conversion_enums Entry_to_convert) {

    switch (Entry_to_convert) {

    case e_Full_State_Vector:

        memcpy(Global_Vec_to_Convert, Local_Vec_to_Convert, e_Full_state_size * sizeof(double));
        break;

    case e_Contravariant_vector:

        memcpy(Global_Vec_to_Convert, Local_Vec_to_Convert, 4 * sizeof(double));
        break;

    case e_Covariant_vector:
        memcpy(Global_Vec_to_Convert, Local_Vec_to_Convert, 4 * sizeof(double));
        break;

    case e_Coordinates:
        memcpy(Global_Vec_to_Convert, Local_Vec_to_Convert, 4 * sizeof(double));
        break;

    default:

        throw std::runtime_error("Unsupported coordinate conversion type. Something Broke in Convert_local_to_global_coords()!");

    }

}