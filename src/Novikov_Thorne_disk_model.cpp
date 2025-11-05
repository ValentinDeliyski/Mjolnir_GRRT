#include "Novikov_Thorne_Model.h"

Novikov_Thorne_Model_class::Novikov_Thorne_Model_class(Simulation_Context_type* p_Sim_Context) {

    this->r_in  = p_Sim_Context->p_Init_Conditions->Disk_params.Novikov_Thorne_params.r_in;
    this->r_out = p_Sim_Context->p_Init_Conditions->Disk_params.Novikov_Thorne_params.r_out;

    this->flux_integral_accuracy = p_Sim_Context->p_Init_Conditions->Integrator_params.Simpson_accuracy;

    this->p_Spacetime = p_Sim_Context->p_Spacetime;
    this->e_Spacetime = p_Sim_Context->p_Init_Conditions->Metric_parameters.e_Spacetime;

    this->e_Mag_field_geometry = p_Sim_Context->p_Init_Conditions->Disk_params.e_Mag_field_geometry;

    this->current_flux_integration_step = 0;
    this->max_flux_integration_teps = 500;

    memset(this->Disk_veclovity_vector, 0, 4 * sizeof(double));
    memset(this->Source_polarization_vector, 0, 4 * sizeof(double));
    memcpy(this->Mag_field_geometry, p_Sim_Context->p_Init_Conditions->Disk_params.Mag_field_geometry, 3 * sizeof(double));

}

double Novikov_Thorne_Model_class::Keplerian_angular_velocity(const double* const State_Vector) {

    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_metric(State_Vector);

    double Angular_velocity = (-s_dr_Metric.Metric[e_t][e_phi] + sqrt(s_dr_Metric.Metric[e_t][e_phi] * s_dr_Metric.Metric[e_t][e_phi] - s_dr_Metric.Metric[e_t][e_t] * s_dr_Metric.Metric[e_phi][e_phi])) / s_dr_Metric.Metric[e_phi][e_phi];

    if (isnan(Angular_velocity) || isinf(Angular_velocity) || isnan(Angular_velocity) || isinf(Angular_velocity)) {

        throw std::runtime_error(std::format("Invalid Novikov-Thorne angular velocity at: r = {}: Omega = {}", State_Vector[e_r], Angular_velocity));

    }

    return Angular_velocity;

}

double Novikov_Thorne_Model_class::dr_Keplerian_angular_velocity(const double* const State_Vector) {

    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_metric(State_Vector);
    Metric_type s_d2r_Metric = this->p_Spacetime->get_d2r_metric(State_Vector);

    double root = sqrt(s_dr_Metric.Metric[e_t][e_phi] * s_dr_Metric.Metric[e_t][e_phi] - s_dr_Metric.Metric[e_t][e_t] * s_dr_Metric.Metric[e_phi][e_phi]);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

    double dr_Kepler = -Kepler / s_dr_Metric.Metric[e_phi][e_phi] * s_d2r_Metric.Metric[e_phi][e_phi] + (-s_d2r_Metric.Metric[e_t][e_phi]
                     + 1.0 / root / 2 * (2 * s_dr_Metric.Metric[e_t][e_phi] * s_d2r_Metric.Metric[e_t][e_phi] - s_dr_Metric.Metric[e_t][e_t] * s_d2r_Metric.Metric[e_phi][e_phi]
                     - s_d2r_Metric.Metric[e_t][e_t] * s_dr_Metric.Metric[e_phi][e_phi])) / s_dr_Metric.Metric[e_phi][e_phi];

    if (isnan(dr_Kepler) || isinf(dr_Kepler) || isnan(dr_Kepler) || isinf(dr_Kepler)) {

        throw std::runtime_error(std::format("Invalid Novikov-Thorne angular velocity derivative at: r = {}: dr_Omega = {}", State_Vector[e_r], dr_Kepler));

    }

    return dr_Kepler;
}

double* Novikov_Thorne_Model_class::get_disk_velocity_vector(const double* const State_Vector) {

    const double& r_source = State_Vector[e_r];
    const double& theta_source = State_Vector[e_theta];

    Metric_type s_Metric_source = this->p_Spacetime->get_metric(State_Vector);

    double Angular_velocity = this->Keplerian_angular_velocity(State_Vector);
    double u_t = 1 / sqrt(-s_Metric_source.Metric[e_t][e_t] - 2 * s_Metric_source.Metric[e_t][e_phi] * Angular_velocity - s_Metric_source.Metric[e_phi][e_phi] * Angular_velocity * Angular_velocity);

    if (isnan(u_t) || isinf(u_t) || isnan(Angular_velocity) || isinf(Angular_velocity)) {

        throw std::runtime_error(std::format("Invalid Novikov-Thorne disk 4-velocity at r = {}: u_t = {}, Angular velocity = {}", State_Vector[e_r], u_t, Angular_velocity));

    }

    this->Disk_veclovity_vector[e_t] = u_t;
    this->Disk_veclovity_vector[e_r] = 0.0;
    this->Disk_veclovity_vector[e_theta] = 0.0;
    this->Disk_veclovity_vector[e_phi] = u_t * Angular_velocity;

    return this->Disk_veclovity_vector;

}

double Novikov_Thorne_Model_class::disk_Energy(const double* const State_Vector) {

    Metric_type s_Metric_source = this->p_Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

    double root = sqrt(-s_Metric_source.Metric[e_t][e_t] - 2 * s_Metric_source.Metric[e_t][e_phi] * Kepler - s_Metric_source.Metric[e_phi][e_phi] * Kepler * Kepler);
    double Disk_Energy = -(s_Metric_source.Metric[e_t][e_t] + s_Metric_source.Metric[e_t][e_phi] * Kepler) / root;

    if (isnan(Disk_Energy) || isinf(Disk_Energy) || isnan(Disk_Energy) || isinf(Disk_Energy)) {

        throw std::runtime_error(std::format("Invalid Novikov-Thorne disk energy at r = {}: E = {}", State_Vector[e_r], Disk_Energy));

    }

    return Disk_Energy;

}

double Novikov_Thorne_Model_class::disk_Angular_Momentum(const double* const State_Vector) {

    Metric_type s_Metric_source = this->p_Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

    double root = sqrt(-s_Metric_source.Metric[e_t][e_t] - 2 * s_Metric_source.Metric[e_t][e_phi] * Kepler - s_Metric_source.Metric[e_phi][e_phi] * Kepler * Kepler);
    double Disk_angular_momentum = (s_Metric_source.Metric[e_phi][e_phi] * Kepler + s_Metric_source.Metric[e_t][e_phi]) / root;

    if (isnan(Disk_angular_momentum) || isinf(Disk_angular_momentum) || isnan(Disk_angular_momentum) || isinf(Disk_angular_momentum)) {

        throw std::runtime_error(std::format("Invalid Novikov-Thorne disk angular momentum at r = {}: L_z = {}", State_Vector[e_r], Disk_angular_momentum));

    }

    return  Disk_angular_momentum;

}

double Novikov_Thorne_Model_class::Flux_integrand(const double* const State_Vector) {

    Metric_type s_Metric = this->p_Spacetime->get_metric(State_Vector);
    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);
    double dr_Kepler = this->dr_Keplerian_angular_velocity(State_Vector);

    double root = sqrt(-s_Metric.Metric[e_t][e_t] - 2 * s_Metric.Metric[e_t][e_phi] * Kepler - s_Metric.Metric[e_phi][e_phi] * Kepler * Kepler);
    double dr_root = (-s_dr_Metric.Metric[e_t][e_t] - 2 * (s_dr_Metric.Metric[e_t][e_phi] * Kepler + s_Metric.Metric[e_t][e_phi] * dr_Kepler)
        - s_dr_Metric.Metric[e_phi][e_phi] * Kepler * Kepler - 2 * s_Metric.Metric[e_phi][e_phi] * Kepler * dr_Kepler);

    double E = this->disk_Energy(State_Vector);
    double L = this->disk_Angular_Momentum(State_Vector);

    double dr_L = (s_dr_Metric.Metric[e_phi][e_phi] * Kepler + s_Metric.Metric[e_phi][e_phi] * dr_Kepler + s_dr_Metric.Metric[e_t][e_phi]) / root - L / root / root / 2 * dr_root;

    double Flux_integrand = (E - Kepler * L) * dr_L;

    if (isnan(Flux_integrand) || isinf(Flux_integrand) || isnan(Flux_integrand) || isinf(Flux_integrand)) {

        throw std::runtime_error(std::format("Invalid Novikov-Thorne disk flux integrand at r = {}: Integrand = {}", State_Vector[e_r], Flux_integrand));

    }

    return Flux_integrand;

}

double Novikov_Thorne_Model_class::solve_Flux_integral(double r_in, const double* const State_Vector, double tolerance) {

    if (this->current_flux_integration_step > this->max_flux_integration_teps) {

        throw std::runtime_error("Novikov-Thorne flux integral not converging!");

    }

    this->current_flux_integration_step += 1;

    const double& lower_bound = r_in;
    const double& upper_bound = State_Vector[e_r];

    double mid_point          = (lower_bound + upper_bound) / 2;
    double left_of_mid_point  = (lower_bound + mid_point) / 2;
    double right_of_mid_point = (mid_point + upper_bound) / 2;

    double lower_bound_state_vector[4]{};
    double mid_point_state_vector[4]{};
    double left_of_mid_point_state_vector[4]{};
    double right_of_mid_point_state_vector[4]{};

    memcpy(lower_bound_state_vector, State_Vector, 4 * sizeof(double));
    memcpy(mid_point_state_vector, State_Vector, 4 * sizeof(double));
    memcpy(left_of_mid_point_state_vector, State_Vector, 4 * sizeof(double));
    memcpy(right_of_mid_point_state_vector, State_Vector, 4 * sizeof(double));

    lower_bound_state_vector[e_r] = r_in;
    mid_point_state_vector[e_r] = mid_point;
    left_of_mid_point_state_vector[e_r] = left_of_mid_point;
    right_of_mid_point_state_vector[e_r] = right_of_mid_point;

    double F_lower_bound = this->Flux_integrand(lower_bound_state_vector);
    double F_mid_point   = this->Flux_integrand(mid_point_state_vector);
    double F_upper_bound = this->Flux_integrand(State_Vector);

    double F_left_mid = this->Flux_integrand(left_of_mid_point_state_vector);
    double F_right_mid = this->Flux_integrand(right_of_mid_point_state_vector);

    double S_left = (mid_point - lower_bound) / 6 * (F_lower_bound + 4 * F_left_mid + F_mid_point);
    double S_right = (upper_bound - mid_point) / 6 * (F_mid_point + 4 * F_right_mid + F_upper_bound);

    double S_2 = S_left + S_right;
    double S_1 = (upper_bound - lower_bound) / 6 * (F_lower_bound + 4 * F_mid_point + F_upper_bound);

    double integral;

    if (fabs(S_2 - S_1) < 15 * tolerance) {

        integral = S_2 + (S_2 - S_1) / 15;


    }
    else {

        double L_value = this->solve_Flux_integral(lower_bound, mid_point_state_vector, tolerance / 2);
        double R_value = this->solve_Flux_integral(mid_point, State_Vector, tolerance / 2);

        integral = L_value + R_value;

    }

    return integral;

}

double Novikov_Thorne_Model_class::get_flux(const double* const State_Vector) {

    Metric_type s_Metric = this->p_Spacetime->get_metric(State_Vector);

    double metric_det = get_eq_induced_metric_det(s_Metric.Metric);
    double E_disk = disk_Energy(State_Vector);
    double L_disk = disk_Angular_Momentum(State_Vector);

    double Kepler = Keplerian_angular_velocity(State_Vector);
    double dr_Kepler = dr_Keplerian_angular_velocity(State_Vector);

    double Flux_coeff = -dr_Kepler / ((E_disk - Kepler * L_disk) * (E_disk - Kepler * L_disk)) / (4 * M_PI * sqrt(-metric_det));

    double Flux_integral = solve_Flux_integral(this->r_in, State_Vector, this->flux_integral_accuracy);

    this->current_flux_integration_step = 0;

    return Flux_coeff * Flux_integral;

}

double* Novikov_Thorne_Model_class::Construct_coord_polarization_vector(const double* const State_Vector) {

    Metric_type s_Metric = this->p_Spacetime->get_metric(State_Vector);

    /* NOTE: This is the covariant momentum */
    const double* const &Photon_coordinate_momentum_covariant = State_Vector + e_p_t;

    double Photon_coordinate_momentum_contravariant[4]{};
    Manipulate_index(&s_Metric, Photon_coordinate_momentum_covariant, Photon_coordinate_momentum_contravariant, Raise_index);

    double Photon_ZAMO_momentum[4]{};
    Contravariant_coord_to_ZAMO(&s_Metric, Photon_coordinate_momentum_contravariant, Photon_ZAMO_momentum);

    double Disk_ZAMO_velocity[4]{};
    Contravariant_coord_to_ZAMO(&s_Metric, this->get_disk_velocity_vector(State_Vector), Disk_ZAMO_velocity);

    double Boost_matrix[4][4]{};  
    get_Lorentz_boost_matrix(Boost_matrix, Disk_ZAMO_velocity, false);

    double Photon_plasma_momentum[4]{};
    mat_vec_multiply_4D(Boost_matrix, Photon_ZAMO_momentum, Photon_plasma_momentum);

    /* =================================== The polarization vector =================================== */

    double Polarization_vector_plasma[4]{};

    // Offset with e_r so I get the spatial components
    cross_product(Photon_plasma_momentum + e_r, this->Mag_field_geometry, Polarization_vector_plasma + e_r);

    for (int idx = e_r; idx <= e_phi; idx++) {

        Polarization_vector_plasma[idx] /= vector_norm(Photon_plasma_momentum + e_r, 3);

    }
    
    double inv_Boost_matrix[4][4]{};
    get_Lorentz_boost_matrix(inv_Boost_matrix, Disk_ZAMO_velocity, true);

    double Polarization_vector_ZAMO[4]{};
    mat_vec_multiply_4D(inv_Boost_matrix, Polarization_vector_plasma, Polarization_vector_ZAMO);

    ZAMO_to_Contravariant_coord(&s_Metric, Polarization_vector_ZAMO, this->Source_polarization_vector);

    return this->Source_polarization_vector;

}
