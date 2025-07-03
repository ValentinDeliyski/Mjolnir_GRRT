#include "Page_Thorne_Model.h"

Page_Thorne_Model_class::Page_Thorne_Model_class(Simulation_Context_type* p_Sim_Context) {

    this->r_in  = p_Sim_Context->p_Init_Conditions->Disk_params.Page_Thorne_params.r_in;
    this->r_out = p_Sim_Context->p_Init_Conditions->Disk_params.Page_Thorne_params.r_out;

    this->flux_integral_accuracy = p_Sim_Context->p_Init_Conditions->Integrator_params.Simpson_accuracy;

    this->p_Spacetime = p_Sim_Context->p_Spacetime;
    this->e_Spacetime = p_Sim_Context->p_Init_Conditions->Metric_parameters.e_Spacetime;

};

double Page_Thorne_Model_class::Keplerian_angular_velocity(const double* const State_Vector) {

    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_metric(State_Vector);

    return (-s_dr_Metric.Metric[e_t][e_phi] - sqrt(s_dr_Metric.Metric[e_t][e_phi] * s_dr_Metric.Metric[e_t][e_phi] - s_dr_Metric.Metric[e_t][e_t] * s_dr_Metric.Metric[e_phi][e_phi])) / s_dr_Metric.Metric[e_phi][e_phi];

}

double Page_Thorne_Model_class::dr_Keplerian_angular_velocity(const double* const State_Vector) {

    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_metric(State_Vector);
    Metric_type s_d2r_Metric = this->p_Spacetime->get_d2r_metric(State_Vector);

    double root = sqrt(s_dr_Metric.Metric[e_t][e_phi] * s_dr_Metric.Metric[e_t][e_phi] - s_dr_Metric.Metric[e_t][e_t] * s_dr_Metric.Metric[e_phi][e_phi]);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

    return  - Kepler / s_dr_Metric.Metric[e_phi][e_phi] * s_d2r_Metric.Metric[e_phi][e_phi] + (-s_d2r_Metric.Metric[e_t][e_phi]
            + 1.0 / root / 2 * (2 * s_dr_Metric.Metric[e_t][e_phi] * s_d2r_Metric.Metric[e_t][e_phi] - s_dr_Metric.Metric[e_t][e_t] * s_d2r_Metric.Metric[e_phi][e_phi]
            - s_d2r_Metric.Metric[e_t][e_t] * s_dr_Metric.Metric[e_phi][e_phi])) / s_dr_Metric.Metric[e_phi][e_phi];

}

double Page_Thorne_Model_class::Redshift(const double* const State_Vector, Observer_class* const p_Observer) {

    const double& r_source  = State_Vector[e_r];
    const double& theta_source = State_Vector[e_theta];

    const double* U_obs = p_Observer->get_obs_velocity();

    /* Get the source 4-velocity */
    Metric_type s_Metric_source = this->p_Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);
    double Gamma  = 1 / sqrt(-s_Metric_source.Metric[0][0] - 2 * s_Metric_source.Metric[0][3] * Kepler - s_Metric_source.Metric[3][3] * Kepler * Kepler);

    if (isnan(Gamma) || isinf(Gamma) || isnan(Kepler) || isinf(Kepler)) {

        std::cout << "Invalid NT disk 4-velocity: "
            << "Gamma = "
            << Gamma
            << "\n"
            << "Kepler = "
            << "\n"
            << Kepler
            << "\n";

        exit(ERROR);

    }

    double U_source[4] = { Gamma, 0, 0, Gamma * Kepler };

    return  (-U_obs[0] + U_obs[3] * State_Vector[e_p_phi]) / (-U_source[0] + U_source[3] * State_Vector[e_p_phi]);


}

double Page_Thorne_Model_class::disk_Energy(const double* const State_Vector) {

    Metric_type s_Metric_source = this->p_Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

    double root = sqrt(-s_Metric_source.Metric[e_t][e_t] - 2 * s_Metric_source.Metric[e_t][e_phi] * Kepler - s_Metric_source.Metric[e_phi][e_phi] * Kepler * Kepler);

    return  -(s_Metric_source.Metric[e_t][e_t] + s_Metric_source.Metric[e_t][e_phi] * Kepler) / root;

}

double Page_Thorne_Model_class::disk_Angular_Momentum(const double* const State_Vector) {

    Metric_type s_Metric_source = this->p_Spacetime->get_metric(State_Vector);

    double Kepler = this->Keplerian_angular_velocity(State_Vector);

    double root = sqrt(-s_Metric_source.Metric[e_t][e_t] - 2 * s_Metric_source.Metric[e_t][e_phi] * Kepler - s_Metric_source.Metric[e_phi][e_phi] * Kepler * Kepler);

    return  (s_Metric_source.Metric[e_phi][e_phi] * Kepler + s_Metric_source.Metric[e_t][e_phi]) / root;

}

double Page_Thorne_Model_class::Flux_integrand(const double* const State_Vector) {

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

    return (E - Kepler * L) * dr_L;

}

double Page_Thorne_Model_class::solve_Flux_integral(double r_in, const double* const State_Vector, double tolerance) {

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

    double error;

    if (S_2 >= S_1) {

        error = S_2 - S_1;

    }
    else {

        error = S_1 - S_2;

    }

    double integral;

    if (error < 15 * tolerance) {

        integral = S_2 + (S_2 - S_1) / 15;


    }
    else {

        double L_value = this->solve_Flux_integral(lower_bound, mid_point_state_vector, tolerance / 2);
        double R_value = this->solve_Flux_integral(mid_point, State_Vector, tolerance / 2);

        integral = L_value + R_value;

    }

    return integral;
}

double Page_Thorne_Model_class::get_flux(const double* const State_Vector) {

    Metric_type s_Metric = this->p_Spacetime->get_metric(State_Vector);

    double metric_det = get_eq_induced_metric_det(s_Metric.Metric);
    double E_disk = disk_Energy(State_Vector);
    double L_disk = disk_Angular_Momentum(State_Vector);

    double Kepler = Keplerian_angular_velocity(State_Vector);
    double dr_Kepler = dr_Keplerian_angular_velocity(State_Vector);

    double Flux_coeff = -dr_Kepler / ((E_disk - Kepler * L_disk) * (E_disk - Kepler * L_disk)) / (4 * M_PI * sqrt(-metric_det));

    double Flux_integral = solve_Flux_integral(this->r_in, State_Vector, this->flux_integral_accuracy);

    return Flux_coeff * Flux_integral;

}
