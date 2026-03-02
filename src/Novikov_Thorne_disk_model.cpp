#include "Novikov_Thorne_Model.h"

static double NT_flux_integrand_wrapper(double r, void* params) {

    /* This is a substitute for the actual state vector - only the radial coordinate is relevant - the functions that are called do not index anything else. */
    double Local_State_Vector[3] = { 0, r, M_PI / 2 };

    double Disk_Energy = static_cast<Novikov_Thorne_Model_class*>(params)->get_Disk_Energy(Local_State_Vector);
    double Disk_Ang_Velocity = static_cast<Novikov_Thorne_Model_class*>(params)->get_Disk_Angular_Velocity(Local_State_Vector);
    double Disk_Ang_Momentum = static_cast<Novikov_Thorne_Model_class*>(params)->get_Disk_Angular_Momentum(Local_State_Vector);
    double Disk_dr_Ang_Momentum = static_cast<Novikov_Thorne_Model_class*>(params)->get_dr_Disk_Angular_Momentum(Local_State_Vector);

    return (Disk_Energy - Disk_Ang_Velocity * Disk_Ang_Momentum) * Disk_dr_Ang_Momentum;

}

Novikov_Thorne_Model_class::Novikov_Thorne_Model_class(Simulation_Context_type* p_Sim_Context) {

    this->r_in = p_Sim_Context->p_Init_Conditions->Disk_params.Novikov_Thorne_params.r_in;
    this->r_out = p_Sim_Context->p_Init_Conditions->Disk_params.Novikov_Thorne_params.r_out;

    this->p_Spacetime = p_Sim_Context->p_Spacetime;
    this->e_Spacetime = p_Sim_Context->p_Init_Conditions->Metric_parameters.e_Spacetime;

    this->e_Mag_field_geometry = p_Sim_Context->p_Init_Conditions->Disk_params.e_Mag_field_geometry;

    this->Flux_integral_fucntion_struct.function = &NT_flux_integrand_wrapper;
    this->Flux_integral_workspace = gsl_integration_cquad_workspace_alloc(20000);

    memset(this->Disk_veclovity_vector, 0, 4 * sizeof(double));
    memset(this->Source_polarization_vector, 0, 4 * sizeof(double));
    memcpy(this->Mag_field_geometry, p_Sim_Context->p_Init_Conditions->Disk_params.Mag_field_geometry, 3 * sizeof(double));

    const int Flux_integral_interpolat_size = 1500;

    this->Flux_integral_spline_instance = gsl_spline_alloc(gsl_interp_cspline, Flux_integral_interpolat_size);
    this->Flux_integral_accelerator = gsl_interp_accel_alloc();

    this->Flux_integral = new double[Flux_integral_interpolat_size];
    this->Flux_r_coords = new double[Flux_integral_interpolat_size];

    for (int idx = 0; idx < Flux_integral_interpolat_size; idx++) {

        this->Flux_r_coords[idx] = this->r_in * (1 - double(idx) / (Flux_integral_interpolat_size - 1)) + double(idx) / (Flux_integral_interpolat_size - 1) * this->r_out;

        double Local_State_Vector[e_Full_state_size]{};
        Local_State_Vector[e_r] = this->Flux_r_coords[idx];
        Local_State_Vector[e_theta] = M_PI_2;

        this->Flux_integral[idx] = this->get_Flux(Local_State_Vector);

    }

    gsl_spline_init(this->Flux_integral_spline_instance, this->Flux_r_coords, Flux_integral, Flux_integral_interpolat_size);

}

Novikov_Thorne_Model_class::~Novikov_Thorne_Model_class(){

    gsl_integration_cquad_workspace_free(this->Flux_integral_workspace);
    gsl_spline_free(this->Flux_integral_spline_instance);
    gsl_interp_accel_free(this->Flux_integral_accelerator);

    free(this->Flux_integral);
    free(this->Flux_r_coords);

}

double Novikov_Thorne_Model_class::get_Disk_Angular_Velocity(const double* const Local_State_Vector) const {

    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_local_metric(Local_State_Vector);

    double Angular_velocity = (-s_dr_Metric.Metric[e_t][e_phi] + sqrt(s_dr_Metric.Metric[e_t][e_phi] * s_dr_Metric.Metric[e_t][e_phi] - s_dr_Metric.Metric[e_t][e_t] * s_dr_Metric.Metric[e_phi][e_phi])) / s_dr_Metric.Metric[e_phi][e_phi];

    if (isnan(Angular_velocity) or isinf(Angular_velocity) or isnan(Angular_velocity) or isinf(Angular_velocity)) {

        throw std::runtime_error(std::format("Invalid Novikov-Thorne angular velocity at: r = {}: Omega = {}", Local_State_Vector[e_r], Angular_velocity));

    }

    return Angular_velocity;

}

double Novikov_Thorne_Model_class::get_dr_Disk_Angular_Velocity(const double* const Local_State_Vector) const {

    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_local_metric(Local_State_Vector);
    Metric_type s_d2r_Metric = this->p_Spacetime->get_d2r_local_metric(Local_State_Vector);

    double root = sqrt(s_dr_Metric.Metric[e_t][e_phi] * s_dr_Metric.Metric[e_t][e_phi] - s_dr_Metric.Metric[e_t][e_t] * s_dr_Metric.Metric[e_phi][e_phi]);

    double Kepler = this->get_Disk_Angular_Velocity(Local_State_Vector);

    double dr_Kepler = -Kepler / s_dr_Metric.Metric[e_phi][e_phi] * s_d2r_Metric.Metric[e_phi][e_phi] + (-s_d2r_Metric.Metric[e_t][e_phi]
                     + 1.0 / root / 2 * (2 * s_dr_Metric.Metric[e_t][e_phi] * s_d2r_Metric.Metric[e_t][e_phi] - s_dr_Metric.Metric[e_t][e_t] * s_d2r_Metric.Metric[e_phi][e_phi]
                     - s_d2r_Metric.Metric[e_t][e_t] * s_dr_Metric.Metric[e_phi][e_phi])) / s_dr_Metric.Metric[e_phi][e_phi];

    if (isnan(dr_Kepler) or isinf(dr_Kepler) or isnan(dr_Kepler) or isinf(dr_Kepler)) {

        throw std::runtime_error(std::format("Invalid Novikov-Thorne angular velocity derivative at: r = {}: dr_Omega = {}", Local_State_Vector[e_r], dr_Kepler));

    }

    return dr_Kepler;
}

double* Novikov_Thorne_Model_class::get_Disk_Velocity_Vector(const double* const Local_State_Vector) {

    Metric_type s_Metric_source = this->p_Spacetime->get_local_metric(Local_State_Vector);

    double Angular_velocity = this->get_Disk_Angular_Velocity(Local_State_Vector);
    double u_t = 1 / sqrt(-s_Metric_source.Metric[e_t][e_t] - 2 * s_Metric_source.Metric[e_t][e_phi] * Angular_velocity - s_Metric_source.Metric[e_phi][e_phi] * Angular_velocity * Angular_velocity);

    if (isnan(u_t) or isinf(u_t) or isnan(Angular_velocity) or isinf(Angular_velocity)) {

        throw std::runtime_error(std::format("Invalid Novikov-Thorne disk 4-velocity at r = {}: u_t = {}, Angular velocity = {}", Local_State_Vector[e_r], u_t, Angular_velocity));

    }

    this->Disk_veclovity_vector[e_t] = u_t;
    this->Disk_veclovity_vector[e_r] = 0.0;
    this->Disk_veclovity_vector[e_theta] = 0.0;
    this->Disk_veclovity_vector[e_phi] = u_t * Angular_velocity;

    return this->Disk_veclovity_vector;

}

double Novikov_Thorne_Model_class::get_Disk_Energy(const double* const Local_State_Vector) const {

    Metric_type s_Metric_source = this->p_Spacetime->get_local_metric(Local_State_Vector);

    double Kepler = this->get_Disk_Angular_Velocity(Local_State_Vector);

    double root = sqrt(-s_Metric_source.Metric[e_t][e_t] - 2 * s_Metric_source.Metric[e_t][e_phi] * Kepler - s_Metric_source.Metric[e_phi][e_phi] * Kepler * Kepler);
    double Disk_Energy = -(s_Metric_source.Metric[e_t][e_t] + s_Metric_source.Metric[e_t][e_phi] * Kepler) / root;

    if (isnan(Disk_Energy) or isinf(Disk_Energy) or isnan(Disk_Energy) or isinf(Disk_Energy)) {

        throw std::runtime_error(std::format("Invalid Novikov-Thorne disk energy at r = {}: E = {}", Local_State_Vector[e_r], Disk_Energy));

    }

    return Disk_Energy;

}

double Novikov_Thorne_Model_class::get_Disk_Angular_Momentum(const double* const Local_State_Vector) const {

    Metric_type s_Metric_source = this->p_Spacetime->get_local_metric(Local_State_Vector);

    double Kepler = this->get_Disk_Angular_Velocity(Local_State_Vector);

    double root = sqrt(-s_Metric_source.Metric[e_t][e_t] - 2 * s_Metric_source.Metric[e_t][e_phi] * Kepler - s_Metric_source.Metric[e_phi][e_phi] * Kepler * Kepler);
    double Disk_angular_momentum = (s_Metric_source.Metric[e_phi][e_phi] * Kepler + s_Metric_source.Metric[e_t][e_phi]) / root;

    if (isnan(Disk_angular_momentum) or isinf(Disk_angular_momentum) or isnan(Disk_angular_momentum) or isinf(Disk_angular_momentum)) {

        throw std::runtime_error(std::format("Invalid Novikov-Thorne disk angular momentum at r = {}: L_z = {}", Local_State_Vector[e_r], Disk_angular_momentum));

    }

    return  Disk_angular_momentum;

}

double Novikov_Thorne_Model_class::get_dr_Disk_Angular_Momentum(const double* const Local_State_Vector) const {

    Metric_type s_Metric = this->p_Spacetime->get_local_metric(Local_State_Vector);
    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_local_metric(Local_State_Vector);

    double Disk_Ang_Velocity = this->get_Disk_Angular_Velocity(Local_State_Vector);
    double Disk_dr_Ang_Velocity = this->get_dr_Disk_Angular_Velocity(Local_State_Vector);
    double Disk_Ang_Momentum = this->get_Disk_Angular_Momentum(Local_State_Vector);

    double root = sqrt(-s_Metric.Metric[e_t][e_t] - 2 * s_Metric.Metric[e_t][e_phi] * Disk_Ang_Velocity - s_Metric.Metric[e_phi][e_phi] * Disk_Ang_Velocity * Disk_Ang_Velocity);
    double dr_root = (-s_dr_Metric.Metric[e_t][e_t] - 2 * (s_dr_Metric.Metric[e_t][e_phi] * Disk_Ang_Velocity + s_Metric.Metric[e_t][e_phi] * Disk_dr_Ang_Velocity)
                     - s_dr_Metric.Metric[e_phi][e_phi] * Disk_Ang_Velocity * Disk_Ang_Velocity - 2 * s_Metric.Metric[e_phi][e_phi] * Disk_Ang_Velocity * Disk_dr_Ang_Velocity);

    double Disk_dr_angular_momentum = (s_dr_Metric.Metric[e_phi][e_phi] * Disk_Ang_Velocity + s_Metric.Metric[e_phi][e_phi] * Disk_dr_Ang_Velocity + s_dr_Metric.Metric[e_t][e_phi]) / root - Disk_Ang_Momentum / root / root / 2 * dr_root;


    if (isnan(Disk_Ang_Momentum) or isinf(Disk_Ang_Momentum) or isnan(Disk_Ang_Momentum) or isinf(Disk_Ang_Momentum)) {

        throw std::runtime_error(std::format("Invalid Novikov-Thorne disk angular momentum derivative at r = {}: dr_L_z = {}", Local_State_Vector[e_r], Disk_Ang_Momentum));

    }

    return Disk_dr_angular_momentum;

}

double Novikov_Thorne_Model_class::get_Interpolated_Flux(const double* const Local_State_Vector) const {

    return gsl_spline_eval(this->Flux_integral_spline_instance, Local_State_Vector[e_r], this->Flux_integral_accelerator);

}

double Novikov_Thorne_Model_class::get_Flux(double* Local_State_Vector) {

    double Equatorial_metric_det = get_eq_induced_metric_det(this->p_Spacetime->get_local_metric(Local_State_Vector).Metric);

    double Disk_Energy = this->get_Disk_Energy(Local_State_Vector);
    double Disk_Ang_Momentum = this->get_Disk_Angular_Momentum(Local_State_Vector);

    double Disk_Ang_Velocity = this->get_Disk_Angular_Velocity(Local_State_Vector);
    double Disk_dr_Ang_Velocity = this->get_dr_Disk_Angular_Velocity(Local_State_Vector);

    double Flux_coeff = -Disk_dr_Ang_Velocity / ((Disk_Energy - Disk_Ang_Velocity * Disk_Ang_Momentum) * (Disk_Energy - Disk_Ang_Velocity * Disk_Ang_Momentum)) / (4 * M_PI * sqrt(-Equatorial_metric_det));

    double Flux_integral{};
    double Error_estimate{};

    this->Flux_integral_fucntion_struct.params = static_cast<void*>(this);
       
    gsl_integration_cquad(&this->Flux_integral_fucntion_struct,
                           this->r_in, 
                           Local_State_Vector[e_r], 
                           0,
                           1e-6,
                           this->Flux_integral_workspace,
                          &Flux_integral,
                          &Error_estimate,
                          nullptr);

    return Flux_coeff * Flux_integral;

}

double* Novikov_Thorne_Model_class::Construct_coord_polarization_vector(const double* const State_at_event_local) {

    Metric_type s_Metric = this->p_Spacetime->get_local_metric(State_at_event_local);

    /* NOTE: This is the covariant momentum */
    const double* const &Photon_coordinate_momentum_covariant = State_at_event_local + e_p_t;

    double Photon_coordinate_momentum_contravariant[4]{};
    Manipulate_index(&s_Metric, Photon_coordinate_momentum_covariant, Photon_coordinate_momentum_contravariant, Raise_index);

    double Photon_ZAMO_momentum[4]{};
    Contravariant_coord_to_ZAMO(&s_Metric, Photon_coordinate_momentum_contravariant, Photon_ZAMO_momentum);

    double Disk_ZAMO_velocity[4]{};
    Contravariant_coord_to_ZAMO(&s_Metric, this->get_Disk_Velocity_Vector(State_at_event_local), Disk_ZAMO_velocity);

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
        Polarization_vector_plasma[idx] *= sqrt(std::abs(Photon_plasma_momentum[e_t] / Photon_plasma_momentum[e_theta]));
    }
    
    double inv_Boost_matrix[4][4]{};
    get_Lorentz_boost_matrix(inv_Boost_matrix, Disk_ZAMO_velocity, true);

    double Polarization_vector_ZAMO[4]{};
    mat_vec_multiply_4D(inv_Boost_matrix, Polarization_vector_plasma, Polarization_vector_ZAMO);

    ZAMO_to_Contravariant_coord(&s_Metric, Polarization_vector_ZAMO, this->Source_polarization_vector);

    return this->Source_polarization_vector;

}
