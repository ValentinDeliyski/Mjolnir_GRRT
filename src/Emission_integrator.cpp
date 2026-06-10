#include "Emission_integrator.h"
#include "Emission_Models.h"

Emission_Integrator_class::Emission_Integrator_class(const Simulation_Context_type* p_Sim_Context, Results_type* const p_Ray_results) {

    this->e_Active_rad_transfer_integrator = p_Sim_Context->p_Init_Conditions->Integrator_params.e_Radiative_transfer_integrator;
    this->e_Active_parallel_transport_integrator = p_Sim_Context->p_Init_Conditions->Integrator_params.e_Parallel_transport_integrator;

    this->p_Sim_Context = p_Sim_Context;
    this->Ray_log_length = p_Ray_results->Ray_log_struct.Log_length;

    for (int idx = I; idx < e_Stokes_param_num; idx ++) {

        this->Emission_log[idx] = p_Ray_results->Ray_log_struct.Ray_emission_log[idx];

    }

    for (int idx = e_x; idx <= e_y; idx++) {

        this->Polarization_log[idx] = p_Ray_results->Ray_log_struct.Ray_polarization_log[idx];
        this->PW_Constant_log[idx] = p_Ray_results->Polarization_debug_log.PW_constant[idx];
    }

    /* --------------------------------------------------- Init the counters -------------------------------------------------- */

    this->Current_emission_log_idx = 0;
    this->Current_polarization_log_idx = 0;

    /* ------------------------------------------ Construct the geodesic ray splines ------------------------------------------ */

    std::unique_ptr<double[]> Affine_param_log = std::make_unique<double[]>(p_Ray_results->Ray_log_struct.Log_length);

    std::unique_ptr<double[]> Ray_Log[e_Dynamic_state_size];
    std::unique_ptr<double[]> Radial_Ray_log_global[2];
    
    for (size_t idx = 0; idx < e_Dynamic_state_size; idx++) {

        Ray_Log[idx] = std::make_unique<double[]>(p_Ray_results->Ray_log_struct.Log_length);

    }

    Radial_Ray_log_global[0] = std::make_unique<double[]>(p_Ray_results->Ray_log_struct.Log_length);
    Radial_Ray_log_global[1] = std::make_unique<double[]>(p_Ray_results->Ray_log_struct.Log_length);


    for (size_t log_idx = 0; log_idx < p_Ray_results->Ray_log_struct.Log_length; log_idx++) {

        for (int component_idx = 0; component_idx < e_Dynamic_state_size; component_idx++) {

            Ray_Log[component_idx][p_Ray_results->Ray_log_struct.Log_length - 1 - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_local[component_idx + log_idx * e_Full_state_size];

        }

        Radial_Ray_log_global[0][p_Ray_results->Ray_log_struct.Log_length - 1 - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_global[e_r + log_idx * e_Full_state_size];
        Radial_Ray_log_global[1][p_Ray_results->Ray_log_struct.Log_length - 1 - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_global[e_p_r + log_idx * e_Full_state_size];

        Affine_param_log[p_Ray_results->Ray_log_struct.Log_length - 1 - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_local[e_ray_affine_param + log_idx * e_Full_state_size];
    }

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        this->p_Ray_spline_instance[idx] = gsl_spline_alloc(gsl_interp_cspline, p_Ray_results->Ray_log_struct.Log_length);
        gsl_spline_init(this->p_Ray_spline_instance[idx], Affine_param_log.get(), Ray_Log[idx].get(), p_Ray_results->Ray_log_struct.Log_length);

    }

    /* ---------------------------------------- Radial spline in global coordinates ---------------------------------------- */

    this->p_Ray_spline_instance[e_Dynamic_state_size + 0] = gsl_spline_alloc(gsl_interp_cspline, p_Ray_results->Ray_log_struct.Log_length);
    gsl_spline_init(this->p_Ray_spline_instance[e_Dynamic_state_size + 0], Affine_param_log.get(), Radial_Ray_log_global[0].get(), p_Ray_results->Ray_log_struct.Log_length);

    this->p_Ray_spline_instance[e_Dynamic_state_size + 1] = gsl_spline_alloc(gsl_interp_cspline, p_Ray_results->Ray_log_struct.Log_length);
    gsl_spline_init(this->p_Ray_spline_instance[e_Dynamic_state_size + 1], Affine_param_log.get(), Radial_Ray_log_global[1].get(), p_Ray_results->Ray_log_struct.Log_length);

    this->p_Spline_accelerator = gsl_interp_accel_alloc();

}

Emission_Integrator_class::~Emission_Integrator_class() {

    for (int idx = 0; idx < e_Dynamic_state_size + 2; idx++) {

        gsl_spline_free(this->p_Ray_spline_instance[idx]);

    }

    gsl_interp_accel_free(this->p_Spline_accelerator);

}

void Emission_Integrator_class::Run_Runge_Kutta_Stokes_Vector(const double Start_Affine_Param, const double End_Affine_Param) {

    if (RK78_Fehlberg != e_Active_rad_transfer_integrator and RK78_DP != e_Active_rad_transfer_integrator and RK54 != e_Active_rad_transfer_integrator) {

        throw std::runtime_error("Wrong active integrator in Emission_Integrator_class::Run_Explicit_Runge_Kutta()!");

    }

    memset(this->Rad_Transfer_RHS_log, 0, sizeof(double) * RK78_size * e_Stokes_param_num);

    int RK_size = RK78_size;
    auto Stage_coeff = this->RK78_DP_Coeff_deriv;
    auto Main_solution_coeff = this->RK78_DP_Coeff_sol_main;
    auto Embedded_solution_coeff = this->RK78_DP_Coeff_sol_embeded;
    auto Affine_param_coeff = this->RK78_DP_Coeff_affine_param;

    if (RK78_Fehlberg == this->e_Active_rad_transfer_integrator) {

        Stage_coeff = this->RK78_Fhelberg_Coeff_deriv;
        Main_solution_coeff = this->RK78_Fhelberg_Coeff_sol_main;
        Embedded_solution_coeff = this->RK78_Fhelberg_Coeff_sol_embeded;
        Affine_param_coeff = this->RK78_Fhelberg_Coeff_affine_param;

    }
    else if (RK54 == this->e_Active_rad_transfer_integrator) {

        Stage_coeff = this->RK54_Coeff_deriv;
        Main_solution_coeff = this->RK54_Coeff_sol_main;
        Embedded_solution_coeff = this->RK54_Coeff_test_embeded;
        Affine_param_coeff = this->RK54_Coeff_affine_param;

        RK_size = RK54_size;

    }

    double Temp_affine_param{};

    const double Geometric_Step = End_Affine_Param - Start_Affine_Param;
    const double CGS_Step = Geometric_Step * MASS_TO_CM * this->p_Sim_Context->p_Init_Conditions->central_object_mass;
    this->Current_Affine_Param = Start_Affine_Param;

    for (int RK_stage = 0; RK_stage < RK_size; RK_stage++) {

        Temp_affine_param = this->Current_Affine_Param + Affine_param_coeff[RK_stage] * Geometric_Step;

        /* ------------------------------------------------------------ Get the radiative transfer RHS ------------------------------------------------------------ */

        Transfer_functions_type Total_Transfer_Functions{};

        for (int emission_medium = Disk; emission_medium < e_Emission_medium_number; emission_medium++) {

            Transfer_functions_type Temp_Transfer_functions{};

            this->p_Sim_Context->p_Emission_Model->get_radiative_transfer_functions(this->get_ray_Local_State_Vector(Temp_affine_param),
                                                                                    static_cast<Emission_medium_enums>(emission_medium),
                                                                                    &Temp_Transfer_functions);

            add_vectors(Temp_Transfer_functions.Absorbtion_functions, Total_Transfer_Functions.Absorbtion_functions, e_Stokes_param_num, Total_Transfer_Functions.Absorbtion_functions);
            add_vectors(Temp_Transfer_functions.Emission_functions, Total_Transfer_Functions.Emission_functions, e_Stokes_param_num, Total_Transfer_Functions.Emission_functions);
            add_vectors(Temp_Transfer_functions.Faradey_functions, Total_Transfer_Functions.Faradey_functions, e_Stokes_param_num, Total_Transfer_Functions.Faradey_functions);

        }

        this->get_Radiative_transfer_RHS(Total_Transfer_Functions.Emission_functions,
                                         Total_Transfer_Functions.Absorbtion_functions,
                                         Total_Transfer_Functions.Faradey_functions,
                                         this->Temp_Stokes_Vector,
                                         this->Rad_Transfer_RHS_log + RK_stage * e_Stokes_param_num);

        /* -------------------------------------------------------- Propagate the intermediate Stokes Vector -------------------------------------------------------- */

        memcpy(this->Temp_Stokes_Vector, this->Current_Stokes_Vector, e_Stokes_param_num * sizeof(double));

        for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {

            for (int derivative_idx = 0; derivative_idx < RK_stage; derivative_idx++) {

                this->Temp_Stokes_Vector[state_idx] += CGS_Step * Stage_coeff[RK_stage][derivative_idx] * this->Rad_Transfer_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];

            }
        }

        /* ---------------------------------------------------------------------------------------------------------------------------------------------------------- */

    }

    this->Current_Affine_Param += Geometric_Step;

    for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {

        for (int derivative_idx = 0; derivative_idx < RK_size; derivative_idx++) {

            this->Current_Stokes_Vector[state_idx] += CGS_Step * Main_solution_coeff[derivative_idx] * this->Rad_Transfer_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];

        }

    }

    this->Update_emission_log();

}

void Emission_Integrator_class::get_Parallel_Transport_RHS(const double* const Global_State_Vector,
                                                           std::complex<double>* Vector_to_transport,
                                                           Tensor_type_enums e_Vec_type,
                                                           std::complex<double>* const RHS) {

    Metric_type s_Metric = this->p_Sim_Context->p_Spacetime->get_global_metric(Global_State_Vector);
    Metric_type s_dr_Metric = this->p_Sim_Context->p_Spacetime->get_dr_global_metric(Global_State_Vector);
    Metric_type s_dtheta_Metric = this->p_Sim_Context->p_Spacetime->get_dtheta_global_metric(Global_State_Vector);

    double inv_Metric[4][4]{};
    invert_metric(inv_Metric, s_Metric.Metric);

    double Connection_Coefficients[4][4][4]{};
    get_connection_coefficients(s_Metric, s_dr_Metric, s_dtheta_Metric, Connection_Coefficients);

    double Photon_momentum_contravariant[4]{};
    Manipulate_index(&s_Metric, Global_State_Vector + e_p_t, Photon_momentum_contravariant, Raise_index);

    /* ========================== Compute the derivative of the polarization vector from the parallel transport ========================== */

    for (int derivative_index = 0; derivative_index < e_Stokes_param_num; derivative_index++) {

        for (int polarization_index = 0; polarization_index < e_Stokes_param_num; polarization_index++) {

            for (int wave_vector_index = 0; wave_vector_index < e_Stokes_param_num; wave_vector_index++) {

                if (e_Vec_type == Contravariant) {

                    RHS[derivative_index] += -Connection_Coefficients[derivative_index][wave_vector_index][polarization_index] * Photon_momentum_contravariant[wave_vector_index] * Vector_to_transport[polarization_index];

                }
                else {

                    RHS[derivative_index] += Connection_Coefficients[polarization_index][derivative_index][wave_vector_index] * Photon_momentum_contravariant[wave_vector_index] * Vector_to_transport[polarization_index];

                }

            }

        }

    }

};

void Emission_Integrator_class::get_Radiative_transfer_RHS(const double* const Emission_Functions,
                                                           const double* const Absorbtion_Functions,
                                                           const double* const Faradey_Functions,
                                                           const double* const Stokes_Vector,
                                                           double* const RHS) {

    double M_matrix[4][4]{};

    M_matrix[0][0] = M_matrix[1][1] = M_matrix[2][2] = M_matrix[3][3] = Absorbtion_Functions[I];
    M_matrix[0][1] = M_matrix[1][0] = Absorbtion_Functions[Q];
    M_matrix[0][2] = M_matrix[2][0] = Absorbtion_Functions[U];
    M_matrix[0][3] = M_matrix[3][0] = Absorbtion_Functions[V];
    M_matrix[1][2] = Faradey_Functions[V];
    M_matrix[2][1] = -Faradey_Functions[V];
    M_matrix[1][3] = -Faradey_Functions[U];
    M_matrix[3][1] = Faradey_Functions[U];
    M_matrix[2][3] = Faradey_Functions[Q];
    M_matrix[3][2] = -Faradey_Functions[Q];

    double M_dot_Stokes[4]{};
    mat_vec_multiply_4D(M_matrix, Stokes_Vector, M_dot_Stokes);
    
    for (int idx = 0; idx < e_Stokes_param_num; idx++) {
    
        RHS[idx] = Emission_Functions[idx] - M_dot_Stokes[idx];
    
    }

}

void Emission_Integrator_class::Update_emission_log() {

    for (int idx = 0; idx < e_Stokes_param_num; idx++) {

        this->Emission_log[idx][this->Current_emission_log_idx] = this->Current_Stokes_Vector[idx];

    }

    this->Current_emission_log_idx += 1;

}

const double* const Emission_Integrator_class::get_ray_Local_State_Vector(const double Affine_param) {

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        this->Current_State_Vector[idx] = gsl_spline_eval(this->p_Ray_spline_instance[idx], Affine_param, this->p_Spline_accelerator);

    }

    return this->Current_State_Vector;

}

const double* const Emission_Integrator_class::get_ray_Global_State_Vector(const double Affine_param) {

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        if (idx == e_r) {

            this->Current_State_Vector[e_r] = gsl_spline_eval(this->p_Ray_spline_instance[e_Dynamic_state_size + 0], Affine_param, this->p_Spline_accelerator);
            continue;

        }

        if (idx == e_p_r) {

            this->Current_State_Vector[e_p_r] = gsl_spline_eval(this->p_Ray_spline_instance[e_Dynamic_state_size + 1], Affine_param, this->p_Spline_accelerator);
            continue;

        }

        this->Current_State_Vector[idx] = gsl_spline_eval(this->p_Ray_spline_instance[idx], Affine_param, this->p_Spline_accelerator);

    }

    return this->Current_State_Vector;

}

const double* const Emission_Integrator_class::get_current_Stokes_Vector() const{

    return this->Current_Stokes_Vector;

}

void Emission_Integrator_class::Propagate_Polarization_Vector(const double Start_Affine_Param, const double End_Affine_Param, const Tensor_type_enums e_Vec_type) {

    if (RK78_Fehlberg != this->e_Active_parallel_transport_integrator and RK78_DP != this->e_Active_parallel_transport_integrator and RK54 != this->e_Active_parallel_transport_integrator) {

        throw std::runtime_error("Wrong active integrator in Emission_Integrator_class::Propagate_Polarization_Vector()!");

    }

    memset(this->Parallel_Transport_RHS_log, 0, sizeof(std::complex<double>) * RK78_size * e_Stokes_param_num);

    int RK_size = RK78_size;
    auto Stage_coeff = this->RK78_DP_Coeff_deriv;
    auto Main_solution_coeff = this->RK78_DP_Coeff_sol_main;
    auto Embedded_solution_coeff = this->RK78_DP_Coeff_sol_embeded;
    auto Affine_param_coeff = this->RK78_DP_Coeff_affine_param;

    if (RK78_Fehlberg == this->e_Active_parallel_transport_integrator) {

        Stage_coeff = this->RK78_Fhelberg_Coeff_deriv;
        Main_solution_coeff = this->RK78_Fhelberg_Coeff_sol_main;
        Embedded_solution_coeff = this->RK78_Fhelberg_Coeff_sol_embeded;
        Affine_param_coeff = this->RK78_Fhelberg_Coeff_affine_param;

    }
    else if (RK54 == this->e_Active_parallel_transport_integrator) {

        Stage_coeff = this->RK54_Coeff_deriv;
        Main_solution_coeff = this->RK54_Coeff_sol_main;
        Embedded_solution_coeff = this->RK54_Coeff_test_embeded;
        Affine_param_coeff = this->RK54_Coeff_affine_param;

        RK_size = RK54_size;

    }

    double Temp_affine_param{};

    const double Step = End_Affine_Param - Start_Affine_Param;
    this->Current_Affine_Param = Start_Affine_Param;

    for (int RK_stage = 0; RK_stage < RK_size; RK_stage++) {

        memcpy(this->Temp_Pol_Vector, this->Current_Pol_Vector, e_Stokes_param_num * sizeof(std::complex<double>));
    
        for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {
    
            for (int derivative_idx = 0; derivative_idx < RK_stage; derivative_idx++) {
    
                this->Temp_Pol_Vector[state_idx] += Step * Stage_coeff[RK_stage][derivative_idx] * this->Parallel_Transport_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];
    
            }
        }

        Temp_affine_param = this->Current_Affine_Param + Affine_param_coeff[RK_stage] * Step;

        /* --------------------------------------- Get the Parallel Transport RHS --------------------------------------- */

        this->get_Parallel_Transport_RHS(this->get_ray_Global_State_Vector(Temp_affine_param),
                                         this->Temp_Pol_Vector,
                                         e_Vec_type,
                                         this->Parallel_Transport_RHS_log + RK_stage * e_Stokes_param_num);

        /* -------------------------------------------------------------------------------------------------------------- */

    }
    
    for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {
    
        for (int derivative_idx = 0; derivative_idx < RK_size; derivative_idx++) {
    
            this->Current_Pol_Vector[state_idx] += Step * Main_solution_coeff[derivative_idx] * this->Parallel_Transport_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];
    
        }
    
    }
    
    this->Current_Affine_Param += Step;
    this->Update_polarization_log();

}

void Emission_Integrator_class::Map_Stokes_to_Polarization_Vector(const double Stokes_Tetrad[e_Stokes_param_num][e_Stokes_param_num], bool Map_Between_Intermediate) {

    double* Stokes_Vector_to_map = nullptr;
    std::complex<double>* Pol_Vector_to_map = nullptr;

    if (Map_Between_Intermediate) {

        Stokes_Vector_to_map = this->Temp_Stokes_Vector;
        Pol_Vector_to_map = this->Temp_Pol_Vector;

    }
    else {

        Stokes_Vector_to_map = this->Current_Stokes_Vector;
        Pol_Vector_to_map = this->Current_Pol_Vector;


    }

    memset(Pol_Vector_to_map, 0, 4 * sizeof(std::complex<double>));

    std::complex<double> Stokes_Basis_Pol_vec[e_Stokes_param_num]{};

    Stokes_Basis_Pol_vec[1] = 1.0 / std::numbers::sqrt2;
    Stokes_Basis_Pol_vec[2] = 1.0 / std::numbers::sqrt2;

    const double Polarized_Intensity = vector_norm(Stokes_Vector_to_map + 1, 3);

    if (Polarized_Intensity < std::numeric_limits<double>::min()) { return; }

    const double Normalized_Q = Stokes_Vector_to_map[Q] / Polarized_Intensity;
    const double Normalized_U = Stokes_Vector_to_map[U] / Polarized_Intensity;
    const double Normalized_V = Stokes_Vector_to_map[V] / Polarized_Intensity;

    if (std::abs(Normalized_Q) > 1 or std::abs(Normalized_V) > 1) {

        throw std::runtime_error("Normalized polarization components > 1 in Map_Stokes_to_Polarization_Vector()!");

    }

    Stokes_Basis_Pol_vec[1] = sqrt((1 + Normalized_Q) / 2);

    if (Stokes_Basis_Pol_vec[1].real() < std::numeric_limits<double>::min()) {

        Stokes_Basis_Pol_vec[2] = 1.0;

    }
    else {

        Stokes_Basis_Pol_vec[2] = (Normalized_U - complex_i * Normalized_V) / (2.0 * Stokes_Basis_Pol_vec[1]);

    }

    for (int coord_idx = 0; coord_idx < 4; coord_idx++) {

        for (int stokes_idx = 0; stokes_idx < 4; stokes_idx++) {

            Pol_Vector_to_map[coord_idx] += Stokes_Tetrad[stokes_idx][coord_idx] * Stokes_Basis_Pol_vec[stokes_idx];

        }

    }

}

void Emission_Integrator_class::Map_Polarization_Vector_to_Stokes(const double inv_Stokes_Tetrad[e_Stokes_param_num][e_Stokes_param_num], bool Map_Between_Intermediate) {

    std::complex<double> Stokes_Basis_Pol_vec[4]{};
    
    double* Stokes_Vector_to_map = nullptr;
    std::complex<double>* Pol_Vector_to_map = nullptr;

    if (Map_Between_Intermediate) {

        Stokes_Vector_to_map = this->Temp_Stokes_Vector;
        Pol_Vector_to_map = this->Temp_Pol_Vector;

    }
    else {

        Stokes_Vector_to_map = this->Current_Stokes_Vector;
        Pol_Vector_to_map = this->Current_Pol_Vector;

    }

    for (int stokes_idx = 0; stokes_idx < 4; stokes_idx++) {
    
        for (int coord_idx = 0; coord_idx < 4 ; coord_idx++) {
    
            Stokes_Basis_Pol_vec[stokes_idx] += inv_Stokes_Tetrad[stokes_idx][coord_idx] * Pol_Vector_to_map[coord_idx];

        }

    }

    Normalize_complex_vector(Stokes_Basis_Pol_vec, Minkowski_Metric, Contravariant);

    double Polarized_Intensity = vector_norm(Stokes_Vector_to_map + 1, 3);

    Stokes_Vector_to_map[Q] = (Polarized_Intensity * (Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[1]) -
                                                      Stokes_Basis_Pol_vec[2] * std::conj(Stokes_Basis_Pol_vec[2]))).real();
    
    Stokes_Vector_to_map[U] = (Polarized_Intensity * (Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[2]) +
                                                      Stokes_Basis_Pol_vec[2] * std::conj(Stokes_Basis_Pol_vec[1]))).real();
    
    Stokes_Vector_to_map[V] = (Polarized_Intensity * (Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[2]) -
                                                      Stokes_Basis_Pol_vec[2] * std::conj(Stokes_Basis_Pol_vec[1]))).imag();

}

void Emission_Integrator_class::Update_polarization_log() {

    std::complex<double> Polarization_vector_ZAMO[4];
    Metric_type s_Metric = this->p_Sim_Context->p_Spacetime->get_global_metric(this->get_ray_Global_State_Vector(this->Current_Affine_Param));

    Contravariant_coord_to_ZAMO(&s_Metric, this->Current_Pol_Vector, Polarization_vector_ZAMO);

    this->Polarization_log[e_x][this->Current_polarization_log_idx] = Polarization_vector_ZAMO[e_phi].real();
    this->Polarization_log[e_y][this->Current_polarization_log_idx] = Polarization_vector_ZAMO[e_theta].real();

    if (Spacetime_enums::Kerr == p_Sim_Context->p_Init_Conditions->Metric_parameters.e_Spacetime) {

        this->PW_Constant_log[e_x][this->Current_polarization_log_idx] = get_Penrose_Walker_constant(this->get_ray_Global_State_Vector(this->Current_Affine_Param), this->p_Sim_Context, this->Current_Pol_Vector).real();
        this->PW_Constant_log[e_y][this->Current_polarization_log_idx] = get_Penrose_Walker_constant(this->get_ray_Global_State_Vector(this->Current_Affine_Param), this->p_Sim_Context, this->Current_Pol_Vector).imag();

    }

    this->Current_polarization_log_idx += 1;

}

const std::complex<double>* const Emission_Integrator_class::get_current_Polarization_Vector() const {

    return this->Current_Pol_Vector;

};

void Emission_Integrator_class::normalize_polarization_vector(const double Affine_param) {

    Metric_type s_Metric = p_Sim_Context->p_Spacetime->get_global_metric(this->get_ray_Global_State_Vector(Affine_param));

    std::complex <double> norm{};

    for (int idx1 = 0; idx1 < 4; idx1++) {

        for (int idx2 = 0; idx2 < 4; idx2++) {

            norm += (s_Metric.Metric[idx1][idx2] * this->Current_Pol_Vector[idx1] * std::conj(this->Current_Pol_Vector[idx2]));

        }

    }

    norm = std::sqrt(norm);

    if (isinf(1.0 / norm.real()) or isnan(1.0 / norm.real())) { return; }

    for (int idx = 0; idx < 4; idx++) {

        this->Current_Pol_Vector[idx] = this->Current_Pol_Vector[idx] / norm;

    }

}

void Emission_Integrator_class::Propagate_Stokes_Vector(const double Start_Affine_Param, const double End_Affine_Param) {

    switch (this->e_Active_rad_transfer_integrator) {

    case Rad_Analytic:

        this->Run_Analytic_Stokes_Vector_Propagator(Start_Affine_Param, End_Affine_Param);

        break;

    default:

        this->Run_Runge_Kutta_Stokes_Vector(Start_Affine_Param, End_Affine_Param);
        break;
         
    }

}

void Emission_Integrator_class::Get_radiative_transfer_operators(const double* const Absorbtion_functions,
                                                                 const double* const Faradey_functions,
                                                                 double const CGS_Step,
                                                                 double Transfer_Operator[e_Stokes_param_num][e_Stokes_param_num],
                                                                 double Integrated_Transfer_Operator[e_Stokes_param_num][e_Stokes_param_num]) {
     
    /* The reference for this implementation is from appendix D in https://arxiv.org/pdf/1602.03184.pdf, originally derived in https://doi.org/10.1007/BF00165988 */

    memset(Transfer_Operator, 0, 16 * sizeof(double));
    memset(Integrated_Transfer_Operator, 0, 16 * sizeof(double));

    // Here I define a bunch of references, because its going to get hairy if I don't...
    auto& alpha = Absorbtion_functions;
    auto& rho = Faradey_functions;

    /* These are the variables defined in D8 - D13, used in calculating the M matricies */
    double const alpha_squared = alpha[Q] * alpha[Q] +
                                 alpha[U] * alpha[U] +
                                 alpha[V] * alpha[V];

    double const rho_squared = rho[Q] * rho[Q] +
                               rho[U] * rho[U] +
                               rho[V] * rho[V];

    double const alpha_rho = alpha[Q] * rho[Q] +
                             alpha[U] * rho[U] +
                             alpha[V] * rho[V];

    // sigma is the sign of the variable alpha_rho
    double const sigma = copysign(1.0, alpha_rho);

    // These quantities can go ever so sligtly negative, which physically should not happen, but nmerically it does.
    // This breaks the sqrt() functions, and so guards have to be put in place
    double const Theta = 2 * std::sqrt((alpha_squared - rho_squared) * (alpha_squared - rho_squared) / 4 + alpha_rho * alpha_rho);

    if (isnan(Theta)) { throw std::runtime_error("Inavlid value for the Theta coefficient in Get_radiative_transfer_operators()! \n"); }

    if (isinf(1 / Theta)) { 

        /* We only end up in here when absorbtion and faradey rotation are negligable. */

        double const M_1[4][4] = { {1.0, 0.0, 0.0, 0.0},
                                   {0.0, 1.0, 0.0, 0.0},
                                   {0.0, 0.0, 1.0, 0.0},
                                   {0.0, 0.0, 0.0, 1.0} };
        
        const double exp_term = exp(-alpha[I] * CGS_Step);

        for (int row_idx = 0; row_idx < e_Stokes_param_num; row_idx++) {

            for (int colum_idx = 0; colum_idx < e_Stokes_param_num; colum_idx++) {

                Transfer_Operator[row_idx][colum_idx] = M_1[row_idx][colum_idx] * exp_term;

                if (not isnan(1 / alpha[I]) and not isinf(1 / alpha[I])) {

                    Integrated_Transfer_Operator[row_idx][colum_idx] = M_1[row_idx][colum_idx] / alpha[I] * (1. - exp_term);

                }
                else {

                    /* ------- This is the zero absorbtion case ------- */

                    Integrated_Transfer_Operator[row_idx][colum_idx] = M_1[row_idx][colum_idx] * CGS_Step;

                }

            }

        }
        
        return; 
    
    }

    double const Lambda[2] = { std::sqrt((Theta / 2 + (alpha_squared - rho_squared) / 2)),
                               std::sqrt((Theta / 2 - (alpha_squared - rho_squared) / 2)) };

    if (isnan(Lambda[0]) or isnan(Lambda[1])) {

        throw std::runtime_error("Inavlid value for the Lambda coefficient in Get_radiative_transfer_operators()! \n");

    }

    /* Thesse are used in the "scaling factors" infront of the M matricies */
    double exp_I = exp(-alpha[I] * CGS_Step);

    /* At high optical depths, the exp * sinh and exp * cosh terms end up evaluating to 0 * inf, which breaks the code.
       this problem is solved by substituting in the exponential definition of the hyberbaulic functions and analytically
       multipling by the other exponential. */
    double exp_sinh = (exp((-alpha[I] + Lambda[0]) * CGS_Step) - exp((-alpha[I] - Lambda[0]) * CGS_Step)) / 2;
    double exp_cosh = (exp((-alpha[I] + Lambda[0]) * CGS_Step) + exp((-alpha[I] - Lambda[0]) * CGS_Step)) / 2;

    double const cos_term = cos(Lambda[1] * CGS_Step);
    double const sin_term = sin(Lambda[1] * CGS_Step);

    /* ========================== M_1 Matrix calculation ========================== */

    double const M_1_scale_factor = (exp_cosh + exp_I * cos_term) / 2;

    double const M_1[4][4] = { {1.0, 0.0, 0.0, 0.0},
                               {0.0, 1.0, 0.0, 0.0},
                               {0.0, 0.0, 1.0, 0.0},
                               {0.0, 0.0, 0.0, 1.0} };

    /* ========================== M_2 Matrix calculation ========================== */

    const double M_2_scale_factor = -exp_I * sin_term;

    if (isinf(M_2_scale_factor) or isnan(M_2_scale_factor)) {

        throw std::runtime_error("Inavlid value for M_2_scale_factor in Get_radiative_transfer_operators()!");

    }

    double M_2[4][4] = { {                         0,                          (Lambda[1] * alpha[Q] - sigma * Lambda[0] * rho[Q]),  ( Lambda[1] * alpha[U] - sigma * Lambda[0] * rho[U]), ( Lambda[1] * alpha[V] - sigma * Lambda[0] * rho[V])},
                         {(Lambda[1] * alpha[Q] - sigma * Lambda[0] * rho[Q]),                           0,                          ( sigma * Lambda[0] * alpha[V] + Lambda[1] * rho[V]), (-sigma * Lambda[0] * alpha[U] - Lambda[1] * rho[U])},
                         {(Lambda[1] * alpha[U] - sigma * Lambda[0] * rho[U]), (-sigma * Lambda[0] * alpha[V] - Lambda[1] * rho[V]),                           0,                          ( sigma * Lambda[0] * alpha[Q] + Lambda[1] * rho[Q])},
                         {(Lambda[1] * alpha[V] - sigma * Lambda[0] * rho[V]), ( sigma * Lambda[0] * alpha[U] + Lambda[1] * rho[U]), (-sigma * Lambda[0] * alpha[Q] - Lambda[1] * rho[Q]),                           0                         } };

    /* ========================== M_3 Matrix calculation ========================== */

    double const M_3_scale_factor = -exp_sinh;

    if (isinf(M_3_scale_factor) or isnan(M_3_scale_factor)) {

        throw std::runtime_error("Inavlid value for M_3_scale_factor in Get_radiative_transfer_operators()!");

    }

    double M_3[4][4] = { {						 0,							   ( Lambda[0] * alpha[Q] + sigma * Lambda[1] * rho[Q]), ( Lambda[0] * alpha[U] + sigma * Lambda[1] * rho[Q]),   ( Lambda[0] * alpha[V] + sigma * Lambda[1] * rho[V])},
                         {(Lambda[0] * alpha[Q] + sigma * Lambda[1] * rho[Q]),	 		               0,                            (-sigma * Lambda[1] * alpha[V] + Lambda[0] * rho[V]), ( sigma * Lambda[1] * alpha[U] - Lambda[0] * rho[U])},
                         {(Lambda[0] * alpha[U] + sigma * Lambda[1] * rho[U]), ( sigma * Lambda[1] * alpha[V] - Lambda[0] * rho[V]),	                         0,	                         (-sigma * Lambda[1] * alpha[Q] + Lambda[0] * rho[Q])},
                         {(Lambda[0] * alpha[V] + sigma * Lambda[1] * rho[V]), (-sigma * Lambda[1] * alpha[U] + Lambda[0] * rho[U]), ( sigma * Lambda[1] * alpha[Q] - Lambda[0] * rho[Q]),                           0                         } };

    /* ========================== M_4 Matrix calculation ========================== */

    double const M_4_scale_factor = (exp_cosh - exp_I * cos_term) / 2;

    if (isinf(M_4_scale_factor) or isnan(M_4_scale_factor)) {

        throw std::runtime_error("Inavlid value for M_4_scale_factor in Get_radiative_transfer_operators()!");

    }

    double M_4[4][4] = { {   (alpha_squared + rho_squared) / 2,                      (alpha[V] * rho[U] - alpha[U] * rho[V]),                                     (alpha[Q] * rho[V] - alpha[V] * rho[Q]),                                     (alpha[U] * rho[Q] - alpha[Q] * rho[U])},
                         {(alpha[U] * rho[V] - alpha[V] * rho[U]), (alpha[Q] * alpha[Q] + rho[Q] * rho[Q] - (alpha_squared + rho_squared) / 2),                   (alpha[Q] * alpha[U] + rho[Q] * rho[U]),                                     (alpha[V] * alpha[Q] + rho[V] * rho[Q])},
                         {(alpha[V] * rho[Q] - alpha[Q] * rho[V]),                   (alpha[Q] * alpha[U] + rho[Q] * rho[U]),                   (alpha[U] * alpha[U] + rho[U] * rho[U] - (alpha_squared + rho_squared) / 2),                   (alpha[U] * alpha[V] + rho[U] * rho[V])},
                         {(alpha[Q] * rho[U] - alpha[U] * rho[Q]),                   (alpha[V] * alpha[Q] + rho[V] * rho[Q]),                                     (alpha[U] * alpha[V] + rho[U] * rho[V]),                   (alpha[V] * alpha[V] + rho[V] * rho[V] - (alpha_squared + rho_squared) / 2)} };

    for (int left_idx = 0; left_idx < 4; left_idx++) {

        for (int right_idx = 0; right_idx < 4; right_idx++) {

            M_2[left_idx][right_idx] /= Theta;
            M_3[left_idx][right_idx] /= Theta;
            M_4[left_idx][right_idx] /= Theta / 2;

        }

    }

    /* ========================== This is the formal operator O(s,s') - the solution to D1 ========================== */

    for (int row_idx = 0; row_idx < e_Stokes_param_num; row_idx++) {

        for (int colum_idx = 0; colum_idx < e_Stokes_param_num; colum_idx++) {

            Transfer_Operator[row_idx][colum_idx] = M_1_scale_factor * M_1[row_idx][colum_idx] +
                                                    M_2_scale_factor * M_2[row_idx][colum_idx] +
                                                    M_3_scale_factor * M_3[row_idx][colum_idx] +
                                                    M_4_scale_factor * M_4[row_idx][colum_idx];

        }

    }

    /* ========================== The intergral of O(s,s') for constant M matricies ========================== */

    /* This part of the implementation is adapted from equation (24) of https://academic.oup.com/mnras/article/475/1/43/4712230 */

    double const f_1 = 1.0 / (alpha[I] * alpha[I] - Lambda[0] * Lambda[0]);
    double const f_2 = 1.0 / (alpha[I] * alpha[I] + Lambda[1] * Lambda[1]);

    if (isinf(f_1) or isnan(f_1) or isinf(f_2) or isnan(f_2)) {

        throw std::runtime_error("Inavlid values for f_1 and f_2 in Get_radiative_transfer_operators()!");

    }

    for (int row_idx = 0; row_idx < e_Stokes_param_num; row_idx++) {

        for (int colum_idx = 0; colum_idx < e_Stokes_param_num; colum_idx++) {


            Integrated_Transfer_Operator[row_idx][colum_idx] = -Lambda[0] * f_1 * M_3[row_idx][colum_idx] + alpha[I] * f_1 / 2 * (M_1[row_idx][colum_idx] + M_4[row_idx][colum_idx]) +
                                                                Lambda[1] * f_2 * M_2[row_idx][colum_idx] + alpha[I] * f_2 / 2 * (M_1[row_idx][colum_idx] - M_4[row_idx][colum_idx]) -
                                                                ((-Lambda[0] * f_1 * M_3[row_idx][colum_idx] + alpha[I]  * f_1 / 2 * (M_1[row_idx][colum_idx] + M_4[row_idx][colum_idx])) * exp_cosh +
                                                                 (-Lambda[1] * f_2 * M_2[row_idx][colum_idx] + alpha[I]  * f_2 / 2 * (M_1[row_idx][colum_idx] - M_4[row_idx][colum_idx])) * exp_I * cos_term +
                                                                 ( -alpha[I] * f_2 * M_2[row_idx][colum_idx] - Lambda[1] * f_2 / 2 * (M_1[row_idx][colum_idx] - M_4[row_idx][colum_idx])) * exp_I * sin_term -
                                                                 (  alpha[I] * f_1 * M_3[row_idx][colum_idx] - Lambda[0] * f_1 / 2 * (M_1[row_idx][colum_idx] + M_4[row_idx][colum_idx])) * exp_sinh);

        }

    }

}

void Emission_Integrator_class::Run_Analytic_Stokes_Vector_Propagator(const double Start_Affine_Param, const double End_Affine_Param) {

    Transfer_functions_type Total_Transfer_Functions{};

    for (int emission_medium = Disk; emission_medium < e_Emission_medium_number; emission_medium++) {

        Transfer_functions_type Temp_Transfer_functions{};

        this->p_Sim_Context->p_Emission_Model->get_radiative_transfer_functions(this->get_ray_Local_State_Vector(Start_Affine_Param),
                                                                                static_cast<Emission_medium_enums>(emission_medium),
                                                                                &Temp_Transfer_functions);

        add_vectors(Temp_Transfer_functions.Absorbtion_functions, Total_Transfer_Functions.Absorbtion_functions, e_Stokes_param_num, Total_Transfer_Functions.Absorbtion_functions);
        add_vectors(Temp_Transfer_functions.Emission_functions, Total_Transfer_Functions.Emission_functions, e_Stokes_param_num, Total_Transfer_Functions.Emission_functions);
        add_vectors(Temp_Transfer_functions.Faradey_functions, Total_Transfer_Functions.Faradey_functions, e_Stokes_param_num, Total_Transfer_Functions.Faradey_functions);

    }

    const double Geometric_Step = std::abs(End_Affine_Param - Start_Affine_Param);
    const double CGS_Step = Geometric_Step * MASS_TO_CM * this->p_Sim_Context->p_Init_Conditions->central_object_mass;

    double Transfer_operator[4][4]{};
    double Integrated_transfer_operator[4][4];

    this->Get_radiative_transfer_operators(Total_Transfer_Functions.Absorbtion_functions, Total_Transfer_Functions.Faradey_functions, CGS_Step, Transfer_operator, Integrated_transfer_operator);

    double Transfered_emission_vector[e_Stokes_param_num]{};

    // Placeholder vector for use in the mat_vec_multiply_4D() function
    double Temp_Intensity[e_Stokes_param_num]{};
    memcpy(Temp_Intensity, this->Current_Stokes_Vector, e_Stokes_param_num * sizeof(double));

    mat_vec_multiply_4D(Integrated_transfer_operator, Total_Transfer_Functions.Emission_functions, Transfered_emission_vector);
    mat_vec_multiply_4D(Transfer_operator, Temp_Intensity, this->Current_Stokes_Vector);

    for (int index = 0; index < e_Stokes_param_num; index++) {

        this->Current_Stokes_Vector[index] += Transfered_emission_vector[index];

    }

    this->Current_Affine_Param += Geometric_Step;
    this->Update_emission_log();

}