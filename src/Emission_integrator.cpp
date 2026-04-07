#include "Emission_integrator.h"
#include "Emission_Models.h"

Emission_Integrator_class::Emission_Integrator_class(const Simulation_Context_type* p_Sim_Context, Results_type* const p_Ray_results) {

    this->e_Active_integrator = p_Sim_Context->p_Init_Conditions->Integrator_params.e_Radiative_transfer_integrator;

    this->p_Sim_Context = p_Sim_Context;
    this->Ray_log_length = p_Ray_results->Ray_log_struct.Log_length;

    for (int idx = I; idx < e_Stokes_param_num; idx ++) {

        this->Emission_log[idx] = p_Ray_results->Ray_log_struct.Ray_emission_log[idx];

    }

    for (int idx = e_x; idx <= e_y; idx++) {

        this->Polarization_log[idx] = p_Ray_results->Ray_log_struct.Ray_polarization_log[idx];

    }

    /* --------------------------------------------------- Init the counters -------------------------------------------------- */

    this->Current_emission_log_idx = 0;
    this->Current_polarization_log_idx = 0;

    /* ------------------------------------------ Construct the geodesic ray splines ------------------------------------------ */

    double* Ray_Log[e_Dynamic_state_size]{};

    double* Radial_Ray_log_global[2]{};
    
    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        Ray_Log[idx] = new double[p_Ray_results->Ray_log_struct.Log_length];

    }

    Radial_Ray_log_global[0] = new double[p_Ray_results->Ray_log_struct.Log_length];
    Radial_Ray_log_global[1] = new double[p_Ray_results->Ray_log_struct.Log_length];

    this->Affine_param_log = new double[p_Ray_results->Ray_log_struct.Log_length];

    for (int log_idx = 0; log_idx < p_Ray_results->Ray_log_struct.Log_length; log_idx++) {

        for (int component_idx = 0; component_idx < e_Dynamic_state_size; component_idx++) {

            Ray_Log[component_idx][(p_Ray_results->Ray_log_struct.Log_length - 1) - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_local[component_idx + log_idx * e_Full_state_size];

        }

        Radial_Ray_log_global[0][(p_Ray_results->Ray_log_struct.Log_length - 1) - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_global[e_r + log_idx * e_Full_state_size];
        Radial_Ray_log_global[1][(p_Ray_results->Ray_log_struct.Log_length - 1) - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_global[e_p_r + log_idx * e_Full_state_size];

        this->Affine_param_log[(p_Ray_results->Ray_log_struct.Log_length - 1) - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_local[e_ray_affine_param + log_idx * e_Full_state_size];
    }

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        this->p_Ray_spline_instance[idx] = gsl_spline_alloc(gsl_interp_cspline, p_Ray_results->Ray_log_struct.Log_length);
        gsl_spline_init(this->p_Ray_spline_instance[idx], this->Affine_param_log, Ray_Log[idx], p_Ray_results->Ray_log_struct.Log_length);
        
        delete Ray_Log[idx];

    }

    /* ---------------------------------------- Radial spline in global coordinates ---------------------------------------- */

    this->p_Ray_spline_instance[e_Dynamic_state_size + 0] = gsl_spline_alloc(gsl_interp_cspline, p_Ray_results->Ray_log_struct.Log_length);
    gsl_spline_init(this->p_Ray_spline_instance[e_Dynamic_state_size + 0], this->Affine_param_log, Radial_Ray_log_global[0], p_Ray_results->Ray_log_struct.Log_length);
    delete Radial_Ray_log_global[0];

    this->p_Ray_spline_instance[e_Dynamic_state_size + 1] = gsl_spline_alloc(gsl_interp_cspline, p_Ray_results->Ray_log_struct.Log_length);
    gsl_spline_init(this->p_Ray_spline_instance[e_Dynamic_state_size + 1], this->Affine_param_log, Radial_Ray_log_global[1], p_Ray_results->Ray_log_struct.Log_length);
    delete Radial_Ray_log_global[1];

    this->p_Spline_accelerator = gsl_interp_accel_alloc();

    delete this->Affine_param_log;

}

Emission_Integrator_class::~Emission_Integrator_class() {

    for (int idx = 0; idx < e_Dynamic_state_size + 2; idx++) {

        gsl_spline_free(this->p_Ray_spline_instance[idx]);

    }

    gsl_interp_accel_free(this->p_Spline_accelerator);

}

void Emission_Integrator_class::Propagate_Stokes_Vector(const double Start_Affine_Param, const double End_Affine_Param) {

    if (RK78_Fehlberg != e_Active_integrator and RK78_DP != e_Active_integrator and RK54 != e_Active_integrator) {

        throw std::runtime_error("Wrong active integrator in Emission_Integrator_class::Run_Explicit_Runge_Kutta()!");

    }

    memset(this->Rad_Transfer_RHS_log, 0, sizeof(double) * RK78_size * e_Stokes_param_num);

    int RK_size = RK78_size;
    auto Stage_coeff = this->RK78_DP_Coeff_deriv;
    auto Main_solution_coeff = this->RK78_DP_Coeff_sol_main;
    auto Embedded_solution_coeff = this->RK78_DP_Coeff_sol_embeded;
    auto Affine_param_coeff = this->RK78_DP_Coeff_affine_param;

    if (RK78_Fehlberg == this->e_Active_integrator) {

        Stage_coeff = this->RK78_Fhelberg_Coeff_deriv;
        Main_solution_coeff = this->RK78_Fhelberg_Coeff_sol_main;
        Embedded_solution_coeff = this->RK78_Fhelberg_Coeff_sol_embeded;
        Affine_param_coeff = this->RK78_Fhelberg_Coeff_affine_param;

    }
    else if (RK54 == this->e_Active_integrator) {

        Stage_coeff = this->RK54_Coeff_deriv;
        Main_solution_coeff = this->RK54_Coeff_sol_main;
        Embedded_solution_coeff = this->RK54_Coeff_test_embeded;
        Affine_param_coeff = this->RK54_Coeff_affine_param;

        RK_size = RK54_size;

    }

    double New_Stokes_vector[e_Stokes_param_num]{};
    double Temp_Stokes_Vector[e_Stokes_param_num]{};
    double Temp_affine_param{};

    const double Geometric_Step = End_Affine_Param - Start_Affine_Param;
    const double CGS_Step = Geometric_Step * MASS_TO_CM * this->p_Sim_Context->p_Init_Conditions->central_object_mass;
    this->Current_Affine_Param = Start_Affine_Param;

    for (int RK_stage = 0; RK_stage < RK_size; RK_stage++) {

        memcpy(Temp_Stokes_Vector, this->Current_Stokes_Vector, e_Stokes_param_num * sizeof(double));

        for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {

            for (int derivative_idx = 0; derivative_idx < RK_stage; derivative_idx++) {

                Temp_Stokes_Vector[state_idx] += CGS_Step * Stage_coeff[RK_stage][derivative_idx] * this->Rad_Transfer_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];

            }
        }

        Temp_affine_param = this->Current_Affine_Param + Affine_param_coeff[RK_stage] * Geometric_Step;

        /* ------------------------------------------------------------ Get the radiative transfer RHS ------------------------------------------------------------ */

        Transfer_functions_type Total_Transfer_Functions{};

        for (int emission_medium = Disk; emission_medium < e_Emission_medium_number; emission_medium++) {

            Transfer_functions_type Temp_Transfer_functions{};

            this->p_Sim_Context->p_Emission_Model->get_radiative_transfer_functions(this->get_ray_Local_State_Vector(Temp_affine_param),
                                                                                    this->p_Sim_Context,
                                                                                    static_cast<Emission_medium_enums>(emission_medium),
                                                                                    &Temp_Transfer_functions);

            add_vectors(Temp_Transfer_functions.Absorbtion_functions, Total_Transfer_Functions.Absorbtion_functions, e_Stokes_param_num, Total_Transfer_Functions.Absorbtion_functions);
            add_vectors(Temp_Transfer_functions.Emission_functions, Total_Transfer_Functions.Emission_functions, e_Stokes_param_num, Total_Transfer_Functions.Emission_functions);
            add_vectors(Temp_Transfer_functions.Faradey_functions, Total_Transfer_Functions.Faradey_functions, e_Stokes_param_num, Total_Transfer_Functions.Faradey_functions);

            
        }

        this->get_Radiative_transfer_RHS(Total_Transfer_Functions.Emission_functions,
                                         Total_Transfer_Functions.Absorbtion_functions,
                                         Total_Transfer_Functions.Faradey_functions,
                                         Temp_Stokes_Vector,
                                         this->Rad_Transfer_RHS_log + RK_stage * e_Stokes_param_num);

        /* ------------------------------------------------------------------------------------------------------------------------------------------------------- */

    }

    memcpy(New_Stokes_vector, this->Current_Stokes_Vector, e_Stokes_param_num * sizeof(double));

    for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {

        for (int derivative_idx = 0; derivative_idx < RK_size; derivative_idx++) {

            New_Stokes_vector[state_idx] += CGS_Step * Main_solution_coeff[derivative_idx] * this->Rad_Transfer_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];

        }

    }

    memcpy(this->Current_Stokes_Vector, New_Stokes_vector, e_Stokes_param_num * sizeof(double));
    this->Current_Affine_Param += Geometric_Step;

    this->Update_emission_log(New_Stokes_vector);

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

void Emission_Integrator_class::Update_emission_log(const double* const New_Stokes_Vector) {

    this->Current_emission_log_idx += 1;

    for (int idx = 0; idx < e_Stokes_param_num; idx++) {

        this->Emission_log[idx][this->Current_emission_log_idx] = New_Stokes_Vector[idx];

    }

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

const double* const Emission_Integrator_class::get_current_Stokes_Vector() const {

    return this->Current_Stokes_Vector;

}

void Emission_Integrator_class::Propagate_Polarization_Vector(const double Start_Affine_Param, const double End_Affine_Param, const Tensor_type_enums e_Vec_type) {

    if (RK78_Fehlberg != e_Active_integrator and RK78_DP != e_Active_integrator and RK54 != e_Active_integrator) {

        throw std::runtime_error("Wrong active integrator in Emission_Integrator_class::Propagate_Polarization_Vector()!");

    }

    memset(this->Parallel_Transport_RHS_log, 0, sizeof(std::complex<double>) * RK78_size * e_Stokes_param_num);

    int RK_size = RK78_size;
    auto Stage_coeff = this->RK78_DP_Coeff_deriv;
    auto Main_solution_coeff = this->RK78_DP_Coeff_sol_main;
    auto Embedded_solution_coeff = this->RK78_DP_Coeff_sol_embeded;
    auto Affine_param_coeff = this->RK78_DP_Coeff_affine_param;

    if (RK78_Fehlberg == this->e_Active_integrator) {

        Stage_coeff = this->RK78_Fhelberg_Coeff_deriv;
        Main_solution_coeff = this->RK78_Fhelberg_Coeff_sol_main;
        Embedded_solution_coeff = this->RK78_Fhelberg_Coeff_sol_embeded;
        Affine_param_coeff = this->RK78_Fhelberg_Coeff_affine_param;

    }
    else if (RK54 == this->e_Active_integrator) {

        Stage_coeff = this->RK54_Coeff_deriv;
        Main_solution_coeff = this->RK54_Coeff_sol_main;
        Embedded_solution_coeff = this->RK54_Coeff_test_embeded;
        Affine_param_coeff = this->RK54_Coeff_affine_param;

        RK_size = RK54_size;

    }

    std::complex<double> New_Pol_Vector[e_Stokes_param_num]{};
    std::complex<double> Temp_Pol_Vector[e_Stokes_param_num]{};
    double Temp_affine_param{};

    const double Step = End_Affine_Param - Start_Affine_Param;
    this->Current_Affine_Param = Start_Affine_Param;

    for (int RK_stage = 0; RK_stage < RK_size; RK_stage++) {
    
        memcpy(Temp_Pol_Vector, this->Current_Pol_Vector, e_Stokes_param_num * sizeof(std::complex<double>));
    
        for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {
    
            for (int derivative_idx = 0; derivative_idx < RK_stage; derivative_idx++) {
    
                Temp_Pol_Vector[state_idx] += Step * Stage_coeff[RK_stage][derivative_idx] * this->Parallel_Transport_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];
    
            }
        }
    
        Temp_affine_param = this->Current_Affine_Param + Affine_param_coeff[RK_stage] * Step;
    
        /* --------------------------------------- Get the Parallel Transport RHS --------------------------------------- */
    
        this->get_Parallel_Transport_RHS(this->get_ray_Global_State_Vector(Temp_affine_param),
                                         Temp_Pol_Vector, 
                                         e_Vec_type, 
                                         this->Parallel_Transport_RHS_log + RK_stage * e_Stokes_param_num);
    
        /* -------------------------------------------------------------------------------------------------------------- */
    
    }
    
    memcpy(New_Pol_Vector, this->Current_Pol_Vector, e_Stokes_param_num * sizeof(std::complex<double>));

    for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {
    
        for (int derivative_idx = 0; derivative_idx < RK_size; derivative_idx++) {
    
            New_Pol_Vector[state_idx] += Step * Main_solution_coeff[derivative_idx] * this->Parallel_Transport_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];
    
        }
    
    }
    
    memcpy(this->Current_Pol_Vector, New_Pol_Vector, e_Stokes_param_num * sizeof(std::complex<double>));
    this->Current_Affine_Param += Step;

    this->Update_polarization_log(New_Pol_Vector, End_Affine_Param);

}

void Emission_Integrator_class::Map_Stokes_to_Polarization_Vector(const double Stokes_Tetrad[e_Stokes_param_num][e_Stokes_param_num]) {

    memset(this->Current_Pol_Vector, 0, 4 * sizeof(std::complex<double>));

    std::complex<double> Stokes_Basis_Pol_vec[e_Stokes_param_num];

    double Polarized_Intensity = vector_norm(this->Current_Stokes_Vector + 1, 3);

    Stokes_Basis_Pol_vec[1] = M_SQRT1_2;

    if (!isinf(this->Current_Stokes_Vector[Q] / Polarized_Intensity) and !isnan(this->Current_Stokes_Vector[Q] / Polarized_Intensity)) {

        Stokes_Basis_Pol_vec[1] = sqrt((1 + this->Current_Stokes_Vector[Q] / Polarized_Intensity) / 2);

    }

    Stokes_Basis_Pol_vec[2] = 1.;

    if (!isinf(1. / std::norm(Stokes_Basis_Pol_vec[1] * Polarized_Intensity)) and !isnan(1. / std::norm(Stokes_Basis_Pol_vec[1] * Polarized_Intensity))) {

        Stokes_Basis_Pol_vec[2] = (this->Current_Stokes_Vector[U] - complex_i * this->Current_Stokes_Vector[V]) / (2.0 * Stokes_Basis_Pol_vec[1] * Polarized_Intensity);

    }

    for (int coord_idx = 0; coord_idx < 4; coord_idx++) {

        for (int stokes_idx = 0; stokes_idx < 4; stokes_idx++) {

            this->Current_Pol_Vector[coord_idx] += Stokes_Tetrad[stokes_idx][coord_idx] * Stokes_Basis_Pol_vec[stokes_idx];

        }

    }

}

void Emission_Integrator_class::Map_Polarization_Vector_to_Stokes(const double inv_Stokes_Tetrad[e_Stokes_param_num][e_Stokes_param_num]) {

    std::complex<double> Stokes_Basis_Pol_vec[4]{};
    
    for (int stokes_idx = 0; stokes_idx < 4; stokes_idx++) {
    
        for (int coord_idx = 0; coord_idx < 4 ; coord_idx++) {
    
            Stokes_Basis_Pol_vec[stokes_idx] += inv_Stokes_Tetrad[stokes_idx][coord_idx] * this->Current_Pol_Vector[coord_idx];

        }

    }

    double Polarized_Intensity = vector_norm(this->Current_Stokes_Vector + 1, 3);

    this->Current_Stokes_Vector[Q] = (Polarized_Intensity * (Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[1]) -
                                                             Stokes_Basis_Pol_vec[2] * std::conj(Stokes_Basis_Pol_vec[2]))).real();
    
    this->Current_Stokes_Vector[U] = (Polarized_Intensity * (Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[2]) +
                                                             Stokes_Basis_Pol_vec[2] * std::conj(Stokes_Basis_Pol_vec[1]))).real();
    
    this->Current_Stokes_Vector[V] = (Polarized_Intensity * (Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[2]) - 
                                                             Stokes_Basis_Pol_vec[2] * std::conj(Stokes_Basis_Pol_vec[1]))).imag();

    double Polarized_Intensity_after_mapping = vector_norm(this->Current_Stokes_Vector + 1, 3);

    if (!isnan(1. / Polarized_Intensity_after_mapping) and !isinf(1.0 / Polarized_Intensity_after_mapping)) {

        for (int idx = 1; idx < 4; idx++) {

            this->Current_Stokes_Vector[idx] *= Polarized_Intensity / Polarized_Intensity_after_mapping;

        }

    }

}


void Emission_Integrator_class::Update_polarization_log(const std::complex<double>* const New_Polarization_Vector, const double New_affine_param) {

    this->Current_polarization_log_idx += 1;

    std::complex<double> Polarization_vector_ZAMO[4];
    Metric_type s_Metric = this->p_Sim_Context->p_Spacetime->get_global_metric(this->get_ray_Global_State_Vector(New_affine_param));

    Contravariant_coord_to_ZAMO(&s_Metric, New_Polarization_Vector, Polarization_vector_ZAMO);

    this->Polarization_log[e_x][this->Current_polarization_log_idx] = Polarization_vector_ZAMO[e_phi].real();
    this->Polarization_log[e_y][this->Current_polarization_log_idx] = Polarization_vector_ZAMO[e_theta].real();
}


const std::complex<double>* const Emission_Integrator_class::get_current_Polarization_Vector() const {

    return this->Current_Pol_Vector;

};