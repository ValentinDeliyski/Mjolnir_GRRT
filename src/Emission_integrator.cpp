#include "Emission_integrator.h"
#include "Emission_Models.h"

Emission_Integrator_class::Emission_Integrator_class(const Simulation_Context_type* p_Sim_Context, Results_type* const p_Ray_results) {

    this->e_Active_integrator = p_Sim_Context->p_Init_Conditions->Integrator_params.e_Radiative_transfer_integrator;

    this->p_Step_controller = new Step_controller_class(p_Sim_Context->p_Init_Conditions->Integrator_params.Rad_Transfer_Step_Controller_Params);

    this->p_Sim_Context = p_Sim_Context;
    this->p_Emission_Model = p_Sim_Context->p_Emission_Model;
    
    this->Ray_log_length = p_Ray_results->Ray_log_struct.Log_length;
    this->Emission_log = p_Ray_results->Ray_log_struct.Ray_emission_log;

    /* ----------------------------------------------- Init the flags / counters ---------------------------------------------- */

    this->Current_log_idx = 0;
    this->N_steps_rejected = 0;
    this->NaN_checker_count = 0;
    this->continue_integration = true;

    /* ------------------------------------------ Construct the geodesic ray splines ------------------------------------------ */

    double* Ray_Log[e_Dynamic_state_size]{};

    double* Radial_Ray_log_global[2]{};
    
    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        Ray_Log[idx] = new double[p_Ray_results->Ray_log_struct.Log_length];

    }

    Radial_Ray_log_global[0] = new double[p_Ray_results->Ray_log_struct.Log_length];
    Radial_Ray_log_global[1] = new double[p_Ray_results->Ray_log_struct.Log_length];

    this->Affine_param_log = new double[p_Ray_results->Ray_log_struct.Log_length];
    this->Geodesic_step_log = new double[p_Ray_results->Ray_log_struct.Log_length];

    for (int log_idx = 0; log_idx < p_Ray_results->Ray_log_struct.Log_length; log_idx++) {

        for (int component_idx = 0; component_idx < e_Dynamic_state_size; component_idx++) {

            Ray_Log[component_idx][(p_Ray_results->Ray_log_struct.Log_length - 1) - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_local[component_idx + log_idx * e_Full_state_size];

        }

        Radial_Ray_log_global[0][(p_Ray_results->Ray_log_struct.Log_length - 1) - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_global[e_r + log_idx * e_Full_state_size];
        Radial_Ray_log_global[1][(p_Ray_results->Ray_log_struct.Log_length - 1) - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_global[e_p_r + log_idx * e_Full_state_size];

        this->Affine_param_log[(p_Ray_results->Ray_log_struct.Log_length - 1) - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_local[e_ray_affine_param + log_idx * e_Full_state_size];
        this->Geodesic_step_log[(p_Ray_results->Ray_log_struct.Log_length - 1) - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_local[e_step + log_idx * e_Full_state_size];
    }

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        this->p_Ray_spline_instance[idx] = gsl_spline_alloc(gsl_interp_cspline, p_Ray_results->Ray_log_struct.Log_length);
        gsl_spline_init(this->p_Ray_spline_instance[idx], this->Affine_param_log, Ray_Log[idx], p_Ray_results->Ray_log_struct.Log_length);
        
        delete[] Ray_Log[idx];

    }

    /* ---------------------------------------- Radial spline in global coordinates ---------------------------------------- */

    this->p_Ray_spline_instance[e_Dynamic_state_size + 0] = gsl_spline_alloc(gsl_interp_cspline, p_Ray_results->Ray_log_struct.Log_length);
    gsl_spline_init(this->p_Ray_spline_instance[e_Dynamic_state_size + 0], this->Affine_param_log, Radial_Ray_log_global[0], p_Ray_results->Ray_log_struct.Log_length);

    this->p_Ray_spline_instance[e_Dynamic_state_size + 1] = gsl_spline_alloc(gsl_interp_cspline, p_Ray_results->Ray_log_struct.Log_length);
    gsl_spline_init(this->p_Ray_spline_instance[e_Dynamic_state_size + 1], this->Affine_param_log, Radial_Ray_log_global[1], p_Ray_results->Ray_log_struct.Log_length);

    this->p_Spline_accelerator = gsl_interp_accel_alloc();

}

Emission_Integrator_class::~Emission_Integrator_class() {

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        gsl_spline_free(this->p_Ray_spline_instance[idx]);

    }

    delete[] this->Affine_param_log;
    delete[] this->Geodesic_step_log;

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

    double state_error[e_Stokes_param_num]{};
    double Temp_Stokes_Vector[e_Stokes_param_num]{};
    double Temp_affine_param{};

    double New_Stokes_vector_main[e_Stokes_param_num]{};
    double New_Stokes_vector_embeded[e_Stokes_param_num]{};

    /* ------ Set the initial step to the one used by the geodesic integrator ------ */
    this->p_Step_controller->step = End_Affine_Param - Start_Affine_Param;
    this->Current_Affine_Param = Start_Affine_Param;

    /* -------- The small offset is here because I have a check below that ensures the ray doesnt overshoot the end affine parameter. As a result it only reaches it when it hits the 
                precision limit on the double's holding Current_affine_param and End_Affine_Param. This wastes compute time to integrate over infinitesimal intervals -------- */

    while (this->Current_Affine_Param < End_Affine_Param - this->End_affine_param_offset) {

        double Geometric_controller_step = this->p_Step_controller->step;
        double CGS_controller_step = Geometric_controller_step * MASS_TO_CM * this->p_Sim_Context->p_Init_Conditions->central_object_mass;

        /* ---------------- This ensures that the integrrator wont overshoot the target affine parameter and either break the splines by going out of bounds
                            or simply overshooting the emitting region, which makes the image visibly noisy ---------------- */

        if (this->Current_Affine_Param + Geometric_controller_step > End_Affine_Param) {

            Geometric_controller_step = 0.99 * std::abs(this->Current_Affine_Param - End_Affine_Param);
            CGS_controller_step = Geometric_controller_step * MASS_TO_CM * this->p_Sim_Context->p_Init_Conditions->central_object_mass;

        }

        for (int RK_stage = 0; RK_stage < RK_size; RK_stage++) {

            memcpy(Temp_Stokes_Vector, this->Current_Stokes_Vector, e_Stokes_param_num * sizeof(double));

            for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {

                for (int derivative_idx = 0; derivative_idx < RK_stage; derivative_idx++) {

                    Temp_Stokes_Vector[state_idx] += CGS_controller_step * Stage_coeff[RK_stage][derivative_idx] * this->Rad_Transfer_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];

                }
            }

            Temp_affine_param = this->Current_Affine_Param + Affine_param_coeff[RK_stage] * Geometric_controller_step;

            /* ------------------------------------------------------------ Get the radiative transfer RHS ------------------------------------------------------------ */

            Transfer_functions_type Total_Transfer_Functions{};

            for (int emission_medium = Disk; emission_medium <= Hotspot; emission_medium++) {

                Transfer_functions_type Temp_Transfer_functions{};

                this->p_Emission_Model->get_radiative_transfer_functions(this->get_ray_Local_State_Vector(Temp_affine_param),
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

        for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {

            New_Stokes_vector_main[state_idx] = this->Current_Stokes_Vector[state_idx];
            New_Stokes_vector_embeded[state_idx] = this->Current_Stokes_Vector[state_idx];

            for (int derivative_idx = 0; derivative_idx < RK_size; derivative_idx++) {

                New_Stokes_vector_main[state_idx] += CGS_controller_step * Main_solution_coeff[derivative_idx] * this->Rad_Transfer_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];
                New_Stokes_vector_embeded[state_idx] += CGS_controller_step * Embedded_solution_coeff[derivative_idx] * this->Rad_Transfer_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];

            }

            state_error[state_idx] = New_Stokes_vector_main[state_idx] - New_Stokes_vector_embeded[state_idx];

        }

        if (ERROR == this->Run_NaN_checker(New_Stokes_vector_main, New_Stokes_vector_embeded)) { continue; }

        this->p_Step_controller->update_state_errors(New_Stokes_vector_main, state_error, this->e_Active_integrator, e_Stokes_param_num);

        if (this->p_Step_controller->current_err < 1.0 or !this->p_Step_controller->Parameters.Use_adaptive_step) {

            this->continue_integration = true;

        }
        else {

            this->continue_integration = false;
            this->N_steps_rejected++;

        }

        if (this->continue_integration) {

            memcpy(this->Current_Stokes_Vector, New_Stokes_vector_main, e_Stokes_param_num * sizeof(double));
            this->Current_Affine_Param += Geometric_controller_step;
            this->Update_emission_log(New_Stokes_vector_main);
            //this->Update_debug_log();

            this->N_steps_rejected = 0;
            this->NaN_checker_count = 0;

        }

        this->p_Step_controller->update_step(this->e_Active_integrator);

    }

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

    this->Current_log_idx += 1;

    for (int idx = 0; idx < e_Stokes_param_num; idx++) {

        this->Emission_log[idx][this->Current_log_idx] = New_Stokes_Vector[idx];

    }

    this->Emission_log[e_Stokes_affine_param][this->Current_log_idx] = this->Current_Affine_Param;

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

void static Parallel_Transport_RHS(const double* const Global_State_Vector,
                                   std::complex<double>* Vector_to_transport,
                                   Spacetime_Base_Class* const p_Spacetime,
                                   Tensor_type_enums e_Vec_type,
    std::complex<double>* const Vector_derivative) {

    Metric_type s_Metric = p_Spacetime->get_global_metric(Global_State_Vector);
    Metric_type s_dr_Metric = p_Spacetime->get_dr_global_metric(Global_State_Vector);
    Metric_type s_dtheta_Metric = p_Spacetime->get_dtheta_global_metric(Global_State_Vector);

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

                    Vector_derivative[derivative_index] += -Connection_Coefficients[derivative_index][wave_vector_index][polarization_index] * Photon_momentum_contravariant[wave_vector_index] * Vector_to_transport[polarization_index];

                }
                else {

                    Vector_derivative[derivative_index] += Connection_Coefficients[polarization_index][derivative_index][wave_vector_index] * Photon_momentum_contravariant[wave_vector_index] * Vector_to_transport[polarization_index];

                }

            }

        }

    }

}

void static Parallel_Transport_Vector(const double* const State_Vector_Global,
                                      Spacetime_Base_Class* const p_Spacetime,
                                      Tensor_type_enums const e_Vec_type,
                                      std::complex<double>* const Vector_to_transport) {

    std::complex<double> RHS[Nyström_size * e_Stokes_param_num]{};
    double EOM[Nyström_size * e_Dynamic_state_size]{};

    std::complex<double> Temp_Vector[e_Stokes_param_num]{};
    double Temp_State_Vector[e_Dynamic_state_size]{};

    for (int RK5_stage = 0; RK5_stage < Nyström_size; RK5_stage++) {

        memcpy(Temp_Vector, Vector_to_transport, e_Stokes_param_num * sizeof(std::complex<double>));
        memcpy(Temp_State_Vector, State_Vector_Global, e_Dynamic_state_size * sizeof(double));

        for (int derivative_indexer = 0; derivative_indexer < RK5_stage; derivative_indexer++) {

            for (int idx = 0; idx < 4; idx++) {

                Temp_Vector[idx] += Nyström_Deriv_coeffs[RK5_stage][derivative_indexer] * RHS[idx + derivative_indexer * e_Stokes_param_num] * State_Vector_Global[e_step];

            }

            for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

                Temp_State_Vector[idx] += Nyström_Deriv_coeffs[RK5_stage][derivative_indexer] * EOM[idx + derivative_indexer * e_Dynamic_state_size] * State_Vector_Global[e_step];

            }
        }

        Parallel_Transport_RHS(Temp_State_Vector, Temp_Vector, p_Spacetime, e_Vec_type, RHS + RK5_stage * e_Stokes_param_num);
        p_Spacetime->get_EOM(Temp_State_Vector, EOM + RK5_stage * e_Dynamic_state_size);

    } 

    for (int idx = 0; idx < e_Stokes_param_num; idx++) {

        for (int deriv_idx = 0; deriv_idx < Nyström_size; deriv_idx++) {

            Vector_to_transport[idx] += State_Vector_Global[e_step] * Nyström_Coeff_sol[deriv_idx] * RHS[idx + deriv_idx * e_Stokes_param_num];

        }
    }
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

    std::complex<double> state_error[e_Stokes_param_num]{};
    std::complex<double> Temp_Pol_Vector[e_Stokes_param_num]{};
    double Temp_affine_param{};

    std::complex<double> New_Pol_Vector_main[e_Stokes_param_num]{};
    std::complex<double> New_Pol_Vector_embeded[e_Stokes_param_num]{};

    /* ------ Set the initial step to the one used by the geodesic integrator ------ */
    this->p_Step_controller->step = End_Affine_Param - Start_Affine_Param;
    this->Current_Affine_Param = Start_Affine_Param;

    /* -------- The small offset is here because I have a check below that ensures the ray doesnt overshoot the end affine parameter. As a result it only reaches it when it hits the
                precision limit on the double's holding Current_affine_param and End_Affine_Param. This wastes compute time to integrate over infinitesimal intervals -------- */

    while (this->Current_Affine_Param < End_Affine_Param - this->End_affine_param_offset) {

        /* ---------------- This ensures that the integrrator wont overshoot the target affine parameter and either break the splines by going out of bounds
                            or simply overshooting the emitting region, which makes the image visibly noisy ---------------- */

        if (this->Current_Affine_Param + this->p_Step_controller->step > End_Affine_Param) {

            this->p_Step_controller->step = 0.99 * std::abs(this->Current_Affine_Param - End_Affine_Param);

        }

        for (int RK_stage = 0; RK_stage < RK_size; RK_stage++) {

            memcpy(Temp_Pol_Vector, this->Current_Pol_Vector, e_Stokes_param_num * sizeof(std::complex<double>));

            for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {

                for (int derivative_idx = 0; derivative_idx < RK_stage; derivative_idx++) {

                    Temp_Pol_Vector[state_idx] += this->p_Step_controller->step * Stage_coeff[RK_stage][derivative_idx] * this->Parallel_Transport_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];

                }
            }

            Temp_affine_param = this->Current_Affine_Param + Affine_param_coeff[RK_stage] * this->p_Step_controller->step;

            /* --------------------------------------- Get the Parallel Transport RHS --------------------------------------- */

            this->get_Parallel_Transport_RHS(this->get_ray_Global_State_Vector(Temp_affine_param),
                                             Temp_Pol_Vector, 
                                             e_Vec_type, 
                                             this->Parallel_Transport_RHS_log + RK_stage * e_Stokes_param_num);

            /* -------------------------------------------------------------------------------------------------------------- */

        }

        for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {

            New_Pol_Vector_main[state_idx] = this->Current_Pol_Vector[state_idx];
            New_Pol_Vector_embeded[state_idx] = this->Current_Pol_Vector[state_idx];

            for (int derivative_idx = 0; derivative_idx < RK_size; derivative_idx++) {

                New_Pol_Vector_main[state_idx] += this->p_Step_controller->step * Main_solution_coeff[derivative_idx] * this->Parallel_Transport_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];
                New_Pol_Vector_embeded[state_idx] += this->p_Step_controller->step * Embedded_solution_coeff[derivative_idx] * this->Parallel_Transport_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];

            }

            state_error[state_idx] = New_Pol_Vector_main[state_idx] - New_Pol_Vector_embeded[state_idx];

        }

        if (ERROR == this->Run_NaN_checker(New_Pol_Vector_main, New_Pol_Vector_embeded)) { continue; }

        this->p_Step_controller->update_state_errors(New_Pol_Vector_main, state_error, this->e_Active_integrator, e_Stokes_param_num);

        if (this->p_Step_controller->current_err < 1.0 or !this->p_Step_controller->Parameters.Use_adaptive_step) {

            this->continue_integration = true;

        }
        else {

            this->continue_integration = false;
            this->N_steps_rejected++;

        }

        if (this->continue_integration) {

            memcpy(this->Current_Pol_Vector, New_Pol_Vector_main, e_Stokes_param_num * sizeof(std::complex<double>));
            this->Current_Affine_Param += this->p_Step_controller->step;

            this->N_steps_rejected = 0;
            this->NaN_checker_count = 0;

        }

        this->p_Step_controller->update_step(this->e_Active_integrator);

    }

}

void Emission_Integrator_class::Map_Stokes_to_Polarization_Vector(const double Stokes_Tetrad[e_Stokes_param_num][e_Stokes_param_num]) {

    memset(this->Current_Pol_Vector, 0, 4 * sizeof(std::complex<double>));

    std::complex<double> Stokes_Basis_Pol_vec[4];

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
    
    this->Current_Stokes_Vector[V] = (Polarized_Intensity * (Stokes_Basis_Pol_vec[2] * std::conj(Stokes_Basis_Pol_vec[1]) -
                                                             Stokes_Basis_Pol_vec[1] * std::conj(Stokes_Basis_Pol_vec[2]))).imag();

    double Polarized_Intensity_after_mapping = vector_norm(this->Current_Stokes_Vector + 1, 3);

    if (!isnan(1. / Polarized_Intensity_after_mapping) and !isinf(1.0 / Polarized_Intensity_after_mapping)) {

        for (int idx = 1; idx < 4; idx++) {

            this->Current_Stokes_Vector[idx] *= Polarized_Intensity / Polarized_Intensity_after_mapping;

        }

    }

}