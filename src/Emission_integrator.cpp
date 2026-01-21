#include "Emission_integrator.h"
#include "Emission_Models.h"


bool static Is_inside_emission_medium(const Simulation_Context_type* const p_Sim_Context,
                                        const double* const State_Vector_Local) {

    Emission_medium_state_type s_Hotspot_state{};
    Emission_medium_state_type s_Disk_state{};

    Return_Values Hotspot_velocity_OK = p_Sim_Context->p_Emission_Model->get_plasma_velocity(p_Sim_Context->p_Init_Conditions->Hotspot_params.Position,
        p_Sim_Context,
        p_Sim_Context->p_Init_Conditions->Hotspot_params.Velocity_profile_type,
        p_Sim_Context->p_Init_Conditions->Hotspot_params.Radial_velocity_fraction,
        s_Hotspot_state.Plasma_Velocity);
    bool In_hotspot = false;

    if (OK == Hotspot_velocity_OK) {

        In_hotspot = p_Sim_Context->p_Emission_Model->p_Hotspot_Model->is_inside_hotspot(State_Vector_Local,
            &s_Hotspot_state);
    }

    const bool In_disk = p_Sim_Context->p_Emission_Model->p_Disk_Model->is_inside_disk(State_Vector_Local,
        p_Sim_Context->p_Emission_Model->p_Disk_Model->s_Disk_params.e_Disk_model,
        &s_Disk_state);

    return In_hotspot || In_disk;

}

Emission_Integrator_class::Emission_Integrator_class(const Simulation_Context_type* p_Sim_Context, Results_type* const p_Ray_results) {

    this->Step_too_small = false;
    this->e_Active_integrator = p_Sim_Context->p_Init_Conditions->Integrator_params.e_Default_geodesic_integrator;

    this->p_Step_controller = new Step_controller_class(p_Sim_Context->p_Init_Conditions->Integrator_params);

    this->p_Sim_Context = p_Sim_Context;

    this->p_Emission_Model = p_Sim_Context->p_Emission_Model;

    memset(this-> Current_Stokes_Vector, 0, sizeof(double) * e_Stokes_param_num);

    /* ------------------------------ Construct the initial Stokes vector and init the log array ------------------------------ */

    this->Current_Affine_Param = p_Ray_results->Ray_log_struct.Ray_path_log_local[e_ray_affine_param + (p_Ray_results->Ray_log_struct.Log_length - 1) * e_Full_state_size];
    
    this->Ray_log_length = p_Ray_results->Ray_log_struct.Log_length;

    /* ---------------------------------------------------- Init the flags ---------------------------------------------------- */

    this->N_steps_rejected = 0;
    this->NaN_checker_count = 0;
    this->continue_integration = true;

    /* ------------------------------------------ Construct the geodesic ray splines ------------------------------------------ */

    double* Ray_Log[e_Dynamic_state_size]{};
    
    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        Ray_Log[idx] = new double[p_Ray_results->Ray_log_struct.Log_length];

    }

    this->Affine_param_log = new double[p_Ray_results->Ray_log_struct.Log_length];
    this->Geodesic_step_log = new double[p_Ray_results->Ray_log_struct.Log_length];

    for (int log_idx = 0; log_idx < p_Ray_results->Ray_log_struct.Log_length; log_idx++) {

        for (int component_idx = 0; component_idx < e_Dynamic_state_size; component_idx++) {

            Ray_Log[component_idx][(p_Ray_results->Ray_log_struct.Log_length - 1) - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_local[component_idx + log_idx * e_Full_state_size];

        }

        this->Affine_param_log[(p_Ray_results->Ray_log_struct.Log_length - 1) - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_local[e_ray_affine_param + log_idx * e_Full_state_size];
        this->Geodesic_step_log[(p_Ray_results->Ray_log_struct.Log_length - 1) - log_idx] = p_Ray_results->Ray_log_struct.Ray_path_log_local[e_step + log_idx * e_Full_state_size];
    }

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        this->p_Ray_spline_instance[idx] = gsl_spline_alloc(gsl_interp_cspline, p_Ray_results->Ray_log_struct.Log_length);
        gsl_spline_init(this->p_Ray_spline_instance[idx], this->Affine_param_log, Ray_Log[idx], p_Ray_results->Ray_log_struct.Log_length);
        
        delete[] Ray_Log[idx];

    }

    this->p_Spline_accelerator = gsl_interp_accel_alloc();

}

Emission_Integrator_class::~Emission_Integrator_class() {

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        gsl_spline_free(this->p_Ray_spline_instance[idx]);

    }

    delete[] this->Affine_param_log;
    delete[] this->Geodesic_step_log;

}

void Emission_Integrator_class::Propagate_Stokes_Vector(Geodesic_Integrator_enums e_Active_integrator, const double Start_Affine_Param, const double End_Affine_Param) {

    if (RK78_Fehlberg != e_Active_integrator && RK78_DP != e_Active_integrator && RK54 != e_Active_integrator) {

        throw std::runtime_error("Wrong active integrator in Emission_Integrator_class::Run_Explicit_Runge_Kutta()!");

    }

    memset(this->Intermediate_RHS_log, 0, sizeof(double) * RK78_size * e_Stokes_param_num);

    int RK_size = RK78_size;
    auto Stage_coeff = this->RK78_DP_Coeff_deriv;
    auto Main_solution_coeff = this->RK78_DP_Coeff_sol_main;
    auto Embedded_solution_coeff = this->RK78_DP_Coeff_sol_embeded;
    auto Affine_param_coeff = this->RK78_DP_Coeff_affine_param;

    if (RK78_Fehlberg == e_Active_integrator) {

        Stage_coeff = this->RK78_Fhelberg_Coeff_deriv;
        Main_solution_coeff = this->RK78_Fhelberg_Coeff_sol_main;
        Embedded_solution_coeff = this->RK78_Fhelberg_Coeff_sol_embeded;
        Affine_param_coeff = this->RK78_Fhelberg_Coeff_affine_param;

    }
    else if (RK54 == e_Active_integrator) {

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

                    Temp_Stokes_Vector[state_idx] += CGS_controller_step * Stage_coeff[RK_stage][derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];

                }
            }

            Temp_affine_param = this->Current_Affine_Param + Affine_param_coeff[RK_stage] * Geometric_controller_step;

            /* ------------------------------------------------------------ Get the radiative transfer RHS ------------------------------------------------------------ */

            Transfer_functions_type Total_Transfer_Functions{};

            for (int emission_medium = Disk; emission_medium <= Hotspot; emission_medium++) {

                Transfer_functions_type Temp_Transfer_functions{};

                this->p_Emission_Model->get_radiative_transfer_functions(this->get_ray_State_Vector(Temp_affine_param),
                                                                         this->p_Sim_Context,
                                                                         static_cast<Emission_medium_enums>(emission_medium),
                                                                         &Temp_Transfer_functions);

                add_vectors(Temp_Transfer_functions.Absorbtion_functions, Total_Transfer_Functions.Absorbtion_functions, e_Stokes_param_num, Total_Transfer_Functions.Absorbtion_functions);
                add_vectors(Temp_Transfer_functions.Emission_functions, Total_Transfer_Functions.Emission_functions, e_Stokes_param_num, Total_Transfer_Functions.Emission_functions);
                add_vectors(Temp_Transfer_functions.Faradey_functions, Total_Transfer_Functions.Faradey_functions, e_Stokes_param_num, Total_Transfer_Functions.Faradey_functions);

            }

            Radiative_transfer_RHS(Total_Transfer_Functions.Emission_functions,
                                   Total_Transfer_Functions.Absorbtion_functions,
                                   Total_Transfer_Functions.Faradey_functions,
                                   Temp_Stokes_Vector,
                                   this->Intermediate_RHS_log + RK_stage * e_Stokes_param_num);

            /* ------------------------------------------------------------------------------------------------------------------------------------------------------- */

        }

        for (int state_idx = 0; state_idx < e_Stokes_param_num; state_idx++) {

            New_Stokes_vector_main[state_idx] = this->Current_Stokes_Vector[state_idx];
            New_Stokes_vector_embeded[state_idx] = this->Current_Stokes_Vector[state_idx];

            for (int derivative_idx = 0; derivative_idx < RK_size; derivative_idx++) {

                New_Stokes_vector_main[state_idx] += CGS_controller_step * Main_solution_coeff[derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];
                New_Stokes_vector_embeded[state_idx] += CGS_controller_step * Embedded_solution_coeff[derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Stokes_param_num];

            }

            state_error[state_idx] = New_Stokes_vector_main[state_idx] - New_Stokes_vector_embeded[state_idx];

        }

        //if (ERROR == this->Run_NaN_checker(New_Stokes_vector_main, New_Stokes_vector_embeded)) { return; }

        this->p_Step_controller->update_state_errors(New_Stokes_vector_main, state_error, this->e_Active_integrator, e_Stokes_param_num);

        if (this->p_Step_controller->current_err < 1.0 || !this->p_Step_controller->Parameters.Use_adaptive_step) {

            this->continue_integration = true;

        }
        else {

            this->continue_integration = false;
            this->N_steps_rejected++;

        }

        if (this->continue_integration) {

            memcpy(this->Current_Stokes_Vector, New_Stokes_vector_main, e_Stokes_param_num * sizeof(double));
            this->Current_Affine_Param += Geometric_controller_step;
            //this->Update_emission_log(New_Stokes_vector_main);
            //this->Update_debug_log();

            this->N_steps_rejected = 0;
            this->NaN_checker_count = 0;

        }

        this->p_Step_controller->update_step(New_Stokes_vector_main, this->e_Active_integrator);

    }

}

void Emission_Integrator_class::Radiative_transfer_RHS(const double* const Emission_Functions,
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

const double* const Emission_Integrator_class::get_ray_State_Vector(const double Affine_param) {

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        this->Current_State_Vector[idx] = gsl_spline_eval(this->p_Ray_spline_instance[idx], Affine_param, this->p_Spline_accelerator);

    }

    return this->Current_State_Vector;

}

const double* const Emission_Integrator_class::get_current_Stokes_Vector() const {

    return this->Current_Stokes_Vector;

}