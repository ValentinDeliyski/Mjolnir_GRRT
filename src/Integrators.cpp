#include "Integrators.h"

Step_controller_class::Step_controller_class(const Integrator_parameters_type Integrator_parameters) {

    this->Parameters = Integrator_parameters;

    this->step = this->Parameters.Init_stepzie;
    this->previous_step = this->Parameters.Init_stepzie;

    this->current_err  = this->Parameters.Safety_2;
    this->prev_err     = this->Parameters.Safety_2;
    this->sec_prev_err = this->Parameters.Safety_2;

}

void Step_controller_class::update_state_errors(const double* State_Vector, const double* State_Error_Vector) {

    this->sec_prev_err = this->prev_err;
    this->prev_err = this->current_err;

    double Error_scale[e_Dynamic_state_size]{};

    for (int idx = 1; idx < e_Dynamic_state_size; idx++) {

        Error_scale[idx] = this->Parameters.RK_78_abs_accuracy + std::fabs(State_Vector[idx]) * this->Parameters.RK_78_rel_accuracy;

    }

    double Total_State_Error{};

    for (int idx = 1; idx < e_Dynamic_state_size; idx++) {

        Total_State_Error += (State_Error_Vector[idx] / Error_scale[idx]) * (State_Error_Vector[idx] / Error_scale[idx]);

    }

    Total_State_Error /= (e_Dynamic_state_size - 1);

    this->current_err = std::sqrt(Total_State_Error) + this->Parameters.Safety_2;

}

void Step_controller_class::update_step(const double* const State_Vector) {

    this->previous_step = this->step;

    if (!this->Parameters.Use_adaptive_step) { return; }

    double Rel_step_increase{};

    switch (this->Parameters.Controller_type) {

    case PID:

        Rel_step_increase = this->Parameters.Safety_1 * pow(this->current_err, this->Parameters.PID_gain_I) *
                                                        pow(this->prev_err, this->Parameters.PID_gain_P) *
                                                        pow(this->sec_prev_err, this->Parameters.PID_gain_D);

        break;

    default:

        Rel_step_increase = this->Parameters.Safety_1 * pow(this->current_err, this->Parameters.Gustafsson_k1) *
                                                        pow(this->current_err / this->prev_err, this->Parameters.Gustafsson_k2);

        break;

    }

    Rel_step_increase = std::min(this->Parameters.Max_rel_step_increase, std::max(this->Parameters.Min_rel_step_increase, Rel_step_increase));

    this->step *= Rel_step_increase;

    if (this->step > this->Parameters.Max_stepsize) { this->step = this->Parameters.Max_stepsize; };

}

/* ======================================================================== GSL Wrapper functions ======================================================================== */

static int implicit_method_system_wrapper_f(const gsl_vector* gsl_trial_State_Vector, void* Params, gsl_vector* gsl_System_to_solve) {

    RHS_wrapper_struct* RHS_wrapper_params = (RHS_wrapper_struct*)Params;

    return RHS_wrapper_params->Integrator->get_implicit_method_system(gsl_trial_State_Vector, RHS_wrapper_params->RHS_params, gsl_System_to_solve);

}

static int implicit_method_system_wrapper_df(const gsl_vector* gsl_trial_State_Vector, void* Params, gsl_matrix* gsl_Jacobian) {

    RHS_wrapper_struct* RHS_wrapper_params = (RHS_wrapper_struct*)Params;

    /* Convert the gsl_vector to a normal double* so I can pass it to the EOM functions */

    double trial_State_Vector[e_Dynamic_state_size]{};

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        trial_State_Vector[idx] = gsl_vector_get(gsl_trial_State_Vector, 0);

    }

    double EOM_Jacobian[e_Dynamic_state_size][e_Dynamic_state_size]{};
    RHS_wrapper_params->Integrator->p_Spacetime->get_EOM_Jacobian(trial_State_Vector, EOM_Jacobian);

    for (int idx_1 = 0; idx_1 < e_Dynamic_state_size; idx_1++) {

        for (int idx_2 = 0; idx_2 < e_Dynamic_state_size; idx_2++) {

            double value_to_set = -2 * RHS_wrapper_params->Integrator->p_Step_controller->step / 3 * EOM_Jacobian[idx_1][idx_2];

            if (idx_1 == idx_2) {

                value_to_set += 1.;

            }

            gsl_matrix_set(gsl_Jacobian, idx_1, idx_2, value_to_set);

        }
    }

    return GSL_SUCCESS;

}

static int implicit_method_system_wrapper_fdf(const gsl_vector* gsl_trial_State_Vector, void* Params, gsl_vector* gsl_System_to_solve, gsl_matrix* gsl_Jacobian) {

    RHS_wrapper_struct* RHS_wrapper_params = (RHS_wrapper_struct*)Params;

    /* Convert the gsl_vector to a normal double* so I can pass it to the EOM functions */

    double trial_State_Vector[e_Dynamic_state_size]{};

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        trial_State_Vector[idx] = gsl_vector_get(gsl_trial_State_Vector, 0);

    }

    implicit_method_system_wrapper_f(gsl_trial_State_Vector, Params, gsl_System_to_solve);
    implicit_method_system_wrapper_df(gsl_trial_State_Vector, Params, gsl_Jacobian);

    return GSL_SUCCESS;

}

/* ======================================================================================================================================================================== */

Integrator_class::Integrator_class(const Simulation_Context_type* const p_Sim_Context, Results_type* p_Ray_results) {

    this->In_stiff_region = false;
    this->Implicit_method_init_status = false;
    this->continue_integration = true;
    this->integration_complete = false;

    this->e_Active_integrator = RK78_adaptive_step;

    this->p_Init_conditions = p_Sim_Context->p_Init_Conditions;
    this->p_Spacetime = p_Sim_Context->p_Spacetime;

    this->p_Step_controller = new Step_controller_class(this->p_Init_conditions->Integrator_params);

    this->p_Ray_log_struct = &p_Ray_results->Ray_log_struct;

    /* ------------------------------ Construct the initial state vector ------------------------------ */ 
    double Init_State_Vector[e_Full_state_size]{};

    Init_State_Vector[e_t] = p_Sim_Context->p_Init_Conditions->Observer_params.init_time;
    Init_State_Vector[e_r] = p_Sim_Context->p_Init_Conditions->Observer_params.distance;
    Init_State_Vector[e_theta] = p_Sim_Context->p_Init_Conditions->Observer_params.inclination;
    Init_State_Vector[e_phi] = p_Sim_Context->p_Init_Conditions->Observer_params.azimuth;
    Init_State_Vector[e_p_phi] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_phi];
    Init_State_Vector[e_p_theta] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_theta];
    Init_State_Vector[e_p_r] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_r];
    Init_State_Vector[e_p_t] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_t];
    Init_State_Vector[e_step] = p_Sim_Context->p_Init_Conditions->Integrator_params.Init_stepzie;
    Init_State_Vector[e_affine_param] = 0;

    memcpy(this->p_Ray_log_struct->Ray_path_log, Init_State_Vector, e_Full_state_size * sizeof(double));

    /* --------------- Allocate space for the root finder, and its trial gsl_vector --------------- */

    this->Root_finder = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, e_Dynamic_state_size);
    this->gsl_trial_State_Vector = gsl_vector_alloc(e_Dynamic_state_size);

    this->RHS_Wrapper_params = { this, nullptr };
    this->Function_to_solve = { &implicit_method_system_wrapper_f,
                                e_Dynamic_state_size, 
                               &RHS_Wrapper_params };

}

int Integrator_class::get_implicit_method_system(const gsl_vector* gsl_trial_State_Vector, void* Params, gsl_vector* System_to_solve) {

    /* Convert the gsl_vector to a normal double* so I can pass it to the EOM functions */

    double trial_State_Vector[e_Dynamic_state_size]{};

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        trial_State_Vector[idx] = gsl_vector_get(gsl_trial_State_Vector, idx);

    }

    /* ---------------------------- References for the sake of readability ---------------------------- */
    int& Current_state_idx = this->p_Ray_log_struct->Log_offset;
    double* State_Vector = &this->p_Ray_log_struct->Ray_path_log[Current_state_idx * e_Full_state_size];

    double* Old_State = &this->p_Ray_log_struct->Ray_path_log[(Current_state_idx - 1) * e_Full_state_size];
    /* ------------------------------------------------------------------------------------------------ */

    double RHS[e_Dynamic_state_size]{};

    switch (this->e_Active_integrator) {

    case BDF:

        this->p_Spacetime->get_EOM(trial_State_Vector, RHS);

        for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

            gsl_vector_set(System_to_solve, idx, trial_State_Vector[idx] - 4. / 3 * State_Vector[idx] + 1. / 3 * Old_State[idx] - 2. * (-this->p_Step_controller->step) / 3. * RHS[idx]);

        }

    }

    return GSL_SUCCESS;
}

void Integrator_class::Run_RK78() {

    /* ---------------------------- References for the sake of readability ---------------------------- */
    int& Current_state_idx = this->p_Ray_log_struct->Log_offset;
    double* State_Vector = &this->p_Ray_log_struct->Ray_path_log[Current_state_idx * e_Full_state_size];
    /* ------------------------------------------------------------------------------------------------ */

    double state_error[e_Dynamic_state_size]{};

    double Derivatives[RK78_size * e_Dynamic_state_size]{};

    double temp_State_vector[e_Dynamic_state_size]{};

    double New_State_vector_O8[e_Dynamic_state_size]{};
    double New_State_vector_O9[e_Dynamic_state_size]{};

    for (int iteration = 0; iteration < RK78_size; iteration++) {

        memcpy(temp_State_vector, State_Vector, e_Dynamic_state_size * sizeof(double));

        for (int vector_indexer = 0; vector_indexer < e_Dynamic_state_size; vector_indexer++) {

            for (int derivative_indexer = 0; derivative_indexer < iteration; derivative_indexer++) {

                temp_State_vector[vector_indexer] += -this->p_Step_controller->step * RK78_Coeff_deriv[iteration][derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_Dynamic_state_size];

            }
        }

        this->p_Spacetime->get_EOM(temp_State_vector, &Derivatives[iteration * e_Dynamic_state_size]);

    }

    for (int vector_indexer = 0; vector_indexer < e_Dynamic_state_size; vector_indexer++) {

        New_State_vector_O8[vector_indexer] = State_Vector[vector_indexer];
        New_State_vector_O9[vector_indexer] = State_Vector[vector_indexer];

        for (int derivative_indexer = 0; derivative_indexer < RK78_size; derivative_indexer++) {

            New_State_vector_O8[vector_indexer] += -this->p_Step_controller->step * RK78_Coeff_sol[derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_Dynamic_state_size];
            New_State_vector_O9[vector_indexer] += -this->p_Step_controller->step * RK78_Coeff_test_sol[derivative_indexer] * Derivatives[vector_indexer + derivative_indexer * e_Dynamic_state_size];

        }

        state_error[vector_indexer] = New_State_vector_O8[vector_indexer] - New_State_vector_O9[vector_indexer];

    }

    this->p_Step_controller->update_state_errors(New_State_vector_O8, state_error);

    // The integrator might jump pass surfaces that are singular for the EOM (like the JNW singularity at 2 / gamma)
    // In this case the whole state vector becomes a NaN. I check for this and update the integration step by hand,
    // then set the continue_integration flag to "false" to force the integrator to redo the current iteration with a smaller step.
    if (isnan(this->p_Step_controller->current_err)) {

        this->continue_integration = false;
        this->p_Step_controller->step /= 10.0;

        return;

    }

    this->p_Step_controller->update_step(New_State_vector_O8);

    if (this->p_Step_controller->current_err < 1.0 || !this->p_Step_controller->Parameters.Use_adaptive_step)
    {
        this->continue_integration = true;
    }
    else
    {
        this->continue_integration = false;
    }

    if (this->continue_integration) {

        // For the JNW Naked Singularity, certain photons scatter from very close to the singularity.
        // Close enough that it requires "manual" scattering, by flipping the p_r sign.
        // Otherwise the photons never reach the turning point and the integration grinds to a halt.
        if (this->p_Init_conditions->Metric_parameters.e_Spacetime == Janis_Newman_Winicour && this->p_Init_conditions->Metric_parameters.JNW_Gamma_Parameter < 0.5) {

            if (State_Vector[e_r] - 2.0 / this->p_Init_conditions->Metric_parameters.JNW_Gamma_Parameter < this->p_Init_conditions->Metric_parameters.Min_distance_to_singular_point) {

                State_Vector[e_p_r] *= -1.0;

            }

        }

        this->Update_ray_log(New_State_vector_O8);

    }

}

void Integrator_class::Init_BDF() {

    this->p_Step_controller->Parameters.Use_adaptive_step = false;
    //this->p_Step_controller->step *= 10;

    for (int init_idx = 0; init_idx < 5; init_idx++) {

        this->Run_RK78();
        this->Check_integration_complete_status();

    }

    this->Implicit_method_init_status = true;

}

void print_state(size_t iter, gsl_multiroot_fsolver* s)
{
    printf("iter = %3u x = % .6f % .6f "
        "f(x) = % .6e % .6e\n",
        iter,
        gsl_vector_get(s->x, 0),
        gsl_vector_get(s->x, 1),
        gsl_vector_get(s->f, 0),
        gsl_vector_get(s->f, 1),
        s->state);
}

int Integrator_class::Run_BDF() {

    /* -------------------- Convert the current state vector to a gsl_vector -------------------- */

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        gsl_vector_set(this->gsl_trial_State_Vector, idx, this->get_current_State_Vector()[idx]);

    }

    gsl_multiroot_fsolver_set(this->Root_finder, &this->Function_to_solve, this->gsl_trial_State_Vector);

    int Root_finder_status = GSL_CONTINUE;
    int Root_finder_iteration = 0;

    // && Root_finder_iteration <= this->p_Init_conditions->Integrator_params.BDF_root_finder_max_iterations

    //std::cout << "Current r:" << this->get_current_State_Vector()[e_r] << "\n";

    while (GSL_CONTINUE == Root_finder_status) {

        Root_finder_status = gsl_multiroot_fsolver_iterate(this->Root_finder);

        Root_finder_iteration++;

        if (Root_finder_status)   /* check if solver is stuck */
            return Root_finder_status;

        Root_finder_status = gsl_multiroot_test_residual(this->Root_finder->f, 1e-4);


        //print_state(Root_finder_iteration, this->Root_finder);

    }

    /* --- Convert the new gsl_vector State to a normal double* so I can pass it to the Update_ray_log function --- */

    double Nnew_State_Vector[e_Dynamic_state_size]{};

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        Nnew_State_Vector[idx] = gsl_vector_get(this->Root_finder->x, idx);

    }

    this->Update_ray_log(Nnew_State_Vector);

    return GSL_SUCCESS;
}

void Integrator_class::Update_ray_log(const double* const New_State_vector) {

    this->p_Ray_log_struct->Log_offset += 1;
    int& log_offset = this->p_Ray_log_struct->Log_offset;

    memcpy(&this->p_Ray_log_struct->Ray_path_log[log_offset * e_Full_state_size], New_State_vector, e_Dynamic_state_size * sizeof(double));

    this->p_Ray_log_struct->Ray_path_log[e_step + log_offset * e_Full_state_size] = this->p_Step_controller->previous_step;
    this->p_Ray_log_struct->Ray_path_log[e_affine_param + log_offset * e_Full_state_size] = this->p_Ray_log_struct->Ray_path_log[e_affine_param + (log_offset - 1) * e_Full_state_size] - this->p_Step_controller->previous_step;
    
    // The wormhole metric works with a "global" radial coordinate, that goes negative on the other side of the throat.
    // The emission model can't work with this coordinate, so I log the normal spherical radial coordinate instead. 
    if (Wormhole == this->p_Init_conditions->Metric_parameters.e_Spacetime) {

        const double& R_throat = this->p_Init_conditions->Metric_parameters.R_throat;
        this->p_Ray_log_struct->Ray_path_log[e_r + log_offset * e_Full_state_size] = sqrt(New_State_vector[e_r] * New_State_vector[e_r] + R_throat * R_throat);

    }

}

Stability_return_type Integrator_class::Check_method_stability() {

    /* ---------------------------- References for the sake of readability ---------------------------- */
    int& Current_state_idx = this->p_Ray_log_struct->Log_offset;
    double* State_Vector = &this->p_Ray_log_struct->Ray_path_log[Current_state_idx * e_Full_state_size];
    /* ------------------------------------------------------------------------------------------------ */

    std::complex<double> Stability_param = -this->p_Step_controller->step;
    std::complex<double> Stability_polynomial{};
    Stability_return_type s_Method_stability{};

    /* ------- Init all stability criteria to true ------- */

    s_Method_stability.RK78_stability_status = true;
    s_Method_stability.BDF_stability_status = true;

    /* --------------------------------------------------- */

    if (this->p_Step_controller->step > this->p_Init_conditions->Integrator_params.Step_stability_check_threshold || State_Vector[e_r] > 5) { return s_Method_stability; }

    Stability_param *= this->p_Spacetime->get_largest_EOM_eigenvalue(State_Vector);

    /* ------------------------------------- Evaluate BDF stability ------------------------------------- */

    s_Method_stability.BDF_stability_status = true;

    /* ------------------------------------- Evaluate RK78 stability ------------------------------------ */

    for (int idx = 0; idx < RK78_size; idx++) {

        Stability_polynomial += this->RK78_Stability_pol_coeffs[idx] * std::pow(Stability_param, idx);

    }

    s_Method_stability.RK78_stability_status = std::abs(Stability_polynomial) <= 1.0;

    /* ------------------------------------------------------------------------------------------------- */

    return s_Method_stability;

}

void Integrator_class::Propagate_ray() {

    switch (this->e_Active_integrator) {

    default:

        this->Run_RK78();

        break;

    case BDF:

        if (this->Implicit_method_init_status) { 
            
          int Status = this->Run_BDF();

          //if (GSL_SUCCESS != Status) {

          //    this->Implicit_method_init_status = false;
          //    this->p_Step_controller->Parameters.Use_adaptive_step = true;
          //    this->e_Active_integrator = RK78_adaptive_step;

          //}
        
        }
        else { this->Init_BDF(); }

        break;

    }

    this->Check_integration_complete_status();

    //if (0) {

    //    Stability_return_type s_Method_Stability = this->Check_method_stability();

    //    switch (this->e_Active_integrator) {

    //    case RK78_adaptive_step:

    //        if (!s_Method_Stability.RK78_stability_status) { 
    //            
    //            this->e_Active_integrator = BDF;
    //        
    //        }
    //        break;

    //    case BDF:

    //        if (s_Method_Stability.RK78_stability_status){ 
    //            
    //            this->Implicit_method_init_status = false;
    //            this->p_Step_controller->Parameters.Use_adaptive_step = true;
    //            this->e_Active_integrator = RK78_adaptive_step; 
    //        
    //        }

    //        if (!s_Method_Stability.BDF_stability_status) { 

    //            /* Placeholder - probably just log it in the hypothetical event logger */
    //        
    //        }

    //        break;

    //    }

    //}

}

void Integrator_class::Check_integration_complete_status() {

    /* ---------------------------- References for the sake of readability ---------------------------- */
    int& Current_state_idx = this->p_Ray_log_struct->Log_offset;
    double* State_Vector = &this->p_Ray_log_struct->Ray_path_log[Current_state_idx * e_Full_state_size];
    /* ------------------------------------------------------------------------------------------------ */

    bool Normal_termination_condition = this->p_Spacetime->terminate_integration(State_Vector);
    bool Max_affine_param_reached = std::abs(State_Vector[e_affine_param]) >= this->p_Step_controller->Parameters.Max_affine_param;
    bool Max_integration_count_reached = this->p_Ray_log_struct->Log_offset >= this->p_Step_controller->Parameters.Max_integration_count;

    if (Max_integration_count_reached) { std::cout << "Max iterations reached! \n"; }
    if (Max_affine_param_reached) { std::cout << "Max affine parameter value reached! \n"; };

    this->integration_complete = Normal_termination_condition || Max_affine_param_reached || Max_integration_count_reached;

}

const double* const Integrator_class::get_current_State_Vector() const {

    return &this->p_Ray_log_struct->Ray_path_log[this->p_Ray_log_struct->Log_offset * e_Full_state_size];

}

const double* const Integrator_class::get_previous_State_Vector() const {

    if (this->p_Ray_log_struct->Log_offset > 0) {

        return &this->p_Ray_log_struct->Ray_path_log[(this->p_Ray_log_struct->Log_offset - 1) * e_Full_state_size];
    }
    else {

        return &this->p_Ray_log_struct->Ray_path_log[this->p_Ray_log_struct->Log_offset * e_Full_state_size];

    }
    
}