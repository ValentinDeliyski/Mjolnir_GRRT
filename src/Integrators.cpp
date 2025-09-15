#include "Integrators.h"

Step_controller_class::Step_controller_class(const Integrator_parameters_type Integrator_parameters) {

    this->Parameters = Integrator_parameters;

    this->step = this->Parameters.Init_stepzie;
    this->previous_step = this->Parameters.Init_stepzie;

    this->current_err  = this->Parameters.Safety_2;
    this->prev_err     = this->Parameters.Safety_2;
    this->sec_prev_err = this->Parameters.Safety_2;

}

void Step_controller_class::update_state_errors(const double* State_Vector, const double* State_Error_Vector, Geodesic_Integrator_enums e_Active_integrator) {

    this->sec_prev_err = this->prev_err;
    this->prev_err = this->current_err;

    double Abs_tol{}, Rel_tol{};

    switch (e_Active_integrator) {

    default:

        Abs_tol = this->Parameters.RK_78_abs_accuracy; 
        Rel_tol = this->Parameters.RK_78_rel_accuracy;

        break;

    case ESDIRK54:

        Abs_tol = this->Parameters.ESDIRK54_abs_accuracy;
        Rel_tol = this->Parameters.ESDIRK54_rel_accuracy;

    }

    double Error_scale[e_Dynamic_state_size]{};

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        Error_scale[idx] = Abs_tol + std::fabs(State_Vector[idx]) * Rel_tol;

    }

    double Total_State_Error{};

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        Total_State_Error += (State_Error_Vector[idx] / Error_scale[idx]) * (State_Error_Vector[idx] / Error_scale[idx]);

    }

    Total_State_Error /= (e_Dynamic_state_size - 1);

    this->current_err = std::sqrt(Total_State_Error) + this->Parameters.Safety_2;

}

void Step_controller_class::update_step(const double* const State_Vector, Geodesic_Integrator_enums e_Active_integrator) {

    this->previous_step = this->step;

    double PID_gain_I{}, PID_gain_P{}, PID_gain_D{}, Gustafsson_k1{}, Gustafsson_k2{};

    switch (e_Active_integrator) {

    default:

        PID_gain_P = this->Parameters.RK78_PID_gain_P;
        PID_gain_I = this->Parameters.RK78_PID_gain_I;
        PID_gain_D = this->Parameters.RK78_PID_gain_D;

        Gustafsson_k1 = this->Parameters.RK78_Gustafsson_k1;
        Gustafsson_k2 = this->Parameters.RK78_Gustafsson_k2;

        break;

    case ESDIRK54:

        PID_gain_P = this->Parameters.ESDIRK54_PID_gain_P;
        PID_gain_I = this->Parameters.ESDIRK54_PID_gain_I;
        PID_gain_D = this->Parameters.ESDIRK54_PID_gain_D;

        Gustafsson_k1 = this->Parameters.ESDIRK54_Gustafsson_k1;
        Gustafsson_k2 = this->Parameters.ESDIRK54_Gustafsson_k2;

    }

    if (!this->Parameters.Use_adaptive_step) { return; }

    double Rel_step_increase{};

    switch (this->Parameters.Controller_type) {

    case PID:

        Rel_step_increase = this->Parameters.Safety_1 * pow(this->current_err, PID_gain_I) *
                                                        pow(this->prev_err, PID_gain_P) *
                                                        pow(this->sec_prev_err, PID_gain_D);

        break;

    default:

        Rel_step_increase = this->Parameters.Safety_1 * pow(this->current_err, Gustafsson_k1) *
                                                        pow(this->current_err / this->prev_err, Gustafsson_k2);

        break;

    }

    Rel_step_increase = std::min(this->Parameters.Max_rel_step_increase, std::max(this->Parameters.Min_rel_step_increase, Rel_step_increase));

    this->step *= Rel_step_increase;

    if (this->step > this->Parameters.Max_stepsize) { this->step = this->Parameters.Max_stepsize; };

}

static int implicit_method_system_wrapper_f(const gsl_vector* gsl_trial_State_Vector, void* Params, gsl_vector* gsl_System_to_solve) {

    RHS_wrapper_struct* RHS_wrapper_params = (RHS_wrapper_struct*)Params;

    return RHS_wrapper_params->Integrator->get_implicit_method_system(gsl_trial_State_Vector, RHS_wrapper_params->p_Iteration_number, gsl_System_to_solve);

}

Integrator_class::Integrator_class(const Simulation_Context_type* const p_Sim_Context, Results_type* p_Ray_results) {

    this->continue_integration = true;
    this->integration_complete = false;
    this->Normal_termination_condition = false;
    this->Max_affine_param_reached = false;
    this->Max_integration_count_reached = false;
    this->Step_too_small = false;

    /* ---------- Placeholder ------------ */
    this->e_Active_integrator = RK78_adaptive_step;

    this->p_Init_conditions = p_Sim_Context->p_Init_Conditions;
    this->p_Spacetime = p_Sim_Context->p_Spacetime;

    this->p_Step_controller = new Step_controller_class(this->p_Init_conditions->Integrator_params);

    this->p_Ray_log_struct = &p_Ray_results->Ray_log_struct;

    /* ------------------------------ Construct the initial state vector ------------------------------ */ 

    this->p_Ray_log_struct->Ray_path_log[e_t] = p_Sim_Context->p_Init_Conditions->Observer_params.init_time;
    this->p_Ray_log_struct->Ray_path_log[e_r] = p_Sim_Context->p_Init_Conditions->Observer_params.distance;
    this->p_Ray_log_struct->Ray_path_log[e_theta] = p_Sim_Context->p_Init_Conditions->Observer_params.inclination;
    this->p_Ray_log_struct->Ray_path_log[e_phi] = p_Sim_Context->p_Init_Conditions->Observer_params.azimuth;
    this->p_Ray_log_struct->Ray_path_log[e_p_phi] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_phi];
    this->p_Ray_log_struct->Ray_path_log[e_p_theta] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_theta];
    this->p_Ray_log_struct->Ray_path_log[e_p_r] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_r];
    this->p_Ray_log_struct->Ray_path_log[e_p_t] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_t];
    this->p_Ray_log_struct->Ray_path_log[e_step] = p_Sim_Context->p_Init_Conditions->Integrator_params.Init_stepzie;
    this->p_Ray_log_struct->Ray_path_log[e_affine_param] = 0;

    this->RK_Integrator_debug_log.N_steps_rejected = p_Ray_results->RK_integrator_debug_log.N_steps_rejected;
    this->RK_Integrator_debug_log.State_error_history = p_Ray_results->RK_integrator_debug_log.State_error_history;
    this->N_steps_rejected = 0;
    this->NaN_checker_count = 0;

    /* --------------- Allocate space for the root finder, and its trial gsl_vector --------------- */

    this->Root_finder = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, e_Dynamic_state_size);
    this->gsl_trial_State_Vector = gsl_vector_alloc(e_Dynamic_state_size);

    this->RHS_Wrapper_params = { this, nullptr };
    this->Function_to_solve = { &implicit_method_system_wrapper_f,
                                e_Dynamic_state_size,
                               &RHS_Wrapper_params };
}

Integrator_class::~Integrator_class() {

    free(this->p_Step_controller);
    gsl_multiroot_fsolver_free(this->Root_finder);
    gsl_vector_free(this->gsl_trial_State_Vector);

}

int Integrator_class::get_implicit_method_system(const gsl_vector* gsl_trial_State_Vector, void* p_Iteration_number, gsl_vector* System_to_solve) {

    int Iteration_number = *(int*)p_Iteration_number;

    /* ------------ Convert the gsl_vector to a normal double* so I can pass it to the EOM functions ------------  */
    double trial_State_Vector[e_Dynamic_state_size]{};

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        trial_State_Vector[idx] = gsl_vector_get(gsl_trial_State_Vector, idx);

    }

    /* ----- Updates the Intermediate_RHS_log with the value at the current intermediate state estimate, which I call a "trial" state ----- */
    this->p_Spacetime->get_EOM(trial_State_Vector, &this->Intermediate_RHS_log[Iteration_number * e_Dynamic_state_size]);

    double Func_to_minimize[e_Dynamic_state_size]{};

    for (int state_idx = 0; state_idx < e_Dynamic_state_size; state_idx++) {

        Func_to_minimize[state_idx] = trial_State_Vector[state_idx] - this->get_current_State_Vector()[state_idx];

        for (int derivative_idx = 0; derivative_idx <= Iteration_number; derivative_idx++) {

            Func_to_minimize[state_idx] -= -this->p_Step_controller->step * ESDIRK54_Coeff_deriv[Iteration_number][derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Dynamic_state_size];
            
        }

        gsl_vector_set(System_to_solve, state_idx, Func_to_minimize[state_idx]);

    }

    return GSL_SUCCESS;
}

void Integrator_class::Run_ESDIRK54() {

    double state_error[e_Dynamic_state_size]{};

    double New_State_vector_O6[e_Dynamic_state_size]{};
    double New_State_vector_O5[e_Dynamic_state_size]{};

    /* ---------------------- Set the RHS log to zero (just in case) ---------------------- */
    memset(this->Intermediate_RHS_log, 0, sizeof(double) * RK78_size * e_Dynamic_state_size);

    /* - Compute the RHS as the current position (no point in doing this inside the loop) - */
    this->p_Spacetime->get_EOM(this->get_current_State_Vector(), this->Intermediate_RHS_log);

    /* ----------------- Convert the current state vector to a gsl_vector ----------------- */
    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        gsl_vector_set(this->gsl_trial_State_Vector, idx, this->get_current_State_Vector()[idx]);

    }

    for (int Integrator_stage = 1; Integrator_stage < ESDIRK54_size; Integrator_stage++) {

        this->RHS_Wrapper_params.p_Iteration_number = &Integrator_stage;

        gsl_multiroot_fsolver_set(this->Root_finder, &this->Function_to_solve, this->gsl_trial_State_Vector);

        int Root_finder_status = GSL_CONTINUE;

        while (GSL_CONTINUE == Root_finder_status) {

            Root_finder_status = gsl_multiroot_fsolver_iterate(this->Root_finder);

            if (GSL_SUCCESS != Root_finder_status) {

                /* TODO ERROR HANDLING */

                this->p_Step_controller->step /= 2;
                this->continue_integration = false;

                return;

            }

            Root_finder_status = gsl_multiroot_test_residual(this->Root_finder->f, 1e-10);

        }

        gsl_vector_memcpy(this->gsl_trial_State_Vector, this->Root_finder->x);

    }

    for (int state_idx = 0; state_idx < e_Dynamic_state_size; state_idx++) {

        New_State_vector_O6[state_idx] = this->get_current_State_Vector()[state_idx];
        New_State_vector_O5[state_idx] = this->get_current_State_Vector()[state_idx];

        for (int derivative_idx = 0; derivative_idx < ESDIRK54_size; derivative_idx++) {

            New_State_vector_O6[state_idx] += -this->p_Step_controller->step * ESDIRK54_Coeff_sol[derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Dynamic_state_size];
            New_State_vector_O5[state_idx] += -this->p_Step_controller->step * ESDIRK54_Coeff_test_sol[derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Dynamic_state_size];

        }

        state_error[state_idx] = New_State_vector_O6[state_idx] - New_State_vector_O5[state_idx];

    }

    if (ERROR == this->Run_NaN_checker(New_State_vector_O6, New_State_vector_O5)) { return; }

    this->p_Step_controller->update_state_errors(New_State_vector_O6, state_error, this->e_Active_integrator);
    this->p_Step_controller->update_step(New_State_vector_O6, this->e_Active_integrator);

    if (this->p_Step_controller->current_err < 1.0 || !this->p_Step_controller->Parameters.Use_adaptive_step) {

        this->continue_integration = true;

    }
    else {

        this->continue_integration = false;
        this->N_steps_rejected++;
    }

    if (this->continue_integration) {

        this->Update_ray_log(New_State_vector_O6);
        this->Update_debug_log();

        this->N_steps_rejected = 0;
        this->NaN_checker_count = 0;

    }

}

bool Integrator_class::Run_NaN_checker(const double* const New_State, const double* const New_State_Embeded) {

    for (int idx = 0; idx < e_Full_state_size; idx++) {

        if (isnan(New_State[idx]) || isnan(New_State_Embeded[idx]) || isinf(New_State[idx]) || isinf(New_State_Embeded[idx])) {

            this->continue_integration = false;
            this->p_Step_controller->step /= 10.0;
            this->NaN_checker_count++;

            return ERROR;

        }

    }

    return OK;

}

void Integrator_class::Run_RK78() {

    memset(this->Intermediate_RHS_log, 0, sizeof(double) * RK78_size * e_Dynamic_state_size);

    const double* State_Vector = this->get_current_State_Vector();
 
    double state_error[e_Dynamic_state_size]{};
    double temp_State_vector[e_Dynamic_state_size]{};

    double New_State_vector_O8[e_Dynamic_state_size]{};
    double New_State_vector_O9[e_Dynamic_state_size]{};

    for (int iteration = 0; iteration < RK78_size; iteration++) {

        memcpy(temp_State_vector, State_Vector, e_Dynamic_state_size * sizeof(double));

        for (int state_idx = 0; state_idx < e_Dynamic_state_size; state_idx++) {

            for (int derivative_idx = 0; derivative_idx < iteration; derivative_idx++) {

                temp_State_vector[state_idx] += -this->p_Step_controller->step * RK78_Coeff_deriv[iteration][derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Dynamic_state_size];

            }
        }

        this->p_Spacetime->get_EOM(temp_State_vector, &this->Intermediate_RHS_log[iteration * e_Dynamic_state_size]);

    }

    for (int state_idx = 0; state_idx < e_Dynamic_state_size; state_idx++) {

        New_State_vector_O8[state_idx] = State_Vector[state_idx];
        New_State_vector_O9[state_idx] = State_Vector[state_idx];

        for (int derivative_idx = 0; derivative_idx < RK78_size; derivative_idx++) {

            New_State_vector_O8[state_idx] += -this->p_Step_controller->step * RK78_Coeff_sol[derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Dynamic_state_size];
            New_State_vector_O9[state_idx] += -this->p_Step_controller->step * RK78_Coeff_test_sol[derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Dynamic_state_size];

        }

        state_error[state_idx] = New_State_vector_O8[state_idx] - New_State_vector_O9[state_idx];

    }

    if (ERROR == this->Run_NaN_checker(New_State_vector_O8, New_State_vector_O9)) { return; }

    this->p_Step_controller->update_state_errors(New_State_vector_O8, state_error, this->e_Active_integrator);
    this->p_Step_controller->update_step(New_State_vector_O8, this->e_Active_integrator);

    if (this->p_Step_controller->current_err < 1.0 || !this->p_Step_controller->Parameters.Use_adaptive_step){

        this->continue_integration = true;

    }
    else {

        this->continue_integration = false;
        this->N_steps_rejected++;

    }

    if (this->continue_integration) {

        this->Update_ray_log(New_State_vector_O8);
        this->Update_debug_log();
        
        this->N_steps_rejected = 0;
        this->NaN_checker_count = 0;

    }

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

void Integrator_class::Update_debug_log() {

    this->RK_Integrator_debug_log.N_steps_rejected[this->p_Ray_log_struct->Log_offset] = this->N_steps_rejected;
    this->RK_Integrator_debug_log.State_error_history[this->p_Ray_log_struct->Log_offset] = this->p_Step_controller->current_err;

}

void Integrator_class::Propagate_ray() {

    switch (this->e_Active_integrator) {

    case ESDIRK54:

        this->Run_ESDIRK54();
        this->Check_integration_complete_status();

        break;

    default:

        this->Run_RK78();
        this->Check_integration_complete_status();

        if (this->Max_integration_count_reached || this->Step_too_small || this->NaN_checker_count > 5) {
            
            this->e_Active_integrator = ESDIRK54;
            this->p_Ray_log_struct->Log_offset = 0;
        
            std::cout << "Switching to ESDIRK54... \n";

        }

        break;

    }

}

void Integrator_class::Check_integration_complete_status() {

    this->Normal_termination_condition  = this->p_Spacetime->terminate_integration(this->get_current_State_Vector());
    this->Max_affine_param_reached      = std::abs(this->get_current_State_Vector()[e_affine_param]) >= this->p_Step_controller->Parameters.Max_affine_param;
    this->Max_integration_count_reached = this->p_Ray_log_struct->Log_offset >= this->p_Step_controller->Parameters.Max_integration_count;
    this->Step_too_small                = this->p_Step_controller->step < std::numeric_limits<double>::min();

    if (this->Max_affine_param_reached) { 

        switch (this->e_Active_integrator) {

        case ESDIRK54:

            std::cout << "Max affine parameter value reachedl with ESDIRK54! \n";

            break;

        default:

            std::cout << "Max affine parameter value reached with RK78! \n";

            break;

        }
    }

    if (this->Max_integration_count_reached) { 

        switch (this->e_Active_integrator) {

        case ESDIRK54:

            std::cout << "Max iterations reached with ESDIRK54! \n";

            break;

        default:

            std::cout << "Max iterations reached with RK78! \n";

            break;

        }
    }

    if (this->Step_too_small) {

        switch (this->e_Active_integrator) {

        case ESDIRK54:

            std::cout << "Step too small with ESDIRK54! \n";

            break;

        default:

            std::cout << "Step too small with RK78! \n";

            break;

        }
    }

    this->integration_complete = this->Normal_termination_condition || this->Max_affine_param_reached;

    if (ESDIRK54 == this->e_Active_integrator) {

        this->integration_complete = this->integration_complete || this->Step_too_small || this->Max_integration_count_reached;

    }

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