#include "Integrators.h"

Step_controller_class::Step_controller_class(const Step_Controller_parameters_type Controller_parameters) {

    this->Parameters = Controller_parameters;

    this->step = this->Parameters.Init_stepzie;
    this->previous_step = this->Parameters.Init_stepzie;

    this->current_err  = 1.0;
    this->prev_err     = 1.0;
    this->sec_prev_err = 1.0;

}

double Step_controller_class::get_max_step(const double r) const {

    const double Step_ratio = this->Parameters.Min_upper_stepsize / this->Parameters.Max_upper_stepsize;

    /* This thing should ideally be as close to 2 * Step_ratio / (1 - Step_ratio) as possible, while still being smaller. 
       When it is exactly equal to 2 * Step_ratio / (1 - Step_ratio), arctanh_arg becomes equal to 1. */
    const double c_coeff = 0.999 * 2 * Step_ratio / (1 - Step_ratio);
    const double arctanh_arg = 1 + c_coeff - (2 + c_coeff) * Step_ratio;

    if (std::abs(arctanh_arg) >= 1) {

        return this->Parameters.Max_upper_stepsize;

    }

    const double a_coeff = this->Parameters.Max_step_b_coeff * atanh(arctanh_arg) + this->Parameters.Dist_at_min_upper_stepsize / this->Parameters.Dist_to_Observer;

    return this->Parameters.Max_upper_stepsize * (tanh((r / this->Parameters.Dist_to_Observer - a_coeff) / this->Parameters.Max_step_b_coeff) + 1. + c_coeff) / (2. + c_coeff);

}

void Step_controller_class::update_step(Integrator_enums e_Active_integrator, const double r, const bool Is_inside_emission_medium) {

    this->previous_step = this->step;

    if (!this->Parameters.Use_adaptive_step) { return; }

    double PID_gain_I{}, PID_gain_P{}, PID_gain_D{}, Gustafsson_k1{}, Gustafsson_k2{};

    switch (e_Active_integrator) {

    default:

        PID_gain_P = this->Parameters.RK_PID_gain_P;
        PID_gain_I = this->Parameters.RK_PID_gain_I;
        PID_gain_D = this->Parameters.RK_PID_gain_D;

        Gustafsson_k1 = this->Parameters.RK_Gustafsson_k1;
        Gustafsson_k2 = this->Parameters.RK_Gustafsson_k2;

        break;

    case ESDIRK54:

        PID_gain_P = this->Parameters.ESDIRK54_PID_gain_P;
        PID_gain_I = this->Parameters.ESDIRK54_PID_gain_I;
        PID_gain_D = this->Parameters.ESDIRK54_PID_gain_D;

        Gustafsson_k1 = this->Parameters.ESDIRK54_Gustafsson_k1;
        Gustafsson_k2 = this->Parameters.ESDIRK54_Gustafsson_k2;

    }

    double Rel_step_increase{};

    switch (this->Parameters.Controller_type) {

    case PID:

        Rel_step_increase = this->Parameters.Safety_1 * pow(this->current_err, -PID_gain_I) *
                                                        pow(this->prev_err, PID_gain_P) *
                                                        pow(this->sec_prev_err, -PID_gain_D);

        break;

    default:

        Rel_step_increase = this->Parameters.Safety_1 * pow(this->current_err, -Gustafsson_k1) *
                                                        pow(this->current_err / this->prev_err, Gustafsson_k2);

        break;

    }

    const double max_stepsize = this->get_max_step(r);

    Rel_step_increase = std::min(this->Parameters.Max_rel_step_increase, std::max(this->Parameters.Min_rel_step_increase, Rel_step_increase));

    this->step *= Rel_step_increase;

    if (this->step > max_stepsize) { this->step = max_stepsize; };

    if (Is_inside_emission_medium) { this->step = std::min(this->step, this->Parameters.Max_step_inisde_emission_medium); }

}

static int implicit_method_system_wrapper_f(const gsl_vector* State_Vector, void* Params, gsl_vector* gsl_System_to_solve) {

    RHS_wrapper_struct* RHS_wrapper_params = (RHS_wrapper_struct*)Params;

    return RHS_wrapper_params->Integrator->get_implicit_method_system(State_Vector, RHS_wrapper_params->p_Iteration_number, gsl_System_to_solve);

}

Geodesic_Integrator_class::Geodesic_Integrator_class(const Simulation_Context_type* const p_Sim_Context, Results_type* p_Ray_results) {

    this->e_Active_integrator = p_Sim_Context->p_Init_Conditions->Integrator_params.e_Default_geodesic_integrator;

    /* ---------------------------------------- Init the internal flags ---------------------------------------- */
    this->propagate_optical_depth = p_Sim_Context->p_Init_Conditions->Integrator_params.Propagate_optical_depth;

    this->continue_integration = true;
    /* This currently exists because of the JNW naked singularity. 
       TODO: figure out a more elegant solution */
    this->Force_scatter = true;

    /* --- Integration termination flags --- */
    this->Step_too_small = false;
    this->integration_complete = false;
    this->Max_affine_param_reached = false;
    this->Normal_termination_condition = false;
    this->Max_integration_count_reached = false;

    /* --------- Set the internal pointers to relevant classes / structs that the integrator uses --------- */
    this->p_Init_conditions = p_Sim_Context->p_Init_Conditions;
    this->p_Emission_Model = p_Sim_Context->p_Emission_Model;
    this->p_Ray_log_struct = p_Ray_results->Ray_log_struct;
    this->p_Spacetime = p_Sim_Context->p_Spacetime;

    this->Current_Optical_Depth = &p_Ray_results->Optical_Depth;
    this->Current_Faraday_Q_Depth = &p_Ray_results->Faraday_Q_Depth;
    this->Current_Faraday_V_Depth = &p_Ray_results->Faraday_V_Depth;

    /* -------- This thing is a unique pointer so I dont have to delete it manually when the integrator goes out of scope -------- */
    this->p_Step_controller = std::make_unique<Step_controller_class>(this->p_Init_conditions->Integrator_params.Geodesic_Step_Controller_Params);

    this->Max_affine_param = this->p_Init_conditions->Integrator_params.Max_affine_param;
    this->Max_integration_count = this->p_Init_conditions->Integrator_params.Max_integration_count;

    /* ----------------------------------------- Construct the initial state vector in global coords ----------------------------------------- */ 
    this->p_Ray_log_struct->Ray_path_log_global[e_t] = p_Sim_Context->p_Init_Conditions->Observer_params.init_time;
    this->p_Ray_log_struct->Ray_path_log_global[e_r] = p_Sim_Context->p_Init_Conditions->Observer_params.distance;
    this->p_Ray_log_struct->Ray_path_log_global[e_theta] = p_Sim_Context->p_Init_Conditions->Observer_params.inclination;
    this->p_Ray_log_struct->Ray_path_log_global[e_phi] = p_Sim_Context->p_Init_Conditions->Observer_params.azimuth;
    this->p_Ray_log_struct->Ray_path_log_global[e_p_phi] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_phi];
    this->p_Ray_log_struct->Ray_path_log_global[e_p_theta] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_theta];
    this->p_Ray_log_struct->Ray_path_log_global[e_p_r] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_r];
    this->p_Ray_log_struct->Ray_path_log_global[e_p_t] = p_Sim_Context->p_Init_Conditions->Init_Momentum[e_t];
    this->p_Ray_log_struct->Ray_path_log_global[e_step] = p_Sim_Context->p_Init_Conditions->Integrator_params.Geodesic_Step_Controller_Params.Init_stepzie;
    this->p_Ray_log_struct->Ray_path_log_global[e_ray_affine_param] = 0;

    /* ----------------------------------------- Construct the initial state vector in local coords ----------------------------------------- */
    this->p_Spacetime->Convert_global_to_local_coords(this->p_Ray_log_struct->Ray_path_log_global.get(), 
                                                      this->p_Ray_log_struct->Ray_path_log_global.get(), 
                                                      this->p_Ray_log_struct->Ray_path_log_local.get(),
                                                      e_Full_State_Vector);

    /* -------------------------------------- Set the initial internal dynamic state to the global one -------------------------------------- */
    memcpy(this->Current_Dynamic_state, this->p_Ray_log_struct->Ray_path_log_global.get(), e_Dynamic_state_size * sizeof(double));

    /* -------------------------- Set the internal debug tracker pointers to point to the external results struct --------------------------- */
    this->RK_Integrator_debug_log = &p_Ray_results->RK_integrator_debug_log;

    /* --- Init the per-step debug counters --- */
    this->N_steps_rejected = 0;
    this->NaN_checker_count = 0;

    /* --------------- Allocate space for the root finder, and its trial gsl_vector --------------- */
    this->Root_finder = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, e_Dynamic_state_size);
    this->gsl_trial_State_Vector = gsl_vector_alloc(e_Dynamic_state_size);

    this->RHS_Wrapper_params = { this, nullptr };
    this->Function_to_solve = { &implicit_method_system_wrapper_f,
                                e_Dynamic_state_size,
                                &this->RHS_Wrapper_params };
}

Geodesic_Integrator_class::~Geodesic_Integrator_class() {

    gsl_multiroot_fsolver_free(this->Root_finder);
    gsl_vector_free(this->gsl_trial_State_Vector);

}

int Geodesic_Integrator_class::get_implicit_method_system(const gsl_vector* gsl_State_Vector, void* p_Iteration_number, gsl_vector* System_to_solve) {

    int Iteration_number = *(int*)p_Iteration_number;

    /* ------------ Convert the gsl_vector to a normal double* so I can pass it to the EOM functions ------------  */
    double State_Vector[e_Dynamic_state_size]{};

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        State_Vector[idx] = gsl_vector_get(gsl_State_Vector, idx);

    }

    /* ----- Updates the Intermediate_RHS_log with the value at the current intermediate state estimate, which I call a "trial" state ----- */
    this->p_Spacetime->get_EOM(State_Vector, &this->Intermediate_RHS_log[Iteration_number * e_Dynamic_state_size]);

    double Func_to_minimize[e_Dynamic_state_size]{};

    for (int state_idx = 0; state_idx < e_Dynamic_state_size; state_idx++) {

        Func_to_minimize[state_idx] = State_Vector[state_idx] - this->Current_Dynamic_state[state_idx];

        for (int derivative_idx = 0; derivative_idx <= Iteration_number; derivative_idx++) {

            Func_to_minimize[state_idx] -= -this->p_Step_controller->step * ESDIRK54_Coeff_deriv[Iteration_number][derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Dynamic_state_size];
            
        }

        gsl_vector_set(System_to_solve, state_idx, Func_to_minimize[state_idx]);

    }

    return GSL_SUCCESS;
}

void Geodesic_Integrator_class::Run_ESDIRK54() {

    double state_error[e_Dynamic_state_size]{};

    double New_State_vector_main[e_Dynamic_state_size]{};
    double New_State_vector_embeded[e_Dynamic_state_size]{};

    /* ---------------------- Set the RHS log to zero (just in case) ---------------------- */
    memset(this->Intermediate_RHS_log, 0, sizeof(double) * RK78_size * e_Dynamic_state_size);

    /* - Compute the RHS as the current position (no point in doing this inside the loop) - */
    this->p_Spacetime->get_EOM(this->Current_Dynamic_state, this->Intermediate_RHS_log);

    /* ----------------- Convert the current state vector to a gsl_vector ----------------- */
    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        gsl_vector_set(this->gsl_trial_State_Vector, idx, this->Current_Dynamic_state[idx]);

    }

    for (int Integrator_stage = 1; Integrator_stage < ESDIRK54_size; Integrator_stage++) {

        this->RHS_Wrapper_params.p_Iteration_number = &Integrator_stage;

        gsl_multiroot_fsolver_set(this->Root_finder, &this->Function_to_solve, this->gsl_trial_State_Vector);

        int Root_finder_status = GSL_CONTINUE;

        while (GSL_CONTINUE == Root_finder_status) {

            Root_finder_status = gsl_multiroot_fsolver_iterate(this->Root_finder);

            if (GSL_SUCCESS != Root_finder_status) {

                this->p_Step_controller->step /= 2;
                this->continue_integration = false;

                return;

            }

            Root_finder_status = gsl_multiroot_test_residual(this->Root_finder->f, 1e-10);

        }

        gsl_vector_memcpy(this->gsl_trial_State_Vector, this->Root_finder->x);

    }

    for (int state_idx = 0; state_idx < e_Dynamic_state_size; state_idx++) {

        New_State_vector_main[state_idx] = this->Current_Dynamic_state[state_idx];
        New_State_vector_embeded[state_idx] = this->Current_Dynamic_state[state_idx];

        for (int derivative_idx = 0; derivative_idx < ESDIRK54_size; derivative_idx++) {

            New_State_vector_main[state_idx] += -this->p_Step_controller->step * ESDIRK54_Coeff_sol_main[derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Dynamic_state_size];
            New_State_vector_embeded[state_idx] += -this->p_Step_controller->step * ESDIRK54_Coeff_sol_embeded[derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Dynamic_state_size];

        }

        state_error[state_idx] = New_State_vector_main[state_idx] - New_State_vector_embeded[state_idx];

    }

    if (ERROR == this->Run_NaN_checker(New_State_vector_main, New_State_vector_embeded)) { return; }

    this->p_Step_controller->update_state_errors(New_State_vector_main, state_error, this->e_Active_integrator, e_Dynamic_state_size);

    Emission_medium_state_type Disk_State{};
    bool Is_inside_emission_medium = false;

    if (this->p_Step_controller->current_err < 1.0 or !this->p_Step_controller->Parameters.Use_adaptive_step) {

        if (this->propagate_optical_depth) {

            Is_inside_emission_medium = this->p_Emission_Model->p_Disk_Model->is_inside_disk(New_State_vector_main, &Disk_State);

        }

        if (Is_inside_emission_medium and this->p_Step_controller->step > this->p_Step_controller->Parameters.Max_step_inisde_emission_medium) {

            this->continue_integration = false;

        }
        else {

            this->continue_integration = true;

        }

    }
    else {

        this->continue_integration = false;
        this->N_steps_rejected++;

    }

    this->p_Step_controller->update_step(this->e_Active_integrator, New_State_vector_main[e_r], Is_inside_emission_medium);

    if (this->continue_integration) {

        memcpy(this->Current_Dynamic_state, New_State_vector_main, e_Dynamic_state_size * sizeof(double));

        this->Update_ray_log(New_State_vector_main);
        this->Update_optical_depth();
        this->Update_debug_log();

        this->N_steps_rejected = 0;
        this->NaN_checker_count = 0;


    }

}

void Geodesic_Integrator_class::Update_optical_depth() {

    if (this->propagate_optical_depth and this->p_Ray_log_struct->Log_offset > 0) {

        Transfer_functions_type Total_Transfer_Functions{};

        for (int emission_medium = Disk; emission_medium < e_Emission_medium_number; emission_medium++) {

            Transfer_functions_type Temp_Transfer_functions{};

            this->p_Emission_Model->get_radiative_transfer_functions(this->get_previous_State_Vector_local(),
                                                                     static_cast<Emission_medium_enums>(emission_medium),
                                                                     &Temp_Transfer_functions);

            add_vectors(Temp_Transfer_functions.Faraday_functions, Total_Transfer_Functions.Faraday_functions, e_Stokes_param_num, Total_Transfer_Functions.Faraday_functions);

        }

        *this->Current_Optical_Depth += Total_Transfer_Functions.Absorbtion_functions[I] * this->p_Step_controller->previous_step * this->p_Init_conditions->central_object_mass * MASS_TO_CM;
        *this->Current_Faraday_Q_Depth += Total_Transfer_Functions.Faraday_functions[Q] * this->p_Step_controller->previous_step * this->p_Init_conditions->central_object_mass * MASS_TO_CM;
        *this->Current_Faraday_V_Depth += Total_Transfer_Functions.Faraday_functions[V] * this->p_Step_controller->previous_step * this->p_Init_conditions->central_object_mass * MASS_TO_CM;

    }
}

Return_Values Geodesic_Integrator_class::Run_NaN_checker(const double* const New_State, const double* const New_State_Embeded) {

    for (int idx = 0; idx < e_Full_state_size; idx++) {

        if (isnan(New_State[idx]) or isnan(New_State_Embeded[idx]) or isinf(New_State[idx]) or isinf(New_State_Embeded[idx])) {

            this->continue_integration = false;
            this->p_Step_controller->step /= 10.0;
            this->NaN_checker_count++;

            return ERROR;

        }

    }

    return OK;

}

void Geodesic_Integrator_class::Run_Explicit_Runge_Kutta() {

    if (RK78_Fehlberg != this->e_Active_integrator and RK78_DP != this->e_Active_integrator and RK54 != this->e_Active_integrator) {

        throw std::runtime_error("Wrong active integrator in Geodesic_Geodesic_Integrator_class::Run_Explicit_Runge_Kutta()!");

    }

    memset(this->Intermediate_RHS_log, 0, sizeof(double) * RK78_size * e_Dynamic_state_size);

    int RK_size = RK78_size;
    auto Stage_coeff = this->RK78_DP_Coeff_deriv;
    auto Main_solution_coeff = this->RK78_DP_Coeff_sol_main;
    auto Embedded_solution_coeff = this->RK78_DP_Coeff_sol_embeded;

    if (RK78_Fehlberg == this->e_Active_integrator) {

        Stage_coeff = this->RK78_Fhelberg_Coeff_deriv;
        Main_solution_coeff = this->RK78_Fhelberg_Coeff_sol_main;
        Embedded_solution_coeff = this->RK78_Fhelberg_Coeff_sol_embeded;

    }
    else if (RK54 == this->e_Active_integrator) {

        Stage_coeff = this->RK54_Coeff_deriv;
        Main_solution_coeff = this->RK54_Coeff_sol_main;
        Embedded_solution_coeff = this->RK54_Coeff_test_embeded;

        RK_size = this->RK54_size;

    }

    double state_error[e_Dynamic_state_size]{};
    double temp_State_vector[e_Dynamic_state_size]{};

    double New_State_vector_main[e_Dynamic_state_size]{};
    double New_State_vector_embeded[e_Dynamic_state_size]{};

    for (int iteration = 0; iteration < RK_size; iteration++) {

        memcpy(temp_State_vector, this->Current_Dynamic_state, e_Dynamic_state_size * sizeof(double));

        for (int state_idx = 0; state_idx < e_Dynamic_state_size; state_idx++) {

            for (int derivative_idx = 0; derivative_idx < iteration; derivative_idx++) {

                temp_State_vector[state_idx] += -this->p_Step_controller->step * Stage_coeff[iteration][derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Dynamic_state_size];

            }
        }

        this->p_Spacetime->get_EOM(temp_State_vector, &this->Intermediate_RHS_log[iteration * e_Dynamic_state_size]);

    }

    for (int state_idx = 0; state_idx < e_Dynamic_state_size; state_idx++) {

        New_State_vector_main[state_idx] = this->Current_Dynamic_state[state_idx];
        New_State_vector_embeded[state_idx] = this->Current_Dynamic_state[state_idx];

        for (int derivative_idx = 0; derivative_idx < RK_size; derivative_idx++) {

            New_State_vector_main[state_idx] += -this->p_Step_controller->step * Main_solution_coeff[derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Dynamic_state_size];
            New_State_vector_embeded[state_idx] += -this->p_Step_controller->step * Embedded_solution_coeff[derivative_idx] * this->Intermediate_RHS_log[state_idx + derivative_idx * e_Dynamic_state_size];

        }

        state_error[state_idx] = New_State_vector_main[state_idx] - New_State_vector_embeded[state_idx];

    }

    if (ERROR == this->Run_NaN_checker(New_State_vector_main, New_State_vector_embeded)) { return; }

    this->p_Step_controller->update_state_errors(New_State_vector_main, state_error, this->e_Active_integrator, e_Dynamic_state_size);

    Emission_medium_state_type Disk_State{};
    bool Is_inside_emission_medium = false;

    if (this->p_Step_controller->current_err < 1.0 or !this->p_Step_controller->Parameters.Use_adaptive_step){

        if (this->propagate_optical_depth) {

            Is_inside_emission_medium = this->p_Emission_Model->p_Disk_Model->is_inside_disk(New_State_vector_main, &Disk_State) or
                                        this->p_Emission_Model->p_Hotspot_Model->is_inside_hotspot(New_State_vector_main, &Disk_State);

        }

        if (Is_inside_emission_medium and this->p_Step_controller->step > this->p_Step_controller->Parameters.Max_step_inisde_emission_medium) {

            this->continue_integration = false;

        }
        else {

            this->continue_integration = true;

        }

    }
    else {

        this->continue_integration = false;
        this->N_steps_rejected++;

    }    


    this->p_Step_controller->update_step(this->e_Active_integrator, New_State_vector_main[e_r], Is_inside_emission_medium);

    if (this->p_Init_conditions->Metric_parameters.e_Spacetime == Janis_Newman_Winicour and this->p_Init_conditions->Metric_parameters.JNW_Gamma_Parameter < 0.5) {

        if (New_State_vector_main[e_r] - 2 * this->p_Init_conditions->Metric_parameters.Mass / this->p_Init_conditions->Metric_parameters.JNW_Gamma_Parameter < 1e-2 and this->Force_scatter) {

            New_State_vector_main[e_p_r] = -1 * abs(New_State_vector_main[e_p_r]);

            if (this->continue_integration) { this->Force_scatter = false; }

        }
    }

    if (this->continue_integration) {

        memcpy(this->Current_Dynamic_state, New_State_vector_main, e_Dynamic_state_size * sizeof(double));

        this->Update_ray_log(New_State_vector_main);
        this->Update_optical_depth();
        this->Update_debug_log();
        
        this->N_steps_rejected = 0;
        this->NaN_checker_count = 0;

    }

}

void Geodesic_Integrator_class::Update_ray_log(const double* const New_State_vector) {

    this->p_Ray_log_struct->Log_offset += 1;
    size_t& log_offset = this->p_Ray_log_struct->Log_offset;

    memcpy(&this->p_Ray_log_struct->Ray_path_log_global[log_offset * e_Full_state_size], New_State_vector, e_Dynamic_state_size * sizeof(double));

    this->p_Ray_log_struct->Ray_path_log_global[e_step + log_offset * e_Full_state_size] = this->p_Step_controller->previous_step;
    this->p_Ray_log_struct->Ray_path_log_global[e_ray_affine_param + log_offset * e_Full_state_size] = this->p_Ray_log_struct->Ray_path_log_global[e_ray_affine_param + (log_offset - 1) * e_Full_state_size] - this->p_Step_controller->previous_step;
    
    this->p_Spacetime->Convert_global_to_local_coords(&this->p_Ray_log_struct->Ray_path_log_global[log_offset * e_Full_state_size],
                                                      &this->p_Ray_log_struct->Ray_path_log_global[log_offset * e_Full_state_size],
                                                      &this->p_Ray_log_struct->Ray_path_log_local.get()[log_offset * e_Full_state_size],
                                                      e_Full_State_Vector);

}

void Geodesic_Integrator_class::Update_debug_log() {

    this->RK_Integrator_debug_log->N_steps_rejected[this->p_Ray_log_struct->Log_offset] = this->N_steps_rejected;
    this->RK_Integrator_debug_log->State_error_history[this->p_Ray_log_struct->Log_offset] = this->p_Step_controller->current_err;

}

void Geodesic_Integrator_class::Propagate_ray() {

    switch (this->e_Active_integrator) {

    case ESDIRK54:

        this->Run_ESDIRK54();
        this->Check_integration_complete_status();

        break;

    default:

        this->Run_Explicit_Runge_Kutta();
        this->Check_integration_complete_status();

        if (this->Step_too_small or this->NaN_checker_count > 5) {
            
            this->e_Active_integrator = ESDIRK54;
            this->p_Ray_log_struct->Log_offset = 0;
        
            std::cout << "Switching to ESDIRK54... \n";

        }

        break;

    }

}

bool Geodesic_Integrator_class::Locate_event(Event_detection_enums e_Event, double* const Global_State_at_Event, double* const Local_State_at_Event) {

    int Event_idx{};
    double Event_target{};

    const double* const Current_State = this->get_current_State_Vector_global();
    const double* const Prev_State = this->get_previous_State_Vector_global();

    switch (e_Event) {

    case Equatorial_crossing:

        if ((Current_State[e_theta] - M_PI_2) * (Prev_State[e_theta] - M_PI_2) > 0) { return false; }

        Event_idx = e_theta;
        Event_target = M_PI_2;

        break;

    case Celestial_sphere_crossing:

        if ((abs(Current_State[e_r]) - this->p_Init_conditions->Metric_parameters.Scattering_radius) *
            (abs(Prev_State[e_r]) - this->p_Init_conditions->Metric_parameters.Scattering_radius) > 0) {
            return false;
        }

        Event_idx = e_r;
        Event_target = this->p_Init_conditions->Metric_parameters.Scattering_radius;

        break;

    default:

        throw std::runtime_error("Unsupported event type in Geodesic_Integrator_class::Locate_event()!");

    }

    /* ========== Construct the polynomial coefficients - ax^3 + bx^2 + cx + d ========== */

    /* Not every method is a first-same-as-last method, so the current RHS needs to be evaluated, rather than read off from Intermediate_RHS_log.
       I overwride the last entries in Intermediate_RHS_log. */
    int RK_size = RK78_size;

    if (RK54 == e_Active_integrator) { RK_size = RK54_size; }
    else if (ESDIRK54 == e_Active_integrator) { RK_size = ESDIRK54_size; }

    double Current_RHS[e_Dynamic_state_size]{};
    this->p_Spacetime->get_EOM(Current_State, this->Intermediate_RHS_log + (RK_size - 1) * e_Dynamic_state_size);
    memcpy(Current_RHS, this->Intermediate_RHS_log + (RK_size - 1) * e_Dynamic_state_size, e_Dynamic_state_size * sizeof(double));

    const double* const &Prev_RHS = this->Intermediate_RHS_log;

    const double& step = -this->p_Step_controller->previous_step;

    const double a_coeff = step * (Current_RHS[Event_idx] + Prev_RHS[Event_idx]) - 2 * (Current_State[Event_idx] - Prev_State[Event_idx]);
    const double b_coeff = 3 * (Current_State[Event_idx] - Prev_State[Event_idx]) - step * (2 * Prev_RHS[Event_idx] + Current_RHS[Event_idx]);
    const double c_coeff = step * Prev_RHS[Event_idx];
    const double d_coeff = Prev_State[Event_idx] - Event_target;

    std::complex<double> Cubic_roots[3]{};
    get_cubic_polynomial_roots(a_coeff, b_coeff, c_coeff, d_coeff, Cubic_roots);

    double Event_interp_param = -1;

    for (int idx = 0; idx < 3; idx++) {

        if (abs(Cubic_roots[idx].imag()) < 1e-10 and Cubic_roots[idx].real() >= 0 and Cubic_roots[idx].real() <= 1) {

            Event_interp_param = Cubic_roots[idx].real();
            break;

        }

    }

    /* This check will pass only if no real roots lie in the interval [0, 1], which should never happen, but sometimes does for some reason. */
    if (Event_interp_param < 0 or Event_interp_param > 1) { return false; }

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        Global_State_at_Event[idx] = this->get_dense_output(Event_interp_param, State_enums(idx), true);

    }

    Global_State_at_Event[e_step] = abs(Event_interp_param * step);
    Global_State_at_Event[e_ray_affine_param] = Prev_State[e_ray_affine_param] + Event_interp_param * step;

    this->p_Spacetime->Convert_global_to_local_coords(Global_State_at_Event, Global_State_at_Event, Local_State_at_Event, e_Full_State_Vector);

    return true;

}

const double Geodesic_Integrator_class::get_dense_output(const double Param, const State_enums idx, bool Is_current_RHS_evaluated) const {

    /* The source for this implementation is https://mezbanhabibi.ir/wp-content/uploads/2020/01/ordinary-differential-equations-vol.1.-Nonstiff-problems.pdf, equation (6.7) */

    int RK_size = RK78_size;

    if (RK54 == e_Active_integrator) { RK_size = RK54_size; }
    else if (ESDIRK54 == e_Active_integrator) { RK_size = ESDIRK54_size; }

    const double* const Current_State = this->get_current_State_Vector_global();
    const double* const Prev_State = this->get_previous_State_Vector_global();

    double Current_RHS[e_Dynamic_state_size]{};
    memcpy(Current_RHS, this->Intermediate_RHS_log + (RK_size - 1) * e_Dynamic_state_size, e_Dynamic_state_size * sizeof(double));
    
    if (!Is_current_RHS_evaluated) {

        this->p_Spacetime->get_EOM(Current_State, Current_RHS);

    }

    const double* const &Prev_RHS = this->Intermediate_RHS_log;
    const double& step = -this->p_Step_controller->previous_step;

    return (1. - Param) * Prev_State[idx] + Param * Current_State[idx] + Param * (Param - 1.) * ((1. - 2. * Param) * (Current_State[idx] - Prev_State[idx]) + (Param - 1.) * step * Prev_RHS[idx] + Param * step * Current_RHS[idx]);

}

void Geodesic_Integrator_class::Check_integration_complete_status() {

    /* This needs to use the internal dynamic state, because it is kept in "global coordainates" (which so far only affects the wormhole). */
    this->Normal_termination_condition = this->p_Spacetime->terminate_integration(this->Current_Dynamic_state);

    this->Max_affine_param_reached      = std::abs(this->get_current_State_Vector_global()[e_ray_affine_param]) >= this->Max_affine_param;
    this->Max_integration_count_reached = this->p_Ray_log_struct->Log_offset >= this->Max_integration_count - 1;
    this->Step_too_small                = this->p_Step_controller->step < std::numeric_limits<double>::min();

    if (this->Max_affine_param_reached) { 

        switch (this->e_Active_integrator) {

        case ESDIRK54:

            std::cout << "Max affine parameter value reachedl with ESDIRK54! \n";

            break;

        default:

            std::cout << "Max affine parameter value reached with the explicit Runge-Kutta! \n";

            break;

        }
    }

    if (this->Max_integration_count_reached) { 

        switch (this->e_Active_integrator) {

        case ESDIRK54:

            std::cout << "Max iterations reached with ESDIRK54! \n";

            break;

        default:

            std::cout << "Max iterations reached with the explicit Runge-Kutta! \n";

            break;

        }
    }

    if (this->Step_too_small) {

        switch (this->e_Active_integrator) {

        case ESDIRK54:

            std::cout << "Step too small with ESDIRK54! \n";

            break;

        default:

            std::cout << "Step too small with the explicit Runge-Kutta! \n";

            break;

        }
    }
    
    this->integration_complete = this->Normal_termination_condition or this->Max_affine_param_reached or (*this->Current_Optical_Depth > 100);

    if (ESDIRK54 == this->e_Active_integrator) {

        this->integration_complete = this->integration_complete or this->Step_too_small or this->Max_integration_count_reached;

    }

}

const double* const Geodesic_Integrator_class::get_current_State_Vector_global() const {

    return &this->p_Ray_log_struct->Ray_path_log_global[this->p_Ray_log_struct->Log_offset * e_Full_state_size];

}

const double* const Geodesic_Integrator_class::get_previous_State_Vector_global() const {

    if (this->p_Ray_log_struct->Log_offset > 0) {

        return &this->p_Ray_log_struct->Ray_path_log_global[(this->p_Ray_log_struct->Log_offset - 1) * e_Full_state_size];
    }
    else {

        return &this->p_Ray_log_struct->Ray_path_log_global[this->p_Ray_log_struct->Log_offset * e_Full_state_size];

    }
    
}

const double* const Geodesic_Integrator_class::get_current_State_Vector_local() const {

    return &this->p_Ray_log_struct->Ray_path_log_local.get()[this->p_Ray_log_struct->Log_offset * e_Full_state_size];

}

const double* const Geodesic_Integrator_class::get_previous_State_Vector_local() const {

    if (this->p_Ray_log_struct->Log_offset > 0) {

        return &this->p_Ray_log_struct->Ray_path_log_local.get()[(this->p_Ray_log_struct->Log_offset - 1) * e_Full_state_size];
    }
    else {

        return &this->p_Ray_log_struct->Ray_path_log_local.get()[this->p_Ray_log_struct->Log_offset * e_Full_state_size];

    }

}