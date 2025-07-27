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

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        Error_scale[idx] = this->Parameters.RK_78_accuracy + std::fabs(State_Vector[idx]) * this->Parameters.RK_78_accuracy;

    }

    double Total_State_Error{};

    for (int idx = 0; idx < e_Dynamic_state_size; idx++) {

        Total_State_Error += (State_Error_Vector[idx] / Error_scale[idx]) * (State_Error_Vector[idx] / Error_scale[idx]);

    }

    Total_State_Error /= e_Dynamic_state_size;

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

    // Construct the initial state sectors
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

void Integrator_class::Update_ray_log(const double* const New_State_vector) {

    this->p_Ray_log_struct->Log_offset += 1;
    int& log_offset = this->p_Ray_log_struct->Log_offset;
    const double& R_throat = this->p_Init_conditions->Metric_parameters.R_throat;

    memcpy(&this->p_Ray_log_struct->Ray_path_log[log_offset * e_Full_state_size], New_State_vector, e_Dynamic_state_size * sizeof(double));

    this->p_Ray_log_struct->Ray_path_log[e_step + log_offset * e_Full_state_size] = this->p_Step_controller->previous_step;

    if (log_offset > 0) {
        
        this->p_Ray_log_struct->Ray_path_log[e_affine_param + log_offset * e_Full_state_size] = this->p_Ray_log_struct->Ray_path_log[(log_offset - 1) * e_Full_state_size + e_affine_param] - this->p_Step_controller->previous_step;

    }
    else {

        this->p_Ray_log_struct->Ray_path_log[e_affine_param + log_offset * e_Full_state_size] = -this->p_Init_conditions->Integrator_params.Init_stepzie;

    }
    
    // The wormhole metric works with a "global" radial coordinate, that goes negative on the other side of the throat.
    // The emission model can't work with this coordinate, so I log the normal spherical radial coordinate instead. 
    if (Wormhole == this->p_Init_conditions->Metric_parameters.e_Spacetime) {

        this->p_Ray_log_struct->Ray_path_log[e_r + log_offset * e_Full_state_size] = sqrt(New_State_vector[e_r] * New_State_vector[e_r] + R_throat * R_throat);

    }

}

void Integrator_class::Propagate_ray() {

    this->Run_RK78();
    this->Check_integration_complete_status();

}

void Integrator_class::Check_integration_complete_status() {

    // References for the sake of readability
    const int& log_offset = this->p_Ray_log_struct->Log_offset;
    const double* State_Vector = &this->p_Ray_log_struct->Ray_path_log[log_offset * e_Full_state_size];

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