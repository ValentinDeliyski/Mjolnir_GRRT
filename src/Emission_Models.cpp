#include "Emission_models.h"

double Emission_models_class::get_electron_pitch_angle(const double* const B_field_coord_frame, 
                                                       const double* const Plasma_velocity,
                                                       const double* const Local_State_Vector) {

    Metric_type s_Metric = this->p_Sim_Context->p_Spacetime->get_local_metric(Local_State_Vector);

    double Wave_vec_dot_Plasma_vec = dot_product(Local_State_Vector + e_p_t, Plasma_velocity, 4);
    double Wave_vec_dot_B_field    = dot_product(Local_State_Vector + e_p_t, B_field_coord_frame, 4);

    double B_field_norm_squared   = get_4vec_dot_product(B_field_coord_frame, B_field_coord_frame, s_Metric.Metric, Contravariant);
    double B_field_dot_Plasma_vel = get_4vec_dot_product(B_field_coord_frame, Plasma_velocity, s_Metric.Metric, Contravariant);

    double cos_angle = 1.0; 

    if (!isinf(1.0 / Wave_vec_dot_Plasma_vec) and !isinf(1.0 / B_field_norm_squared)) {

        cos_angle = Wave_vec_dot_B_field + B_field_dot_Plasma_vel * Wave_vec_dot_Plasma_vec;
        cos_angle /= fabs(Wave_vec_dot_Plasma_vec);
        cos_angle /= sqrt(B_field_norm_squared + B_field_dot_Plasma_vel * B_field_dot_Plasma_vel);
           
    }

    if (fabs(cos_angle) <= 1.0) {

        return acos(cos_angle);

    }
    else {

        return acos(copysign(1.0 - 1e-10, cos_angle));

    }
}

/* =============================================== Thermal synchrotron Transfer Functions =============================================== */

void Emission_models_class::get_thermal_synchrotron_transfer_functions(const double* const Local_State_Vector,
                                                                       const Emission_medium_state_type* const p_Emission_medium_state,
                                                                       Transfer_functions_type* const p_Transfer_functions) {

    /* === Zero out the transfer functions just in case === */
    memset(p_Transfer_functions, 0, sizeof(Transfer_functions_type));

    double redshift = get_redshift(Local_State_Vector, p_Emission_medium_state->Plasma_Velocity, this->p_Sim_Context);

    /* Check weather redshift is numerically OK to use in the transfer functions. */
    if (isinf(redshift) or isnan(redshift) or isinf(1.0 / redshift)) { return; }

    /* The dimensionless electron temperature. */
    double const T_electron_dim = BOLTZMANN_CONST_CGS * p_Emission_medium_state->Temperature / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;

    /* The cyclotron frequency. */
    double const f_cyclo = Q_ELECTRON_CGS * p_Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    /* The "averaged" rescaled (by a factor of 27 / 4) critical frequency (without the sin(theta) term.
       That gets added on later from a pre-computed table). */
    double const f_s_no_sin = 2. / 9 * f_cyclo * T_electron_dim * T_electron_dim;

    /* Check weather the rescaled critical frequency f_s is numerically OK to use in the transfer functions. */
    if (isinf(f_s_no_sin) or isnan(f_s_no_sin) or isinf(1.0 / f_s_no_sin)) { return; }

    Thermal_transfer_f_arguments_type Transfer_args_uncorrected{};

    /* Observation Frequency */
    double const obs_frequency = this->p_Sim_Context->p_Init_Conditions->Observer_params.obs_frequency;

    /* Compute all the wierd powers of X outside the pitch angle averaging loop. */
    Transfer_args_uncorrected.X      = obs_frequency / f_s_no_sin / redshift;
    Transfer_args_uncorrected.sqrt_X = sqrt(Transfer_args_uncorrected.X);
    Transfer_args_uncorrected.cbrt_X = cbrt(Transfer_args_uncorrected.X);
    Transfer_args_uncorrected.X_to_0_p_5175 = pow(Transfer_args_uncorrected.X, 0.5175);
    Transfer_args_uncorrected.X_to_0_p_6    = pow(Transfer_args_uncorrected.X, 0.6);
    Transfer_args_uncorrected.X_to_0_p_7515 = pow(Transfer_args_uncorrected.X, 0.7515);

    /* Compute the werid power of the electron temperature outside the pitch angle averaging loop. */
    Transfer_args_uncorrected.T_electron_dim          = T_electron_dim;
    Transfer_args_uncorrected.T_electron_dim_to_24_25 = pow(T_electron_dim, 24. / 25);

    /* We also need the current photon frequency to evalauate the emission functions. */
    Transfer_args_uncorrected.frequency = obs_frequency / redshift;

    /* Reference for the sake of readabiity. */
    int& Num_Samples_to_avg = p_Sim_Context->p_Init_Conditions->Emission_pitch_angle_samples_to_average;

    /* This structcs holds the transfer function args, corrected for the electron pitch angle. I am setting it equal to the uncorrected one
        so they copy over the variables that don't depend on the pitch angle. The rest get corrected inside the averaging loop. */
    Thermal_transfer_f_arguments_type Transfer_args_corrected = Transfer_args_uncorrected;

    if (this->p_Sim_Context->p_Init_Conditions->Average_electron_pitch_angle) {

        /* ============ This loop averages over the emission pitch angle, which it gets from a pre-computed table ============ */

        for (int averaging_idx = 1; averaging_idx < Num_Samples_to_avg; averaging_idx++) {

            /* References for the sake of readabiity. */
            double& sin_pitch_angle = this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[averaging_idx];
            double& cos_pitch_angle = this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[averaging_idx];

            Transfer_args_corrected.X      = Transfer_args_uncorrected.X / sin_pitch_angle;
            Transfer_args_corrected.sqrt_X = Transfer_args_uncorrected.sqrt_X * this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
            Transfer_args_corrected.cbrt_X = Transfer_args_uncorrected.cbrt_X * this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[averaging_idx];
            Transfer_args_corrected.X_to_0_p_5175 = Transfer_args_uncorrected.X_to_0_p_5175 * this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_5175[averaging_idx];
            Transfer_args_corrected.X_to_0_p_6    = Transfer_args_uncorrected.X_to_0_p_6 * this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_6[averaging_idx];
            Transfer_args_corrected.X_to_0_p_7515 = Transfer_args_uncorrected.X_to_0_p_7515 * this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_7515[averaging_idx];

            Transfer_args_corrected.sin_pitch_angle = sin_pitch_angle;
            Transfer_args_corrected.cos_pitch_angle = cos_pitch_angle;

            Transfer_functions_type temp_Transfer_functions{};

            this->get_synchrotron_transfer_fit_functions(e_Thermal_ensamble, p_Emission_medium_state, &Transfer_args_corrected, &temp_Transfer_functions);

            // The U component is 0 by definition
            p_Transfer_functions->Emission_functions[I] += temp_Transfer_functions.Emission_functions[I] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Emission_functions[Q] += temp_Transfer_functions.Emission_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Emission_functions[V] += temp_Transfer_functions.Emission_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

            // The U component is 0 by definition
            p_Transfer_functions->Absorbtion_functions[I] += temp_Transfer_functions.Absorbtion_functions[I] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Absorbtion_functions[Q] += temp_Transfer_functions.Absorbtion_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Absorbtion_functions[V] += temp_Transfer_functions.Absorbtion_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

            // The I and U components are 0 by definition
            p_Transfer_functions->Faraday_functions[Q] += temp_Transfer_functions.Faraday_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Faraday_functions[V] += temp_Transfer_functions.Faraday_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

        }
    }
    else {

        /* The magnetic field is the one measured by a comoving with the plasma observer, but expressed in the cooridante frame */
        double pitch_angle = get_electron_pitch_angle(p_Emission_medium_state->Magnetic_fields.B_field_plasma_frame, p_Emission_medium_state->Plasma_Velocity, Local_State_Vector);
        double sin_pitch_angle = sin(pitch_angle);

        double one_over_sqrt_sin = 1.0 / sqrt(sin_pitch_angle);
        double one_over_cbrt_sin = 1.0 / cbrt(sin_pitch_angle);
        double one_over_sin_to_0_p_5175 = 1.0 / pow(sin_pitch_angle, 0.5175);
        double one_over_sin_to_0_p_6    = 1.0 / pow(sin_pitch_angle, 0.6);
        double one_over_sin_to_0_p_7515 = 1.0 / pow(sin_pitch_angle, 0.7515);

        Transfer_args_corrected.X      = Transfer_args_uncorrected.X / sin_pitch_angle;
        Transfer_args_corrected.sqrt_X = Transfer_args_uncorrected.sqrt_X * one_over_sqrt_sin;
        Transfer_args_corrected.cbrt_X = Transfer_args_uncorrected.cbrt_X * one_over_cbrt_sin;
        Transfer_args_corrected.X_to_0_p_5175 = Transfer_args_uncorrected.X_to_0_p_5175 * one_over_sin_to_0_p_5175;
        Transfer_args_corrected.X_to_0_p_6    = Transfer_args_uncorrected.X_to_0_p_6 * one_over_sin_to_0_p_6;
        Transfer_args_corrected.X_to_0_p_7515 = Transfer_args_uncorrected.X_to_0_p_7515 * one_over_sin_to_0_p_7515;

        Transfer_args_corrected.sin_pitch_angle = sin_pitch_angle;
        Transfer_args_corrected.cos_pitch_angle = cos(pitch_angle);

        this->get_synchrotron_transfer_fit_functions(e_Thermal_ensamble, p_Emission_medium_state, &Transfer_args_corrected, p_Transfer_functions);

    }

    /* Account for the relativistic doppler effect via the redshift. */
    for (int stokes_idx = 0; stokes_idx < e_Stokes_param_num; stokes_idx++) {

        p_Transfer_functions->Emission_functions[stokes_idx] *= redshift * redshift;
        p_Transfer_functions->Faraday_functions[stokes_idx] /= redshift;
        p_Transfer_functions->Absorbtion_functions[stokes_idx] /= redshift;
    }

}
 
/* ========================================== Kappa synchrotron Transfer Functions ========================================== */

void Emission_models_class::get_kappa_synchrotron_transfer_functions(const double* const Local_State_Vector,
                                                                           const Emission_medium_state_type* const p_Emission_medium_state,
                                                                           Transfer_functions_type* const p_Transfer_functions){

    /* === Zero out the transfer functions just in case === */
    memset(p_Transfer_functions, 0, sizeof(Transfer_functions_type));

    const double redshift = get_redshift(Local_State_Vector, p_Emission_medium_state->Plasma_Velocity, this->p_Sim_Context);

    if (isinf(redshift) or isnan(redshift) or isinf(1.0 / redshift)) { return; }

    /* Dimensionless Electron Temperature */
    double T_electron_dim = BOLTZMANN_CONST_CGS * p_Emission_medium_state->Temperature / M_ELECTRON_CGS / C_LIGHT_CGS / C_LIGHT_CGS;

    /* Cyclotron Frequency */
    double f_cyclo = Q_ELECTRON_CGS * p_Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    /* The "averaged" critical frequency (without the sin(theta) term - that gets added on later from a pre-computed table) */
    double f_k_no_sin = f_cyclo * (this->s_Emission_params.Kappa * T_electron_dim) * (this->s_Emission_params.Kappa * T_electron_dim);

    if (isinf(f_k_no_sin) or isnan(f_k_no_sin) or isinf(1.0 / f_k_no_sin)) { return; }

    /* Observation frequency */
    double& obs_frequency = this->p_Sim_Context->p_Init_Conditions->Observer_params.obs_frequency;

    Kappa_transfer_f_arguments_type Transfer_args_uncorrected{};

    /* Populate the angle-uncorrected transfer arguments struct. */
    Transfer_args_uncorrected.X              = obs_frequency / f_k_no_sin / redshift;
    Transfer_args_uncorrected.sqrt_X         = sqrt(Transfer_args_uncorrected.X);
    Transfer_args_uncorrected.cbrt_X         = cbrt(Transfer_args_uncorrected.X);
    Transfer_args_uncorrected.X_to_7_over_20 = pow(Transfer_args_uncorrected.X, 7. / 20);
    Transfer_args_uncorrected.kappa          = this->s_Emission_params.Kappa;
    Transfer_args_uncorrected.T_electron_dim = T_electron_dim;
    Transfer_args_uncorrected.frequency      = obs_frequency / redshift;

    int& Num_Samples_to_avg = this->p_Sim_Context->p_Init_Conditions->Emission_pitch_angle_samples_to_average;

    /* This structcs holds the transfer function args, corrected for the electron pitch angle. I am setting it equal to the uncorrected one
        so they copy over the variables that don't depend on the pitch angle. The rest get corrected inside the averaging loop. */
    Kappa_transfer_f_arguments_type Transfer_args_corrected = Transfer_args_uncorrected;

    if (this->p_Sim_Context->p_Init_Conditions->Average_electron_pitch_angle) {

        /* ============ This loop averages over the emission pitch angle, which it gets from a pre-computed table ============ */

        for (int averaging_idx = 1; averaging_idx < Num_Samples_to_avg; averaging_idx++) {

            double& sin_pitch_angle = this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[averaging_idx];
            double& cos_pitch_angle = this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[averaging_idx];

            Transfer_args_corrected.X                  = Transfer_args_uncorrected.X / sin_pitch_angle;
            Transfer_args_corrected.sqrt_X             = Transfer_args_uncorrected.sqrt_X * this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[averaging_idx];
            Transfer_args_corrected.cbrt_X             = Transfer_args_uncorrected.cbrt_X * this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[averaging_idx];
            Transfer_args_corrected.X_to_7_over_20     = Transfer_args_uncorrected.X_to_7_over_20 * this->s_Precomputed_e_pitch_angles.one_over_sin_to_7_over_20[averaging_idx];
            Transfer_args_corrected.sin_emission_angle = sin_pitch_angle;
            Transfer_args_corrected.cos_emission_angle = cos_pitch_angle;

            Transfer_functions_type temp_Transfer_functions{};

            this->get_synchrotron_transfer_fit_functions(e_Kappa_ensamble, p_Emission_medium_state, &Transfer_args_corrected, &temp_Transfer_functions);

            // The U component is 0 by definition
            p_Transfer_functions->Emission_functions[I] += temp_Transfer_functions.Emission_functions[I] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Emission_functions[Q] += temp_Transfer_functions.Emission_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Emission_functions[V] += temp_Transfer_functions.Emission_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

            // The U component is 0 by definition
            p_Transfer_functions->Absorbtion_functions[I] += temp_Transfer_functions.Absorbtion_functions[I] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Absorbtion_functions[Q] += temp_Transfer_functions.Absorbtion_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Absorbtion_functions[V] += temp_Transfer_functions.Absorbtion_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

            // The I and U components are 0 by definition
            p_Transfer_functions->Faraday_functions[Q] += temp_Transfer_functions.Faraday_functions[Q] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;
            p_Transfer_functions->Faraday_functions[V] += temp_Transfer_functions.Faraday_functions[V] * sin_pitch_angle * M_PI / Num_Samples_to_avg / 2;

        }
    }
    else {

        /* The magnetic field is the one measured by a comoving with the plasma observer, but expressed in the cooridante frame */
        double pitch_angle = get_electron_pitch_angle(p_Emission_medium_state->Magnetic_fields.B_field_plasma_frame, p_Emission_medium_state->Plasma_Velocity, Local_State_Vector);
        double sin_pitch_angle = sin(pitch_angle);

        double one_over_sqrt_sin    = 1. / sqrt(sin_pitch_angle);
        double one_over_cbrt_sin    = 1. / cbrt(sin_pitch_angle);
        double one_over_7_to_20_sin = 1. / pow(sin_pitch_angle, 7. / 20);

        Transfer_args_corrected.X                  = Transfer_args_uncorrected.X / sin_pitch_angle;
        Transfer_args_corrected.sqrt_X             = Transfer_args_uncorrected.sqrt_X * one_over_sqrt_sin;
        Transfer_args_corrected.cbrt_X             = Transfer_args_uncorrected.cbrt_X * one_over_cbrt_sin;
        Transfer_args_corrected.X_to_7_over_20     = Transfer_args_uncorrected.X_to_7_over_20 * one_over_7_to_20_sin;
        Transfer_args_corrected.sin_emission_angle = sin_pitch_angle;
        Transfer_args_corrected.cos_emission_angle = cos(pitch_angle);

        this->get_synchrotron_transfer_fit_functions(e_Kappa_ensamble, p_Emission_medium_state, &Transfer_args_corrected, p_Transfer_functions);

    }

    /* Account for the relativistic doppler effect via the redshift */
    for (int stokes_idx = 0; stokes_idx < e_Stokes_param_num; stokes_idx++) {

        p_Transfer_functions->Emission_functions[stokes_idx] *= redshift * redshift;
        p_Transfer_functions->Faraday_functions[stokes_idx]  /= redshift;
        p_Transfer_functions->Absorbtion_functions[stokes_idx] /= redshift;
    }

}

/* ========================================== Phenomenological synchrotron Transfer Functions ========================================== */

void Emission_models_class::get_phenomenological_synchrotron_functions(const double* const Local_State_Vector,
                                                                       const Emission_medium_state_type* const p_Emission_medium_state,
                                                                       Transfer_functions_type* const p_Transfer_functions) {

    /* === Zero out the transfer functions just in case. === */
    memset(p_Transfer_functions, 0, sizeof(Transfer_functions_type));

    Phenomenological_transfer_f_arguments_type Transfer_args{};

    Transfer_args.redshift = get_redshift(Local_State_Vector, p_Emission_medium_state->Plasma_Velocity, this->p_Sim_Context);

    if (isinf(Transfer_args.redshift) or isnan(Transfer_args.redshift) or isinf(1.0 / Transfer_args.redshift)) { return; }

    Transfer_args.frequency = this->p_Sim_Context->p_Init_Conditions->Observer_params.obs_frequency / Transfer_args.redshift;
    Transfer_args.f_cyclo = Q_ELECTRON_CGS * p_Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    this->get_synchrotron_transfer_fit_functions(e_Phenomenological_ensamble, p_Emission_medium_state, &Transfer_args, p_Transfer_functions);

    /* Account for the relativistic doppler effect via the redshift. */
    for (int stokes_idx = 0; stokes_idx < e_Stokes_param_num; stokes_idx++) {

        p_Transfer_functions->Emission_functions[stokes_idx] *= Transfer_args.redshift * Transfer_args.redshift;
        p_Transfer_functions->Faraday_functions[stokes_idx] /= Transfer_args.redshift;
        p_Transfer_functions->Absorbtion_functions[stokes_idx] /= Transfer_args.redshift;
    }
}

/* ============================================ Main "Selector" For The Transfer Functions ============================================ */

void Emission_models_class::get_radiative_transfer_functions(const double* const Local_State_Vector,
                                                             const Emission_medium_enums Emission_medium,
                                                             Transfer_functions_type* const p_Transfer_functions) {

    /* === Zero out the transfer functions just in case. === */
    memset(p_Transfer_functions, 0, sizeof(Transfer_functions_type));

    if (e_Novikov_Thorne == this->p_Disk_Model->s_Disk_params.e_Disk_model) { return; }

    Emission_medium_state_type Emission_medium_state{};

    /* This variable exist for the case where the hotspot and disk are in "Thermalized" mode. */
    Emission_medium_state_type Hotspot_state{};

    Metric_type Metric = this->p_Sim_Context->p_Spacetime->get_local_metric(Local_State_Vector);

    bool Is_inside_hotspot = false;
    bool Is_inside_disk = false;

    /* This function call populates the density and temperature values for the hotspot - this is why they are not populated along with the magnetic field parameters. */
    Is_inside_hotspot = this->p_Hotspot_Model->is_inside_hotspot(Local_State_Vector, &Hotspot_state);

    switch (Emission_medium) {

    case Disk:

       /* This function call populates the density and temperature values for the disk - this is why they are not populated along with the magnetic field parameters. */
        Is_inside_disk = this->p_Disk_Model->is_inside_disk(Local_State_Vector, &Emission_medium_state);

        if (!Is_inside_disk) { return; };

        memcpy(Emission_medium_state.Plasma_Velocity, this->p_Disk_Model->get_disk_velocity(Local_State_Vector), 4 * sizeof(double));

        if (Is_inside_hotspot) {

            /* -------- We are inside the hotspot - set the magnetic field properties of the disk to be equal to those of the hotspot. -------- */
            this->p_Hotspot_Model->get_magnetic_field(Local_State_Vector, &Metric, &Emission_medium_state);

        }
        else {

            /* -------- We are outside the hotspot - set the magnetic field properties of the disk to be equal to their "background" values. -------- */
            this->p_Disk_Model->get_magnetic_field(Local_State_Vector, &Metric, &Emission_medium_state);
         
        }

        Emission_medium_state.Ensamble_type = this->p_Disk_Model->s_Disk_params.Ensamble_type;

        break;

    case Hotspot:

        if (!Is_inside_hotspot) { return; };

        Emission_medium_state.Density = Hotspot_state.Density;
        Emission_medium_state.Temperature = Hotspot_state.Temperature;
        Emission_medium_state.Ensamble_type = this->p_Hotspot_Model->s_Hotspot_params.Ensamble_type;

        memcpy(Emission_medium_state.Plasma_Velocity, this->p_Hotspot_Model->get_hotspot_velocity(false, Local_State_Vector), 4 * sizeof(double));
        this->p_Hotspot_Model->get_magnetic_field(Local_State_Vector, &Metric, &Emission_medium_state);

        break;

    default:

        throw std::runtime_error("Unsupported emissison medium - something broke in the get_radiative_transfer_functions function!");

    }
 
    switch (Emission_medium_state.Ensamble_type) {

    case(e_Phenomenological_ensamble):

        this->get_phenomenological_synchrotron_functions(Local_State_Vector, &Emission_medium_state, p_Transfer_functions);
        break;

    case(e_Kappa_ensamble):

        this->get_kappa_synchrotron_transfer_functions(Local_State_Vector, &Emission_medium_state, p_Transfer_functions);
        break;

    case(e_Debug_constant_functions):

        this->get_debug_synchrotron_functions(p_Transfer_functions);
        break;

    default:

        this->get_thermal_synchrotron_transfer_functions(Local_State_Vector, &Emission_medium_state, p_Transfer_functions);
        break;
    }

}

void Emission_models_class::get_debug_synchrotron_functions(Transfer_functions_type* p_Transfer_functions) const {

    p_Transfer_functions->Emission_functions[I] = this->s_Emission_params.Debug_j_I_value;
    p_Transfer_functions->Emission_functions[Q] = this->s_Emission_params.Debug_j_Q_value;
    p_Transfer_functions->Emission_functions[U] = this->s_Emission_params.Debug_j_U_value;
    p_Transfer_functions->Emission_functions[V] = this->s_Emission_params.Debug_j_V_value;

    p_Transfer_functions->Absorbtion_functions[I] = this->s_Emission_params.Debug_alpha_I_value;
    p_Transfer_functions->Absorbtion_functions[Q] = this->s_Emission_params.Debug_alpha_Q_value;
    p_Transfer_functions->Absorbtion_functions[U] = this->s_Emission_params.Debug_alpha_U_value;
    p_Transfer_functions->Absorbtion_functions[V] = this->s_Emission_params.Debug_alpha_V_value;

    p_Transfer_functions->Faraday_functions[I] = this->s_Emission_params.Debug_rho_I_value;
    p_Transfer_functions->Faraday_functions[Q] = this->s_Emission_params.Debug_rho_Q_value;
    p_Transfer_functions->Faraday_functions[U] = this->s_Emission_params.Debug_rho_U_value;
    p_Transfer_functions->Faraday_functions[V] = this->s_Emission_params.Debug_rho_V_value;

}

void Emission_models_class::get_synchrotron_transfer_fit_functions(const Ensamble_enums e_Ensamble_type,
                                                                           const Emission_medium_state_type* const p_Emission_medium_state,
                                                                           const void* const p_Transfer_args,
                                                                           Transfer_functions_type* const p_Transfer_functions) const {

    /* The dimensionless frequency needs to get extracted from the p_Transfer_args pointer, but it first needs to be recast to not-void.
       This happens in the scopes of the switch statemeent below, so I create a variable here to store it. */
    double  frequency{};

    switch (e_Ensamble_type) {

    case e_Thermal_ensamble:

        this->get_thermal_synchrotron_emission_fit_functions(static_cast<const Thermal_transfer_f_arguments_type*>(p_Transfer_args), p_Transfer_functions->Emission_functions);
        this->get_thermal_synchrotron_absorbtion_fit_functions(static_cast<const Thermal_transfer_f_arguments_type*>(p_Transfer_args), p_Emission_medium_state, p_Transfer_functions->Emission_functions, p_Transfer_functions->Absorbtion_functions);
        this->get_thermal_synchrotron_Faraday_fit_functions(static_cast<const Thermal_transfer_f_arguments_type*>(p_Transfer_args), p_Transfer_functions->Faraday_functions);

        frequency = static_cast<const Thermal_transfer_f_arguments_type*>(p_Transfer_args)->frequency;

        break;

    case e_Kappa_ensamble:

        this->get_kappa_synchrotron_emission_fit_functions(static_cast<const Kappa_transfer_f_arguments_type*>(p_Transfer_args), p_Transfer_functions->Emission_functions);
        this->get_kappa_synchrotron_absorbtion_fit_functions(static_cast<const Kappa_transfer_f_arguments_type*>(p_Transfer_args), p_Transfer_functions->Absorbtion_functions);
        this->get_kappa_synchrotron_Faraday_fit_functions(static_cast<const Kappa_transfer_f_arguments_type*>(p_Transfer_args), p_Transfer_functions->Faraday_functions);

        frequency = static_cast<const Kappa_transfer_f_arguments_type*>(p_Transfer_args)->frequency;

        break;

    case e_Phenomenological_ensamble:

        this->get_phenomenological_synchrotron_fit_functions(static_cast<const Phenomenological_transfer_f_arguments_type*>(p_Transfer_args), p_Transfer_functions);

        frequency = static_cast<const Phenomenological_transfer_f_arguments_type*>(p_Transfer_args)->frequency;

        break;

    default:

        std::cout << "\n" << "Error! Unsupported emission model - something broke in the evaluate_synchrotron_transfer_functions function!" << "\n";
        exit(ERROR);

    }

    const double f_cyclo = Q_ELECTRON_CGS * p_Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm / (2 * M_PI * M_ELECTRON_CGS * C_LIGHT_CGS);

    /* ================================================ The emission functions ================================================ */

    p_Transfer_functions->Emission_functions[I] *= p_Emission_medium_state->Density * f_cyclo * Q_ELECTRON_CGS * Q_ELECTRON_CGS / C_LIGHT_CGS;
    p_Transfer_functions->Emission_functions[Q] *= p_Emission_medium_state->Density * f_cyclo * Q_ELECTRON_CGS * Q_ELECTRON_CGS / C_LIGHT_CGS;
    p_Transfer_functions->Emission_functions[V] *= p_Emission_medium_state->Density * f_cyclo * Q_ELECTRON_CGS * Q_ELECTRON_CGS / C_LIGHT_CGS;

    /* ================================================ The absorbtion functions ================================================ */

    p_Transfer_functions->Absorbtion_functions[I] *= p_Emission_medium_state->Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / frequency / M_ELECTRON_CGS / C_LIGHT_CGS;
    p_Transfer_functions->Absorbtion_functions[Q] *= p_Emission_medium_state->Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / frequency / M_ELECTRON_CGS / C_LIGHT_CGS;
    p_Transfer_functions->Absorbtion_functions[V] *= p_Emission_medium_state->Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS / frequency / M_ELECTRON_CGS / C_LIGHT_CGS;

    /* ================================================ The Faraday functions ================================================ */
    /* Originally derived in https://iopscience.iop.org/article/10.1086/592326/pdf - expressions 25, 26 and 33. */

    p_Transfer_functions->Faraday_functions[Q] *= -p_Emission_medium_state->Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS * f_cyclo * f_cyclo / M_ELECTRON_CGS / C_LIGHT_CGS / frequency / frequency / frequency;
    p_Transfer_functions->Faraday_functions[V] *= 2 * p_Emission_medium_state->Density * Q_ELECTRON_CGS * Q_ELECTRON_CGS * f_cyclo / M_ELECTRON_CGS / C_LIGHT_CGS / frequency / frequency;

}

/* ========================================================== Misc Functions ========================================================== */

void Emission_models_class::precompute_electron_pitch_angles(Initial_conditions_type* p_Init_Conditions) {

    // ====================================================== Allocate memory for the arrays ====================================================== //

    this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];
    this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];

    this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];
    this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];

    this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_5175 = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];
    this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_6    = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];
    this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_7515 = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];

    this->s_Precomputed_e_pitch_angles.one_over_sin_to_7_over_20 = new double[p_Init_Conditions->Emission_pitch_angle_samples_to_average];

    // =========================================================================================================================================== //

    for (int index = 0; index < p_Init_Conditions->Emission_pitch_angle_samples_to_average; index++) {

        double pitch_angle = double(index) / p_Init_Conditions->Emission_pitch_angle_samples_to_average * M_PI;
        this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] = sin(pitch_angle);
        this->s_Precomputed_e_pitch_angles.cos_electron_pitch_angles[index] = cos(pitch_angle);

        if (this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index] != 0) {

            // Used in the thermal and kappa synchrotron emission functions.

            this->s_Precomputed_e_pitch_angles.one_over_sqrt_sin[index] = 1. / sqrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);
            this->s_Precomputed_e_pitch_angles.one_over_cbrt_sin[index] = 1. / cbrt(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index]);

            // Used in the thermal synchrotron Faraday functions.

            this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_5175[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 0.5175);
            this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_6[index]    = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 0.6);
            this->s_Precomputed_e_pitch_angles.one_over_sin_to_0_p_7515[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 0.7515);

            // Used in the kappa synchrotron emission functions.

            this->s_Precomputed_e_pitch_angles.one_over_sin_to_7_over_20[index] = 1. / pow(this->s_Precomputed_e_pitch_angles.sin_electron_pitch_angles[index], 7. / 20);

        }
    }
}

Emission_models_class::Emission_models_class(Simulation_Context_type* p_Sim_Context) {

    this->Num_samples_to_avg = p_Sim_Context->p_Init_Conditions->Emission_pitch_angle_samples_to_average;
    this->Include_polarization = p_Sim_Context->p_Init_Conditions->Observer_params.include_polarization;

    this->p_Sim_Context = p_Sim_Context;

    this->p_Disk_Model = new Disk_model_type(p_Sim_Context);
    this->p_Hotspot_Model = new Hotspot_model_type(p_Sim_Context);

    if (nullptr != p_Sim_Context) {

        this->s_Emission_params = p_Sim_Context->p_Init_Conditions->Emission_params;

    }
    else {

        throw std::runtime_error("Could not load the emission models parameter struct! \n");
        
    }

}