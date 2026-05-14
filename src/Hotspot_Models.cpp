#include "Hotspot_Models.h"

Hotspot_model_type::Hotspot_model_type(Simulation_Context_type* p_Sim_Context) {

    if (nullptr != p_Sim_Context) {

        this->s_Hotspot_params = p_Sim_Context->p_Init_Conditions->Hotspot_params;
        this->p_Spacetime = p_Sim_Context->p_Spacetime;
    }
    else { throw std::runtime_error("Could not load the hotspot parameter struct! \n"); }

    memset(this->Current_Velocity, 0, 4 * sizeof(double));

    this->Current_Position.Distance = this->s_Hotspot_params.Init_Position[e_r];
    this->Current_Position.Inclination = this->s_Hotspot_params.Init_Position[e_theta];
    this->Current_Position.Azimuth = this->s_Hotspot_params.Init_Position[e_phi];

}

double* Hotspot_model_type::get_hotspot_velocity(bool Eval_at_hotspot_center, const double* const Local_State_Vector) {

    double Hotspot_State_Vector[e_Dynamic_state_size]{};

    Hotspot_State_Vector[e_r] = this->Current_Position.Distance;
    Hotspot_State_Vector[e_theta] = this->Current_Position.Inclination;
    Hotspot_State_Vector[e_phi] = this->Current_Position.Azimuth;

    Metric_type s_Metric{};

    if (Eval_at_hotspot_center) {

        s_Metric = this->p_Spacetime->get_local_metric(Hotspot_State_Vector);

    }
    else {

        s_Metric = this->p_Spacetime->get_local_metric(Local_State_Vector);

    }

    Metric_type s_dr_Metric = this->p_Spacetime->get_dr_local_metric(Hotspot_State_Vector);

    /* = References for the sake of readability = */

    const auto& g = s_Metric.Metric;
    const auto& dr_g = s_dr_Metric.Metric;

    /* ========================================== */

    if (e_Circular_fixed_rate != this->s_Hotspot_params.Velocity_profile_type) {

        throw std::runtime_error("Invalid hotspot velocity profile!");

    }

    /* This is the angular velocity of a geodesic with zero a radial velocity component. */
    const double Omega = (-dr_g[e_t][e_phi] + sqrt(dr_g[e_t][e_phi] * dr_g[e_t][e_phi] - dr_g[e_t][e_t] * dr_g[e_phi][e_phi])) / dr_g[e_phi][e_phi];
    const double Normalization = g[e_t][e_t] + 2 * g[e_t][e_phi] * Omega + g[e_phi][e_phi] * Omega * Omega;

    this->Current_Velocity[e_t]     = sqrt(-1.0 / Normalization);
    this->Current_Velocity[e_r]     = 0.0;
    this->Current_Velocity[e_theta] = 0.0;
    this->Current_Velocity[e_phi]   = this->Current_Velocity[e_t] * Omega;

    if (isnan(this->Current_Velocity[e_t]) or
        isinf(this->Current_Velocity[e_t]) or
        isnan(this->Current_Velocity[e_phi]) or
        isinf(this->Current_Velocity[e_phi])) {

        throw std::runtime_error("Invalid hotspot velocity!");

    }

    if (nullptr == this->Current_Velocity) { throw std::runtime_error("Hotspot velocity is a null pointer!"); }

    return this->Current_Velocity;

}

Hotspot_position_type Hotspot_model_type::get_hotspot_position(const double* const Local_State_Vector) {

    /* The entire hotspot is supposed to move with the angular velocity of its center, 
       so here I just get the velocity of the center. The call to get_hotspot_velocity() updates
       the velocity internally in the this->Current_Velocity variable. */

    this->get_hotspot_velocity(true, Local_State_Vector);

    const double Hotspot_ang_velocity = this->Current_Velocity[e_phi] / this->Current_Velocity[e_t];

    this->Current_Position.Distance    = this->s_Hotspot_params.Init_Position[e_r];
    this->Current_Position.Inclination = this->s_Hotspot_params.Init_Position[e_theta];
    this->Current_Position.Azimuth     = this->s_Hotspot_params.Init_Position[e_phi] + Hotspot_ang_velocity * Local_State_Vector[e_t];

   double sin_hotspot_inclination = sin(this->Current_Position.Inclination);

   this->Current_Position.x = this->Current_Position.Distance * sin_hotspot_inclination * cos(this->Current_Position.Azimuth);
   this->Current_Position.y = this->Current_Position.Distance * sin_hotspot_inclination * sin(this->Current_Position.Azimuth);
   this->Current_Position.z = this->Current_Position.Distance * cos(this->Current_Position.Inclination);

    return this->Current_Position;

}

double Hotspot_model_type::get_hotspot_profile(const Hotspot_profile_parameters_type* const p_Profile_parameters,
                                               Profile_enums e_Profile_type) const {

    double exponent_term{};

    switch (e_Profile_type) {

    case e_Gaussian:

        return exp(-pow((p_Profile_parameters->Gaussian_variable - p_Profile_parameters->Gaussian_mean) / p_Profile_parameters->Gaussian_spread, 2) / 2);

    case e_Spherical:

        if (p_Profile_parameters->Distance_from_sphere_center < p_Profile_parameters->Sphere_radius) { return 1.0; }
        else { return 0.0; }

    case e_Hybrid_power_gaussian:

        exponent_term = exp(-pow((p_Profile_parameters->Gaussian_variable - p_Profile_parameters->Gaussian_mean) / p_Profile_parameters->Gaussian_spread, 2) / 2);

        return pow(p_Profile_parameters->Power_law_scale / p_Profile_parameters->Power_law_variable, p_Profile_parameters->Power_law_power) * exponent_term;

    default: throw std::runtime_error("Unsupported hotspot profile type! \n");

    }

};

void Hotspot_model_type::get_density_and_temperature(const double* const State_Vector,
                                                     Emission_medium_state_type* const p_Emission_medium_state) {

    Hotspot_position_type Hotspot_position = this->get_hotspot_position(State_Vector);

    const double& photon_r  = State_Vector[e_r];
    double sin_photon_theta = sin(State_Vector[e_theta]);

    double x_photon = photon_r * sin_photon_theta * cos(State_Vector[e_phi]);
    double y_photon = photon_r * sin_photon_theta * sin(State_Vector[e_phi]);
    double z_photon = photon_r * cos(State_Vector[e_theta]);

    double Distance_to_hotspot_center = sqrt((Hotspot_position.x - x_photon) * (Hotspot_position.x - x_photon)
                                           + (Hotspot_position.y - y_photon) * (Hotspot_position.y - y_photon)
                                           + (Hotspot_position.z - z_photon) * (Hotspot_position.z - z_photon));

    Hotspot_profile_parameters_type Profile_prameters{};

    /* ======================= The spatial part of the density profile ======================= */

    Profile_prameters.Gaussian_variable = Distance_to_hotspot_center;

    Profile_prameters.Gaussian_spread   = this->s_Hotspot_params.Profile_params.Density_gaussian_spread;
    Profile_prameters.Gaussian_mean     = 0.0;

    Profile_prameters.Sphere_radius = this->s_Hotspot_params.Profile_params.Radius;
    Profile_prameters.Distance_from_sphere_center = Distance_to_hotspot_center;

    Profile_prameters.Power_law_variable = State_Vector[e_r];
    Profile_prameters.Power_law_power = this->s_Hotspot_params.Profile_params.Density_power_law_power;
    Profile_prameters.Power_law_scale = this->s_Hotspot_params.Profile_params.Density_power_law_scale;

    double Spatial_profile = this->get_hotspot_profile(&Profile_prameters, this->s_Hotspot_params.Density_profile_type);

    /* ======================= The temporal part of the density profile ======================= */

    Profile_prameters.Gaussian_variable = State_Vector[e_t];
    Profile_prameters.Gaussian_spread = this->s_Hotspot_params.Profile_params.Temporal_gaussian_spread;
    Profile_prameters.Gaussian_mean   = -this->s_Hotspot_params.Profile_params.Coord_time_offset;

    double Temporal_profile = this->get_hotspot_profile(&Profile_prameters, e_Gaussian);

    /* ======================= =============================== ======================= */

    p_Emission_medium_state->Density = this->s_Hotspot_params.Electron_density_scale * Spatial_profile * Temporal_profile;
    
    if (isnan(p_Emission_medium_state->Density) or isinf(p_Emission_medium_state->Density) or p_Emission_medium_state->Density < 0.0) {

        throw std::runtime_error(std::format("Invalid hotspot density profile: {} \n", p_Emission_medium_state->Density));

    }

    /* ======================= The spatial part of the temperature profile ======================= */

    Profile_prameters.Gaussian_variable = Distance_to_hotspot_center;
    Profile_prameters.Gaussian_spread   = this->s_Hotspot_params.Profile_params.Temperature_gaussian_spread;
    Profile_prameters.Gaussian_mean     = 0.0;

    Profile_prameters.Power_law_variable = State_Vector[e_r];
    Profile_prameters.Power_law_power    = this->s_Hotspot_params.Profile_params.Temperature_power_law_power;
    Profile_prameters.Power_law_scale    = this->s_Hotspot_params.Profile_params.Temperature_power_law_scale;

    Spatial_profile = this->get_hotspot_profile(&Profile_prameters, this->s_Hotspot_params.Temperature_profile_type);

    /* ======================= The temporal part of the temperaure profile ======================= */

    /* Currently the temporal profile is the same for both the density and the temperature, 
       so I don't update the param struct. */

    /* ======================= =========================================== ======================= */

    p_Emission_medium_state->Temperature = this->s_Hotspot_params.Electron_temperature_scale * Spatial_profile * Temporal_profile;

    if (isnan(p_Emission_medium_state->Temperature) or isinf(p_Emission_medium_state->Temperature) or p_Emission_medium_state->Temperature < 0) {

        throw std::runtime_error(std::format("Invalid hotspot temperature profile: {} \n", p_Emission_medium_state->Temperature));

    }

}

bool Hotspot_model_type::is_inside_hotspot(const double* const State_Vector, Emission_medium_state_type* const Hotspot_State) {

    if (!this->s_Hotspot_params.Enable_flag) {

        return false;

    }

    this->get_density_and_temperature(State_Vector, Hotspot_State);

    return (Hotspot_State->Density / this->s_Hotspot_params.Electron_density_scale > this->s_Hotspot_params.Threshold_relative_density);

}

void Hotspot_model_type::get_magnetic_field(const double* const Local_State_Vector,
                                             const Metric_type* const p_Metric,
                                             Emission_medium_state_type* const Emission_medium_state) const {

    /* ================================================================================================================================================================ *|
    |                                                                                                                                                                    |
    |  The reference for this implementation is https://arxiv.org/pdf/2404.13824v1, expressions (1.54). The desired megnetic field geometry is specified for an          |
    |  Eualrian observer, with covarian 4-velocity n_mu = (-Lapse, 0, 0, 0). Writing the dual Maxwell tensor in terms of the magnetic 4-vector measured by a comoving    |
    |  with  the plasma observer, and his 4-velocity (1.20 - but they have a overall missing minus sign for some reason), one can express the Eulerian magnetic field by |
    |  projecting the *F^mu^nu onto n_mu. Inverting this expression, one obtains the magnetic 4-vector measured by the comoving observer in terms of the one measured by |
    |  the Eulerian observer.                                                                                                                                            |
    |                                                                                                                                                                    |
    * ================================================================================================================================================================= */

    /* ======================= References for the sake of readability ======================= */

    double (&B_eulerian)[4] = Emission_medium_state->Magnetic_fields.B_field_eulerian_frame;
    double (&B_plasma)[4] = Emission_medium_state->Magnetic_fields.B_field_plasma_frame;
    double& B_plasma_norm = Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm;

    const double& r = Local_State_Vector[e_r];

    const double& Power = this->s_Hotspot_params.Mag_field_power;
    const double& B_0 = this->s_Hotspot_params.Mag_field_B_0;
    const double& r_0 = this->s_Hotspot_params.Mag_field_r_0;

    /* ====================================================================================== */

    B_eulerian[e_t] = 0.0;

    switch (this->s_Hotspot_params.e_Mag_field_geometry) {

    case Constant:

        B_eulerian[e_r]     = Emission_medium_state->Magnetic_fields.Mag_field_geometry_vector[e_r - 1];
        B_eulerian[e_theta] = Emission_medium_state->Magnetic_fields.Mag_field_geometry_vector[e_theta - 1];
        B_eulerian[e_phi]   = Emission_medium_state->Magnetic_fields.Mag_field_geometry_vector[e_phi - 1];

        break;

    case Vertical:

        B_eulerian[e_r]     = cos(Local_State_Vector[e_theta]);
        B_eulerian[e_theta] = -sin(Local_State_Vector[e_theta]);
        B_eulerian[e_phi]   = 0;

        break;

    case Toroidal:

        B_eulerian[e_r]     = 0;
        B_eulerian[e_theta] = 0;
        B_eulerian[e_phi]   = 1;

        break;

    default:

        throw std::runtime_error("Unsupported magnetic field geometry!");

    }

    /* --------------------------------------------- Normalize the magnetic vector in the Eularian frame. --------------------------------------------- */

    double Mag_field_eularian_norm{};

    for (int idx = 0; idx < 4; idx++) {

        Mag_field_eularian_norm += B_eulerian[idx] * B_eulerian[idx];

    }

    for (int idx = 0; idx < 4; idx++) {

        B_eulerian[idx] /= sqrt(Mag_field_eularian_norm);

    }

    /* ------------------------------------------------------------------------------------------------------------------------------------------------ */

    const double Lorentz_factor = this->Current_Velocity[e_t] * p_Metric->Lapse_function;

    for (int left_idx = 0; left_idx < 4; left_idx++) {

        for (int right_idx = 0; right_idx < 4; right_idx++) {

            B_plasma[e_t] += p_Metric->Metric[left_idx][right_idx] * this->Current_Velocity[left_idx] * B_eulerian[right_idx] / p_Metric->Lapse_function;
        }

    }

    for (int index = 0; index < 4; index++) {

        B_plasma[index] = (B_eulerian[index] + p_Metric->Lapse_function * B_plasma[e_t] * this->Current_Velocity[index]) / Lorentz_factor;

    }

    switch (this->s_Hotspot_params.e_Mag_field_magnitude_profile) {

    case Magnetization_based:

        B_plasma_norm = sqrt(this->s_Hotspot_params.Magnetization * C_LIGHT_CGS * C_LIGHT_CGS * Emission_medium_state->Density * M_PROTON_CGS * 4 * std::numbers::pi);
        break;

    case Power_law_based:

        B_plasma_norm = B_0 * pow(r_0 / r, Power);
        break;

    default:

        throw std::runtime_error("Unsupported magnetic field magnitude profile for the hotspot!");

    }

}