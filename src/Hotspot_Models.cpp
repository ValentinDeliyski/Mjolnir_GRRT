#include "Hotspot_Models.h"

Hotspot_model_type::Hotspot_model_type(Simulation_Context_type* p_Sim_Context) {

    if (nullptr != p_Sim_Context) {

        this->s_Hotspot_params = p_Sim_Context->p_Init_Conditions->Hotspot_params;

    }
    else { throw std::runtime_error("Could not load the hotspot parameter struct! \n"); }

}

Hotspot_position_type Hotspot_model_type::get_hotspot_position(const double* const State_Vector,
                                                               const double* const Hotspot_Velocity) const {

    Hotspot_position_type Hotspot_position{};
    double Hotspot_ang_velocity{};

    if (nullptr != Hotspot_Velocity) { Hotspot_ang_velocity = Hotspot_Velocity[e_phi] / Hotspot_Velocity[e_t]; }

    Hotspot_position.Distance    = this->s_Hotspot_params.Position[e_r];
    Hotspot_position.Inclination = M_PI_2;
    Hotspot_position.Azimuth     = this->s_Hotspot_params.Position[e_phi] + Hotspot_ang_velocity * State_Vector[e_t];

   double sin_hotspot_inclination = sin(Hotspot_position.Inclination);

    Hotspot_position.x = Hotspot_position.Distance * sin_hotspot_inclination * cos(Hotspot_position.Azimuth);
    Hotspot_position.y = Hotspot_position.Distance * sin_hotspot_inclination * sin(Hotspot_position.Azimuth);
    Hotspot_position.z = Hotspot_position.Distance * cos(Hotspot_position.Inclination);

    return Hotspot_position;

}

double Hotspot_model_type::get_hotspot_profile(const Hotspot_profile_parameters_type* const p_Profile_parameters,
                                               Profile_enums e_Profile_type) const {

    double exponent_term{};

    switch (e_Profile_type) {

    case e_Gaussian:

        return exp(-pow((p_Profile_parameters->Gaussian_variable - p_Profile_parameters->Gaussian_mean) / p_Profile_parameters->Gaussian_spread, 2) / 2);

    case e_Spherical:

        if (p_Profile_parameters->Distance_from_sphere_center < p_Profile_parameters->Sphere_radius) {
            return 1.0; 
        }
        else { return 0.0; }

    case e_Hybrid_power_gaussian:

        exponent_term = exp(-pow((p_Profile_parameters->Gaussian_variable - p_Profile_parameters->Gaussian_mean) / p_Profile_parameters->Gaussian_spread, 2) / 2);

        return pow(p_Profile_parameters->Power_law_scale / p_Profile_parameters->Power_law_variable, p_Profile_parameters->Power_law_power) * exponent_term;

    default: throw std::runtime_error("Unsupported hotspot profile type! \n");

    }

};

void Hotspot_model_type::get_density_and_temperature(const double* const State_Vector,
                                                     Emission_medium_state_type* const p_Emission_medium_state) const {

    Hotspot_position_type Hotspot_position = this->get_hotspot_position(State_Vector, p_Emission_medium_state->Plasma_Velocity);

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

    if (Distance_to_hotspot_center < 3.) {

        int test{};

    }

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
    
    if (isnan(p_Emission_medium_state->Density) || isinf(p_Emission_medium_state->Density) || p_Emission_medium_state->Density < 0) {

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

    if (isnan(p_Emission_medium_state->Temperature) || isinf(p_Emission_medium_state->Temperature) || p_Emission_medium_state->Temperature < 0) {

        throw std::runtime_error(std::format("Invalid hotspot temperature profile: {} \n", p_Emission_medium_state->Temperature));

    }

}

bool Hotspot_model_type::is_inside_hotspot(const double* const State_Vector, Emission_medium_state_type* const Hotspot_State) const {

    this->get_density_and_temperature(State_Vector, Hotspot_State);

    return (Hotspot_State->Density / this->s_Hotspot_params.Electron_density_scale > this->s_Hotspot_params.Threshold_relative_density);

}