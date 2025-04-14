#include "Disk_Models.h"

Disk_model_type::Disk_model_type(Simulation_Context_type* p_Sim_Context) {

    if (NULL != p_Sim_Context) {

        this->s_Disk_params = p_Sim_Context->p_Init_Conditions->Disk_params;

    }
    else {

        std::cout << "Could not load the disk parameter struct! \n";
        exit(ERROR);

    }

}

double Disk_model_type::get_disk_profile(const Disk_profile_parameters_type* const p_Profile_parameters,
                                          Profile_enums e_Profile_type) const{

    double Profile = 0.0;

    switch (e_Profile_type) {

    case e_Power_law_profile:

        Profile = pow(p_Profile_parameters->r_0 / p_Profile_parameters->r, p_Profile_parameters->power) *
                  exp(-int_power(p_Profile_parameters->z / p_Profile_parameters->rho / p_Profile_parameters->tan_opening_angle, 2) / 2);

        if (p_Profile_parameters->r < p_Profile_parameters->r_cutoff) {

            Profile *= exp(-int_power((p_Profile_parameters->r - p_Profile_parameters->r_cutoff) / p_Profile_parameters->cutoff_scale, 2));

        }
        
        return Profile;

    case e_Exponential_law_profile:

        return exp(-int_power(p_Profile_parameters->r / p_Profile_parameters->exp_radial_scale, 2) / 2 - int_power(p_Profile_parameters->z / p_Profile_parameters->exp_height_scale, 2) / 2);


    default:

        std::cout << "Unsupported disk profile type! \n";
        exit(ERROR);

    }

}

void Disk_model_type::get_density_and_temperature(const double* const State_Vector,
                                                  Emission_medium_state_type* const p_Emission_medium_state) const {

    Disk_profile_parameters_type Profile_prameters{};

    /* ======================================== The density profile ======================================== */

    Profile_prameters.r   = State_Vector[e_r];
    Profile_prameters.z   = State_Vector[e_r] * cos(State_Vector[e_theta]);
    Profile_prameters.rho = State_Vector[e_r] * sin(State_Vector[e_theta]);

    /* The disk opening angle is a common parameter for both density and temperature profiles. */
    Profile_prameters.tan_opening_angle = this->s_Disk_params.Power_law_disk_opening_angle;

    Profile_prameters.r_0          = this->s_Disk_params.Power_law_density_R_0;
    Profile_prameters.r_cutoff     = this->s_Disk_params.Power_law_density_R_cutoff;
    Profile_prameters.cutoff_scale = this->s_Disk_params.Power_law_density_cutoff_scale;
    Profile_prameters.power        = this->s_Disk_params.Power_law_density_radial_power_law;

    Profile_prameters.exp_radial_scale = this->s_Disk_params.Exp_law_density_radial_scale;
    Profile_prameters.exp_height_scale = this->s_Disk_params.Exp_law_density_height_scale;

    p_Emission_medium_state->Density = this->s_Disk_params.Electron_density_scale * this->get_disk_profile(&Profile_prameters, this->s_Disk_params.Density_profile_type);

    if (isnan(p_Emission_medium_state->Density) || isinf(p_Emission_medium_state->Density) || p_Emission_medium_state->Density < 0) {

        std::cout << "Invalid disk density profile: " << p_Emission_medium_state->Density << "\n";

        exit(ERROR);

    }

    /* ====================================== The temperature profile ====================================== */

    Profile_prameters.r            = State_Vector[e_r];
    Profile_prameters.r_0          = this->s_Disk_params.Power_law_temperature_R_0;
    Profile_prameters.r_cutoff     = this->s_Disk_params.Power_law_temperature_R_cutoff;
    Profile_prameters.cutoff_scale = this->s_Disk_params.Power_law_temperature_cutoff_scale;
    Profile_prameters.power        = this->s_Disk_params.Power_law_temperature_radial_power_law;

    /* NOTE: The temperature is modelled to only decrease radially, so I set the opening angle to pi / 2 to get rid of the vertical dependance of the profile.
       This should be offloaded to the configurator and not hardcoded here! */
    Profile_prameters.tan_opening_angle = std::numeric_limits<double>::infinity();

    p_Emission_medium_state->Temperature = this->s_Disk_params.Electron_temperature_scale * this->get_disk_profile(&Profile_prameters, this->s_Disk_params.Temperature_profile_type);

    if (isnan(p_Emission_medium_state->Temperature) || isinf(p_Emission_medium_state->Temperature) || p_Emission_medium_state->Temperature < 0) {

        std::cout << "Invalid disk temperature profile: " << p_Emission_medium_state->Temperature << "\n";

        exit(ERROR);

    }

}
