#include "Disk_Models.h"

Disk_model_type::Disk_model_type(Simulation_Context_type* p_Sim_Context) {

    if (nullptr != p_Sim_Context) {

        this->s_Disk_params = p_Sim_Context->p_Init_Conditions->Disk_params;
        this->p_Sim_Context = p_Sim_Context;

    }
    else { throw std::runtime_error("Could not load the disk parameter struct! \n"); }

    if (this->s_Disk_params.e_Disk_model == Disk_model_enums::e_Numerical_Polytrope) {

        if (this->s_Disk_params.Numerical_disk_params.e_Spline_type == Spline_selection_enums::Custom_cubic) {

            throw std::runtime_error("Custom interpolants are not supported for numerical disks. Use GSL splines. \n"); 
    
        }

        if (this->s_Disk_params.Numerical_disk_params.e_Spline_type == Spline_selection_enums::GSL_linear) {

            this->Spline_instance_density = gsl_spline2d_alloc(gsl_interp2d_bilinear, this->s_Disk_params.Numerical_disk_params.Radial_grid_size, this->s_Disk_params.Numerical_disk_params.Theta_grid_size);

        }
        else {

            this->Spline_instance_density = gsl_spline2d_alloc(gsl_interp2d_bicubic, this->s_Disk_params.Numerical_disk_params.Radial_grid_size, this->s_Disk_params.Numerical_disk_params.Theta_grid_size);

        }

        this->Radial_interp_accelerator = gsl_interp_accel_alloc();
        this->Theta_interp_accelerator = gsl_interp_accel_alloc();

        gsl_spline2d_init(this->Spline_instance_density,
                          this->s_Disk_params.Numerical_disk_params.Radial_grid,
                          this->s_Disk_params.Numerical_disk_params.Theta_grid,
                          this->s_Disk_params.Numerical_disk_params.Raw_density_data,
                          this->s_Disk_params.Numerical_disk_params.Radial_grid_size,
                          this->s_Disk_params.Numerical_disk_params.Theta_grid_size);

    }

}

Disk_model_type::~Disk_model_type() {

    if (this->s_Disk_params.e_Disk_model == Disk_model_enums::e_Numerical_Polytrope) {

        gsl_spline2d_free(this->Spline_instance_density);

        gsl_interp_accel_free(this->Radial_interp_accelerator);
        gsl_interp_accel_free(this->Theta_interp_accelerator);

        delete this->s_Disk_params.Numerical_disk_params.Raw_density_data;

        delete this->s_Disk_params.Numerical_disk_params.Radial_grid;
        delete this->s_Disk_params.Numerical_disk_params.Theta_grid;

    }
}

double Disk_model_type::get_disk_internal_energy(double density) const {

    /* This assumes a ideal fluid, which is undergoing an iso-entropic process (Rezzolla (2.248)).
       The polytropic index of the polytropic EOS (Gamma) is assumed to be equal to the adiabatic index,
       which appears in the ideal fluid thermal EOS (2.228) */

    const double& K = this->s_Disk_params.Thermal_EOS_params.Polytrope_Coeff;
    const double& Gamma = this->s_Disk_params.Thermal_EOS_params.Polytrope_Power;

    return K / (Gamma - 1) * std::pow(density, Gamma - 1);
}

double Disk_model_type::get_disk_profile(const Disk_profile_parameters_type* const p_Profile_parameters,
                                         Profile_enums e_Profile_type) const {

    double Profile{};
    double Exponent_arg{};

    switch (e_Profile_type) {

    case e_Power_law:

        Profile = std::pow(p_Profile_parameters->power_law_scale / p_Profile_parameters->radial_coordinate, p_Profile_parameters->power);
        break;

    case e_Hybrid_power_gaussian:

        Exponent_arg = (p_Profile_parameters->gaussian_variable - p_Profile_parameters->gaussian_mean) / p_Profile_parameters->gaussian_std;

        Profile = std::pow(p_Profile_parameters->power_law_scale / p_Profile_parameters->radial_coordinate, p_Profile_parameters->power) * exp(-std::pow(Exponent_arg, 2) / 2);
        break;

    case e_Gaussian:

        Exponent_arg = (p_Profile_parameters->gaussian_variable - p_Profile_parameters->gaussian_mean) / p_Profile_parameters->gaussian_std;

        Profile = std::exp(-std::pow(Exponent_arg, 2) / 2);
        break;

    default: throw std::runtime_error("Unsupported disk profile type! \n");

    }

    if (p_Profile_parameters->radial_coordinate < p_Profile_parameters->cutoff_radius) {

        double Cutoff_exponent_arg = (p_Profile_parameters->radial_coordinate - p_Profile_parameters->cutoff_radius) / p_Profile_parameters->cutoff_scale;

        Profile *= exp(-std::pow(Cutoff_exponent_arg, 2) / 2);

    }

    return Profile;

}



void Disk_model_type::get_density_and_temperature(const double* const State_Vector,
                                                  Emission_medium_state_type* const p_Emission_medium_state) const {

    Disk_profile_parameters_type Density_profile_params{}, Temperature_profile_params{};

    const double& Gamma = this->s_Disk_params.Thermal_EOS_params.Polytrope_Power;
    int Err_code = 0;

    switch (this->s_Disk_params.e_Disk_model) {

    case e_Numerical_Polytrope:

        /* ------------------------------------------------ Get the density profile ------------------------------------------------ */

        Err_code = gsl_spline2d_eval_e(this->Spline_instance_density,
                                       State_Vector[e_r],
                                       State_Vector[e_theta],
                                       this->Radial_interp_accelerator,
                                       this->Theta_interp_accelerator,
                                       &p_Emission_medium_state->Density);

        if (GSL_EDOM == Err_code) { p_Emission_medium_state->Density = 0.0; }

        /* TODO: Scale this thing so it comes out in g / cm^3 */

        /* ------------------------------------------------ Get the temperature profile ------------------------------------------------ */

        p_Emission_medium_state->Temperature = M_PROTON_SI / BOLTZMANN_CONST_SI * (Gamma - 1) * this->get_disk_internal_energy(p_Emission_medium_state->Density);

        /* TODO: Scale this thing so it comes out in K */

        break;

    case e_Phenom_RIAF_1:

        /* ============= This is the model used in https://arxiv.org/pdf/2206.12066, with an added cutoff exponential. ============= */

        /* ------------------------------------------------ Get the density profile ------------------------------------------------ */

        Density_profile_params.radial_coordinate = State_Vector[e_r];
        Density_profile_params.power_law_scale   = this->s_Disk_params.Common_RIAF_params.Density_power_law_scale;
        Density_profile_params.power             = this->s_Disk_params.Common_RIAF_params.Density_power_law_power;

        Density_profile_params.gaussian_variable = cos(State_Vector[e_theta]) / sin(State_Vector[e_theta]);
        Density_profile_params.gaussian_mean     = 0.0;
        Density_profile_params.gaussian_std      = this->s_Disk_params.Common_RIAF_params.Disk_opening_angle;

        Density_profile_params.cutoff_radius = this->s_Disk_params.Common_RIAF_params.Density_cutoff_radius;
        Density_profile_params.cutoff_scale = this->s_Disk_params.Common_RIAF_params.Density_cutoff_scale;

        p_Emission_medium_state->Density = this->s_Disk_params.Electron_density_scale * this->get_disk_profile(&Density_profile_params, e_Hybrid_power_gaussian);

        /* ---------------------------------------------- Get the temperature profile ---------------------------------------------- */

        Temperature_profile_params.radial_coordinate = State_Vector[e_r];
        Temperature_profile_params.power_law_scale   = this->s_Disk_params.Common_RIAF_params.Temperature_power_law_scale;
        Temperature_profile_params.power             = this->s_Disk_params.Common_RIAF_params.Temperature_power_law_power;

        /* Vertical gaussians are not used for the temperature profile. */

        Temperature_profile_params.gaussian_variable = 0.0;
        Temperature_profile_params.gaussian_mean     = 0.0;
        Temperature_profile_params.gaussian_std      = 0.0;

        Temperature_profile_params.cutoff_radius = this->s_Disk_params.Common_RIAF_params.Temperature_cutoff_radius;
        Temperature_profile_params.cutoff_scale = this->s_Disk_params.Common_RIAF_params.Temperature_cutoff_scale;

        p_Emission_medium_state->Temperature = this->s_Disk_params.Electron_temperature_scale * this->get_disk_profile(&Temperature_profile_params, e_Power_law);

        break;

    case e_Phenom_RIAF_2:

        /* ============= This is the model used in https://arxiv.org/pdf/2209.09931, with an added cutoff exponential. ============= */

        /* ------------------------------------------------ Get the density profile ------------------------------------------------ */

        Density_profile_params.radial_coordinate = State_Vector[e_r];
        Density_profile_params.power_law_scale   = this->s_Disk_params.Common_RIAF_params.Density_power_law_scale;
        Density_profile_params.power             = this->s_Disk_params.Common_RIAF_params.Density_power_law_power;

        Density_profile_params.gaussian_variable = cos(State_Vector[e_theta]);
        Density_profile_params.gaussian_mean     = 0.0;
        Density_profile_params.gaussian_std      = this->s_Disk_params.Common_RIAF_params.Disk_opening_angle;

        Density_profile_params.cutoff_radius = this->s_Disk_params.Common_RIAF_params.Density_cutoff_radius;
        Density_profile_params.cutoff_scale = this->s_Disk_params.Common_RIAF_params.Density_cutoff_scale;

        p_Emission_medium_state->Density = this->s_Disk_params.Electron_density_scale * this->get_disk_profile(&Density_profile_params, e_Hybrid_power_gaussian);

        /* ---------------------------------------------- Get the temperature profile ---------------------------------------------- */

        Temperature_profile_params.radial_coordinate = State_Vector[e_r];
        Temperature_profile_params.power_law_scale   = this->s_Disk_params.Common_RIAF_params.Temperature_power_law_scale;
        Temperature_profile_params.power             = this->s_Disk_params.Common_RIAF_params.Temperature_power_law_power;

        /* Vertical gaussians are not used for the temperature profile. */

        Temperature_profile_params.gaussian_variable = 0.0;
        Temperature_profile_params.gaussian_mean     = 0.0;
        Temperature_profile_params.gaussian_std      = 0.0;

        Temperature_profile_params.cutoff_radius = this->s_Disk_params.Common_RIAF_params.Temperature_cutoff_radius;
        Temperature_profile_params.cutoff_scale = this->s_Disk_params.Common_RIAF_params.Temperature_cutoff_scale;

        p_Emission_medium_state->Temperature = this->s_Disk_params.Electron_temperature_scale * this->get_disk_profile(&Temperature_profile_params, e_Power_law);
        
        break;

    case e_Phenom_RIAF_3:

        /* ============= This is anlagous to https://iopscience.iop.org/article/10.3847/1538-4357/ab96c6, but I offset the radial vairable. ============= */

        Density_profile_params.radial_coordinate = State_Vector[e_r];

        /* ------------------------------------------------ Get the radial density profile ------------------------------------------------ */

        Density_profile_params.gaussian_variable = State_Vector[e_r];
        Density_profile_params.gaussian_mean = this->s_Disk_params.Common_RIAF_params.Density_cutoff_radius;
        Density_profile_params.gaussian_std  = this->s_Disk_params.Common_RIAF_params.Density_cutoff_scale;

        /* -------- This profile is essentially _just_ a cutoff exponent, so I don't want to add on antother one ontop of that -> set the curoff radius to zero,
                    so the check for adding it on never passes -------- */

        Density_profile_params.cutoff_radius = 0.0;
        Density_profile_params.cutoff_scale  = 0.0;

        p_Emission_medium_state->Density = this->s_Disk_params.Electron_density_scale * this->get_disk_profile(&Density_profile_params, e_Gaussian);

        /* ------------------------------------------------- Get the theta density profile ------------------------------------------------- */
        
        Density_profile_params.gaussian_variable = cos(State_Vector[e_theta]);
        Density_profile_params.gaussian_mean = 0.0;
        Density_profile_params.gaussian_std  = this->s_Disk_params.Common_RIAF_params.Disk_opening_angle;

        /* -------- This profile is essentially _just_ a cutoff exponent, so I don't want to add on antother one ontop of that -> set the curoff radius to zero,
                    so the check for adding it on never passes -------- */

        Density_profile_params.cutoff_radius = 0.0;
        Density_profile_params.cutoff_scale  = 0.0;

        p_Emission_medium_state->Density *= this->get_disk_profile(&Density_profile_params, e_Gaussian);

        /* ------------------------------------------------- Get the temperature profile ------------------------------------------------- */

        Temperature_profile_params.radial_coordinate = State_Vector[e_r];

        Temperature_profile_params.gaussian_variable = State_Vector[e_r];
        Temperature_profile_params.gaussian_mean = this->s_Disk_params.Common_RIAF_params.Temperature_cutoff_radius;
        Temperature_profile_params.gaussian_std = this->s_Disk_params.Common_RIAF_params.Temperature_cutoff_scale;

        /* -------- This profile is essentially _just_ a cutoff exponent, so I don't want to add on antother one ontop of that -> set the curoff radius to zero,
                    so the check for adding it on never passes -------- */

        Temperature_profile_params.cutoff_radius = 0.0;
        Temperature_profile_params.cutoff_scale = 0.0;

        p_Emission_medium_state->Temperature = this->s_Disk_params.Electron_temperature_scale * this->get_disk_profile(&Temperature_profile_params, e_Gaussian);

        break;

    case e_Colab_test_1:

        /* =============== This is the model used in https://iopscience.iop.org/article/10.3847/1538-4357/ab96c6/pdf ============== */

        /* ------------------------------------------------ Get the density profile ------------------------------------------------ */

        Density_profile_params.radial_coordinate = 0.0; // This is not used in this profile, so I set it to zero.
        Density_profile_params.gaussian_variable = State_Vector[e_r];
        Density_profile_params.gaussian_mean     = 0.0;
        Density_profile_params.gaussian_std      = this->s_Disk_params.Colab_test_1_params.Radial_scale;

        p_Emission_medium_state->Density = this->s_Disk_params.Electron_density_scale * this->get_disk_profile(&Density_profile_params, e_Gaussian);

        /* The above only evaluates the radial part of the profile. Below we evaluate the vertical part (another gaussian). */

        Density_profile_params.radial_coordinate = 0.0; // This is not used in this profile, so I set it to zero.
        Density_profile_params.gaussian_variable = this->s_Disk_params.Colab_test_1_params.Vertical_scale * cos(State_Vector[e_theta]);
        Density_profile_params.gaussian_mean     = 0.0;
        Density_profile_params.gaussian_std      = 1.0;

        /* Note that the two profiles multiply together. */
        p_Emission_medium_state->Density *= this->get_disk_profile(&Density_profile_params, e_Gaussian);

        /* ---------------------------------------------- Get the temperature profile ---------------------------------------------- */
        /* This model does not specify a temperature profile at all. */

        p_Emission_medium_state->Temperature = 0.0;

        break;

    case e_Debug_constant_density:

        p_Emission_medium_state->Density = this->s_Disk_params.Electron_density_scale;
        p_Emission_medium_state->Temperature = this->s_Disk_params.Electron_temperature_scale;

        break;
        
    default: throw std::runtime_error("Unsupported disk profile type! \n");

    }

    if (isnan(p_Emission_medium_state->Density) or isinf(p_Emission_medium_state->Density) or p_Emission_medium_state->Density < 0) {

        throw std::runtime_error(std::format("Invalid disk density profile: {} \n", p_Emission_medium_state->Density));

    }

    if (isnan(p_Emission_medium_state->Temperature) or isinf(p_Emission_medium_state->Temperature) or p_Emission_medium_state->Temperature < 0) {

        throw std::runtime_error(std::format("Invalid disk temperature profile: {} \n", p_Emission_medium_state->Temperature));

    }

}

bool Disk_model_type::is_inside_disk(const double* const State_Vector, Emission_medium_state_type* const Disk_State) const {

    this->get_density_and_temperature(State_Vector, Disk_State);

    return (Disk_State->Density / this->s_Disk_params.Electron_density_scale > this->s_Disk_params.Threshold_relative_density);

}


void Disk_model_type::get_magnetic_field(const double* const Local_State_Vector,
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

    const double& Power = this->s_Disk_params.Mag_field_power;
    const double& B_0 = this->s_Disk_params.Mag_field_B_0;
    const double& r_0 = this->s_Disk_params.Mag_field_r_0;

    /* ====================================================================================== */

    B_eulerian[e_t] = 0.0;

    switch (this->s_Disk_params.e_Mag_field_geometry) {

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

    const double Lorentz_factor = Emission_medium_state->Plasma_Velocity[e_t] * p_Metric->Lapse_function;

    for (int left_idx = 0; left_idx < 4; left_idx++) {

        for (int right_idx = 0; right_idx < 4; right_idx++) {

            B_plasma[e_t] += p_Metric->Metric[left_idx][right_idx] * Emission_medium_state->Plasma_Velocity[left_idx] * B_eulerian[right_idx] / p_Metric->Lapse_function;
        }

    }

    for (int index = 0; index < 4; index++) {

        B_plasma[index] = (B_eulerian[index] + p_Metric->Lapse_function * B_plasma[e_t] * Emission_medium_state->Plasma_Velocity[index]) / Lorentz_factor;

    }

    switch (this->s_Disk_params.e_Mag_field_magnitude_profile) {

    case Magnetization_based:

        B_plasma_norm = sqrt(this->s_Disk_params.Magnetization * C_LIGHT_CGS * C_LIGHT_CGS * Emission_medium_state->Density * M_PROTON_CGS * 4 * std::numbers::pi);
        break;

    case Power_law_based:

        B_plasma_norm = B_0 * pow(r_0 / r, Power);
        break;

    default:
        std::cout << "Unsupported magnetic field magnitude profile! \n";
        exit(ERROR);

    }

}

const double* const Disk_model_type::get_disk_velocity(const double* const Local_State_Vector) {

    /* The reference for this implementation is https://arxiv.org/pdf/2206.12066. */

    /* === Initialize some variables === */
    double Omega{}, rho{}, ell{}, u_t{}, u_phi{}, Normalization{}, inv_metric[4][4]{};


    Metric_type s_Metric = this->p_Sim_Context->p_Spacetime->get_local_metric(Local_State_Vector);
    Metric_type s_dr_Metric = this->p_Sim_Context->p_Spacetime->get_dr_local_metric(Local_State_Vector);

    /* = References for the sake of readability = */

    const double& r_source = Local_State_Vector[e_r];
    const double& theta_source = Local_State_Vector[e_theta];

    const auto& g = s_Metric.Metric;
    const auto& dr_g = s_dr_Metric.Metric;

    /* ========================================== */

    invert_metric(inv_metric, s_Metric.Metric);

    switch (this->s_Disk_params.Velocity_profile_type) {

    case e_Keplarian:

        /* This velocity profile is defined only for orbit radii > ISCO. */
        if (fabs(r_source) < this->p_Sim_Context->p_Spacetime->get_ISCO()[Inner]) { throw std::runtime_error("Disk with a Keplarian velocity profile extends below the ISCO orbit!"); }

        Omega = (-dr_g[e_t][e_phi] + sqrt(dr_g[e_t][e_phi] * dr_g[e_t][e_phi] - dr_g[e_t][e_t] * dr_g[e_phi][e_phi])) / dr_g[e_phi][e_phi];
        Normalization = g[e_t][e_t] + 2 * g[e_t][e_phi] * Omega + g[e_phi][e_phi] * Omega * Omega;

        this->Disk_Velocity[e_t] = sqrt(-1.0 / Normalization);
        this->Disk_Velocity[e_r] = 0.0;
        this->Disk_Velocity[e_theta] = 0.0;
        this->Disk_Velocity[e_phi] = this->Disk_Velocity[e_t] * Omega;

        break;

    default:

        rho = r_source * fabs(sin(theta_source));
        ell = sqrt(rho * rho * rho) / (1 + rho);

        /* I have noticed that this velocity profile becomes ill-defined in some places for the metric in the below "if" clause.
           I correct this by modifying the angular momentum profile by something that seems reasonable. */
        if (Janis_Newman_Winicour == this->p_Sim_Context->p_Init_Conditions->Metric_parameters.e_Spacetime) {

            double& gamma = p_Sim_Context->p_Init_Conditions->Metric_parameters.JNW_Gamma_Parameter;
            ell *= pow(1. - 2. / r_source / gamma, gamma);

        }
        else if (Wormhole == p_Sim_Context->p_Init_Conditions->Metric_parameters.e_Spacetime) {

            ell *= (1. - p_Sim_Context->p_Init_Conditions->Metric_parameters.R_throat / r_source);

        }

        u_t = -1.0 / sqrt(-(inv_metric[e_t][e_t] - 2 * inv_metric[e_t][e_phi] * ell + inv_metric[e_phi][e_phi] * ell * ell));
        u_phi = -u_t * ell;

        /* Convert U_source to contravariant components to compute the circular velocity profile */
        this->Disk_Velocity[e_t] = inv_metric[e_t][e_t] * u_t + inv_metric[e_t][e_phi] * u_phi;
        this->Disk_Velocity[e_r] = 0.0;
        this->Disk_Velocity[e_theta] = 0.0;
        this->Disk_Velocity[e_phi] = inv_metric[e_phi][e_phi] * u_phi + inv_metric[e_phi][e_t] * u_t;

        break;

    }

    if (isnan(this->Disk_Velocity[e_t]) or
        isinf(this->Disk_Velocity[e_t]) or
        isnan(this->Disk_Velocity[e_phi]) or
        isinf(this->Disk_Velocity[e_phi])) {

        throw std::runtime_error("Invalid disk velocity!");

    }

    return this->Disk_Velocity;

}