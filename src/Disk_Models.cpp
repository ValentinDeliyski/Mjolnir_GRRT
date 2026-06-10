#include "Disk_Models.h"

static double von_Zeipel_cylinder_condition_wrapper(double r_0, void* Params) {

    von_Zeipel_cylinder_condition_wrapper_struct* RHS_wrapper_params = (von_Zeipel_cylinder_condition_wrapper_struct*)Params;

    return RHS_wrapper_params->Disk_model->get_von_Zeipel_cylinder_condition(RHS_wrapper_params->Metric, r_0);

}

Disk_model_type::Disk_model_type(Simulation_Context_type* p_Sim_Context) {

    if (nullptr != p_Sim_Context) {

        this->s_Disk_params = p_Sim_Context->p_Init_Conditions->Disk_params;
        this->p_Sim_Context = p_Sim_Context;

    }
    else { throw std::runtime_error("Could not load the disk parameter struct! \n"); }

    if (this->s_Disk_params.e_Disk_model == Disk_model_enums::e_Numerical) {

        if (this->s_Disk_params.Numerical_disk_params.e_Spline_type == Spline_selection_enums::Custom_cubic) {

            throw std::runtime_error("Custom interpolants are not supported for numerical disks. Use GSL splines. \n"); 
    
        }

        if (this->s_Disk_params.Numerical_disk_params.e_Spline_type == Spline_selection_enums::GSL_linear) {

            this->Spline_instance_density = gsl_spline2d_alloc(gsl_interp2d_bilinear, this->s_Disk_params.Numerical_disk_params.R_coord_grid_size, this->s_Disk_params.Numerical_disk_params.Z_coord_grid_size);

        }
        else {

            this->Spline_instance_density = gsl_spline2d_alloc(gsl_interp2d_bicubic, this->s_Disk_params.Numerical_disk_params.R_coord_grid_size, this->s_Disk_params.Numerical_disk_params.Z_coord_grid_size);

        }

        this->Radial_interp_accelerator = gsl_interp_accel_alloc();
        this->Theta_interp_accelerator = gsl_interp_accel_alloc();

        gsl_spline2d_init(this->Spline_instance_density,
                          this->s_Disk_params.Numerical_disk_params.R_coord_grid.get(),
                          this->s_Disk_params.Numerical_disk_params.Z_coord_grid.get(),
                          this->s_Disk_params.Numerical_disk_params.Density_data.get(),
                          this->s_Disk_params.Numerical_disk_params.R_coord_grid_size,
                          this->s_Disk_params.Numerical_disk_params.Z_coord_grid_size);

    }

    this->Root_finder = gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);

    this->Function_to_solve.function = &von_Zeipel_cylinder_condition_wrapper;

    this->von_Zeipel_cylinder_condition_wrapper_params = { this, {} };

    this->Function_to_solve = { &von_Zeipel_cylinder_condition_wrapper,
                                &this->von_Zeipel_cylinder_condition_wrapper_params };

    this->Geometric_to_cgs_density_convertor = rho_0 / pow(this->p_Sim_Context->p_Init_Conditions->central_object_mass, 2) / M_PROTON_CGS;
}

Disk_model_type::~Disk_model_type() {

    if (this->s_Disk_params.e_Disk_model == Disk_model_enums::e_Numerical) {

        gsl_spline2d_free(this->Spline_instance_density);

        gsl_interp_accel_free(this->Radial_interp_accelerator);
        gsl_interp_accel_free(this->Theta_interp_accelerator);

    }

    gsl_root_fsolver_free(this->Root_finder);

}

double Disk_model_type::get_disk_gas_pressure(const double density) const {

    const double geometric_density = density / this->Geometric_to_cgs_density_convertor;

    const double& Gamma = this->s_Disk_params.Numerical_disk_params.Density_Polytrope_index;
    const double& K = this->s_Disk_params.Numerical_disk_params.Density_Polytrope_coeff;

    return K * std::pow(geometric_density, Gamma);

}

double Disk_model_type::get_disk_mag_pressure(const double density, const double* const Local_State_Vector) const {

    const double geometric_density = density / this->Geometric_to_cgs_density_convertor;

    const double K = 0.0001269061363;
    const double Gamma = 4. / 3;

    const double gas_pressure = this->get_disk_gas_pressure(density);
    const double internal_energy = this->get_disk_internal_energy(density);

    if (isinf(gas_pressure / density)) { return 0.0; }

    const double enthalpy = 1.0 + internal_energy + gas_pressure / density;

    const Metric_type s_Metric = this->p_Sim_Context->p_Spacetime->get_local_metric(Local_State_Vector);
    const double Metric_factor = s_Metric.Metric[e_t][e_phi] * s_Metric.Metric[e_t][e_phi] - s_Metric.Metric[e_t][e_t] * s_Metric.Metric[e_phi][e_phi];

    return K * pow(Metric_factor, Gamma - 1) * pow(geometric_density * enthalpy, Gamma);

}

double Disk_model_type::get_disk_internal_energy(const double density) const {

    const double geometric_density = density / this->Geometric_to_cgs_density_convertor;

    /* This assumes a ideal fluid, which is undergoing an iso-entropic process (Rezzolla (2.248)).
       The polytropic index of the polytropic EOS (Gamma) is assumed to be equal to the adiabatic index,
       which appears in the ideal fluid thermal EOS (2.228) */

    const double& Gamma = this->s_Disk_params.Numerical_disk_params.Density_Polytrope_index;
    const double& K = this->s_Disk_params.Numerical_disk_params.Density_Polytrope_coeff;

    return K / (Gamma - 1) * std::pow(geometric_density, Gamma - 1);
}

double Disk_model_type::get_disk_profile(const Disk_profile_parameters_type* const p_Profile_parameters,
                                         Profile_enums e_Profile_type) const {

    double Profile{};
    double Exponent_arg{};

    switch (e_Profile_type) {

    case Profile_enums::e_Power_law:

        Profile = std::pow(p_Profile_parameters->power_law_scale / p_Profile_parameters->radial_coordinate, p_Profile_parameters->power);
        break;

    case Profile_enums::e_Hybrid_power_gaussian:

        Exponent_arg = (p_Profile_parameters->gaussian_variable - p_Profile_parameters->gaussian_mean) / p_Profile_parameters->gaussian_std;

        Profile = std::pow(p_Profile_parameters->power_law_scale / p_Profile_parameters->radial_coordinate, p_Profile_parameters->power) * exp(-std::pow(Exponent_arg, 2) / 2);
        break;

    case Profile_enums::e_Gaussian:

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

double Disk_model_type::get_von_Zeipel_cylinder_condition(const Metric_type Metric, double r_0) const {

    const double Target_State_Vector[3] = { 0.0, r_0, std::numbers::pi / 2 };

    const Metric_type s_Target_Metric = this->p_Sim_Context->p_Spacetime->get_local_metric(Target_State_Vector);

    const auto& g = Metric.Metric;
    const auto& g_target = s_Target_Metric.Metric;

    const double quadratic_coeff = g[e_t][e_t] * g_target[e_t][e_phi] - g_target[e_t][e_t] - g[e_t][e_phi];
    const double linear_coeff = g[e_t][e_t] * g_target[e_phi][e_phi] - g_target[e_t][e_t] * g[e_phi][e_phi];
    const double free_coeff = g[e_t][e_phi] * g_target[e_phi][e_phi] - g_target[e_t][e_phi] * g[e_phi][e_phi];

    const double eq_ang_momentum_profile = this->get_disk_eq_ang_momentum_profile(r_0);

    return quadratic_coeff * pow(eq_ang_momentum_profile, 2) + linear_coeff * eq_ang_momentum_profile + free_coeff;

}

double Disk_model_type::get_Keplarian_ang_momentum_profile(double r_0) const {

    const double State_Vector[3] = {0.0, r_0, std::numbers::pi / 2 };

    const Metric_type s_Metric = this->p_Sim_Context->p_Spacetime->get_local_metric(State_Vector);
    const Metric_type s_dr_Metric = this->p_Sim_Context->p_Spacetime->get_dr_local_metric(State_Vector);

    const auto& g = s_Metric.Metric;
    const auto& dr_g = s_dr_Metric.Metric;

    const double B_coeff = -dr_g[e_t][e_phi] + sqrt(dr_g[e_t][e_phi] * dr_g[e_t][e_phi] - dr_g[e_t][e_t] * dr_g[e_phi][e_phi]);

    return -(B_coeff * g[e_phi][e_phi] + dr_g[e_phi][e_phi] * g[e_t][e_phi]) / (B_coeff * g[e_t][e_phi] + dr_g[e_phi][e_phi] * g[e_t][e_t]);

}

double Disk_model_type::get_disk_eq_ang_momentum_profile(double r_0) const {

    const double Keplarian_profile = this->get_Keplarian_ang_momentum_profile(r_0);
    const double Keplarian_profile_at_ISCO = this->get_Keplarian_ang_momentum_profile(this->s_Disk_params.r_ISCO);

    if (r_0 > this->s_Disk_params.r_ISCO) {

        return this->s_Disk_params.Ang_momentum_below_ISCO * pow(Keplarian_profile / Keplarian_profile_at_ISCO, this->s_Disk_params.Ang_momentum_exponent);

    }
    else {

        return this->s_Disk_params.Ang_momentum_below_ISCO;

    }

}

double Disk_model_type::get_disk_ang_momentum_profile(const double* const Local_State_Vector) {

    this->von_Zeipel_cylinder_condition_wrapper_params.Metric = this->p_Sim_Context->p_Spacetime->get_local_metric(Local_State_Vector);
   
    const double rho_coord = Local_State_Vector[e_r] * sin(Local_State_Vector[e_theta]);

    double offset_1 = 0.5;
    double offset_2 = 0.5;

    int iteration_num = 0;
    int root_finder_status = GSL_CONTINUE;

    double root_finder_lo_lim_test = this->get_von_Zeipel_cylinder_condition(this->von_Zeipel_cylinder_condition_wrapper_params.Metric, rho_coord - offset_1);
    double root_finder_hi_lim_test = this->get_von_Zeipel_cylinder_condition(this->von_Zeipel_cylinder_condition_wrapper_params.Metric, rho_coord + offset_2);

    while (root_finder_lo_lim_test * root_finder_hi_lim_test > 0 and iteration_num < 20) {

        if (std::abs(root_finder_lo_lim_test) > std::abs(root_finder_hi_lim_test)) {

            offset_2 += 1;

            root_finder_hi_lim_test = this->get_von_Zeipel_cylinder_condition(this->von_Zeipel_cylinder_condition_wrapper_params.Metric, rho_coord + offset_2);

        }
        else {

            offset_1 += 1;

            if (rho_coord - offset_1 < 0) {

                throw std::runtime_error("Could not find an approprite interval for the root finder in get_disk_velocity!");

            }

            root_finder_lo_lim_test = this->get_von_Zeipel_cylinder_condition(this->von_Zeipel_cylinder_condition_wrapper_params.Metric, rho_coord - offset_1);

        };

        iteration_num++;

    }

    /* The roots of this equation are expected to be rather close to r * sin(theta), so I give it an interval around that.
       That interval annoyingly enough can't be static, so the above logic finds an approprite one (the function at the endpoins has opposite signs). */
    gsl_root_fsolver_set(this->Root_finder, &this->Function_to_solve, rho_coord - offset_1, rho_coord + offset_2);

    while (root_finder_status == GSL_CONTINUE) {

        root_finder_status = gsl_root_fsolver_iterate(this->Root_finder);

        if ((root_finder_status != GSL_CONTINUE) and (root_finder_status != GSL_SUCCESS)) {

            throw std::runtime_error(std::format("GSL returned error code {} in get_disk_velocity!", root_finder_status));

        }

        const double gsl_root_lo = gsl_root_fsolver_x_lower(this->Root_finder);
        const double gsl_root_hi = gsl_root_fsolver_x_upper(this->Root_finder);

        root_finder_status = gsl_root_test_interval(gsl_root_lo, gsl_root_hi, 1e-8, 1e-8);

    };

    return this->get_disk_eq_ang_momentum_profile(gsl_root_fsolver_root(this->Root_finder));

}

void Disk_model_type::get_density_and_temperature(const double* const State_Vector,
                                                  Emission_medium_state_type* const p_Emission_medium_state) const {

    Disk_profile_parameters_type Density_profile_params{}, Temperature_profile_params{};

    const double& Upper_r_limit = this->s_Disk_params.Numerical_disk_params.R_coord_grid[this->s_Disk_params.Numerical_disk_params.R_coord_grid_size - 1];
    const double& Lower_r_limit = this->s_Disk_params.Numerical_disk_params.R_coord_grid[0];
    const double& Upper_z_limit = this->s_Disk_params.Numerical_disk_params.Z_coord_grid[this->s_Disk_params.Numerical_disk_params.Z_coord_grid_size - 1];

    const double& Gamma = this->s_Disk_params.Numerical_disk_params.Density_Polytrope_index;

    switch (this->s_Disk_params.e_Disk_model) {

    case Disk_model_enums::e_Numerical:

        /* ------------------------------------------------ Get the density profile ------------------------------------------------ */

        if (State_Vector[e_r] > Upper_r_limit or State_Vector[e_r] < Lower_r_limit or std::abs(State_Vector[e_r] * cos(State_Vector[e_theta])) > Upper_z_limit) {

            p_Emission_medium_state->Density = 0.0;

        }
        else {
            
            p_Emission_medium_state->Density = std::abs(gsl_spline2d_eval(this->Spline_instance_density,
                                                                          State_Vector[e_r],
                                                                          std::abs(State_Vector[e_r] * cos(State_Vector[e_theta])),
                                                                          this->Radial_interp_accelerator,
                                                                          this->Theta_interp_accelerator));

        }


        /* TODO: check this scaling */
        p_Emission_medium_state->Density *= this->Geometric_to_cgs_density_convertor;

        /* ------------------------------------------------ Get the temperature profile ------------------------------------------------ */

        p_Emission_medium_state->Temperature = M_PROTON_CGS / BOLTZMANN_CONST_CGS * (Gamma - 1) * this->get_disk_internal_energy(p_Emission_medium_state->Density);

        /* ---------------------------------------------------- Scale to CGS units ----------------------------------------------------- */

        /* TODO: check this scaling */
        p_Emission_medium_state->Temperature *= P_0 / rho_0;

        break;

    case Disk_model_enums::e_Phenom_RIAF_1:

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

        p_Emission_medium_state->Density = this->s_Disk_params.Common_RIAF_params.Electron_density_scale * this->get_disk_profile(&Density_profile_params, e_Hybrid_power_gaussian);

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

        p_Emission_medium_state->Temperature = this->s_Disk_params.Common_RIAF_params.Electron_temperature_scale * this->get_disk_profile(&Temperature_profile_params, e_Power_law);

        break;

    case Disk_model_enums::e_Phenom_RIAF_2:

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

        p_Emission_medium_state->Density = this->s_Disk_params.Common_RIAF_params.Electron_density_scale * this->get_disk_profile(&Density_profile_params, e_Hybrid_power_gaussian);

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

        p_Emission_medium_state->Temperature = this->s_Disk_params.Common_RIAF_params.Electron_temperature_scale * this->get_disk_profile(&Temperature_profile_params, e_Power_law);
        
        break;

    case Disk_model_enums::e_Phenom_RIAF_3:

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

        p_Emission_medium_state->Density = this->s_Disk_params.Common_RIAF_params.Electron_density_scale * this->get_disk_profile(&Density_profile_params, e_Gaussian);

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

        p_Emission_medium_state->Temperature = this->s_Disk_params.Common_RIAF_params.Electron_temperature_scale * this->get_disk_profile(&Temperature_profile_params, e_Gaussian);

        break;

    case Disk_model_enums::e_Colab_test_1:

        /* =============== This is the model used in https://iopscience.iop.org/article/10.3847/1538-4357/ab96c6/pdf ============== */

        /* ------------------------------------------------ Get the density profile ------------------------------------------------ */

        Density_profile_params.radial_coordinate = 0.0; // This is not used in this profile, so I set it to zero.
        Density_profile_params.gaussian_variable = State_Vector[e_r];
        Density_profile_params.gaussian_mean     = 0.0;
        Density_profile_params.gaussian_std      = this->s_Disk_params.Colab_test_1_params.Radial_scale;

        p_Emission_medium_state->Density = this->s_Disk_params.Colab_test_1_params.Density_scale * this->get_disk_profile(&Density_profile_params, e_Gaussian);

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

    case Disk_model_enums::e_Debug_constant_density:

        p_Emission_medium_state->Density = 1;
        p_Emission_medium_state->Temperature = 1;

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

    if (!this->s_Disk_params.Enable_flag or e_Novikov_Thorne == this->s_Disk_params.e_Disk_model) { return false; }

    this->get_density_and_temperature(State_Vector, Disk_State);

    return (Disk_State->Density / this->s_Disk_params.Max_disk_density > this->s_Disk_params.Threshold_relative_density);

}

void Disk_model_type::get_magnetic_field(const double* const Local_State_Vector,
                                         const Metric_type* const p_Metric,
                                         Emission_medium_state_type* const Emission_medium_state) {

    switch (this->s_Disk_params.e_Disk_model) {

    case Disk_model_enums::e_Numerical:

        this->get_numerical_mag_field(Local_State_Vector, p_Metric, Emission_medium_state);
        break;

    default: 

        this->get_phenomenological_mag_field(Local_State_Vector, p_Metric, Emission_medium_state);
        break;

    }

 }

void Disk_model_type::get_numerical_mag_field(const double* const Local_State_Vector,
                                              const Metric_type* const p_Metric,
                                              Emission_medium_state_type* const Emission_medium_state) {

    const double Mag_pressure = this->get_disk_mag_pressure(Emission_medium_state->Density, Local_State_Vector);

    double Plasma_velocity_covariant[4]{};
    Manipulate_index(p_Metric, Emission_medium_state->Plasma_Velocity, Plasma_velocity_covariant, Lower_index);

    const double Ang_momentum = -Plasma_velocity_covariant[e_phi] / Plasma_velocity_covariant[e_t];

    if (isinf(Ang_momentum) or isnan(Ang_momentum)) {

        throw std::runtime_error("Invalid angular momentum in get_numerical_mag_field()!");

    }

    const double Normalization = p_Metric->Metric[e_t][e_t] * Ang_momentum * Ang_momentum + 2 * p_Metric->Metric[e_t][e_phi] * Ang_momentum + p_Metric->Metric[e_phi][e_phi];

    Emission_medium_state->Magnetic_fields.B_field_plasma_frame[e_phi] = sqrt(2 * Mag_pressure / Normalization);
    Emission_medium_state->Magnetic_fields.B_field_plasma_frame[e_t] = Ang_momentum * Emission_medium_state->Magnetic_fields.B_field_plasma_frame[e_phi];

    /* TODO: Scale this thing so it gomes out in Gauss */

    Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm = sqrt(2 * Mag_pressure * P_0) / this->p_Sim_Context->p_Init_Conditions->central_object_mass;
    //Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm = 1;

}

void Disk_model_type::get_phenomenological_mag_field(const double* const Local_State_Vector,
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

    Phenomenological_mag_field_params_type Mag_field_params{};

    switch (this->s_Disk_params.e_Disk_model) {

    case Disk_model_enums::e_Colab_test_1:

        Mag_field_params = this->s_Disk_params.Colab_test_1_params.Mag_field_params;
        break;

    default:

        Mag_field_params = this->s_Disk_params.Common_RIAF_params.Mag_field_params;
        break;

    }

    /* ======================= References for the sake of readability ======================= */

    double (&B_eulerian)[4] = Emission_medium_state->Magnetic_fields.B_field_eulerian_frame;
    double (&B_plasma)[4] = Emission_medium_state->Magnetic_fields.B_field_plasma_frame;
    double& B_plasma_norm = Emission_medium_state->Magnetic_fields.B_field_plasma_frame_norm;

    const double& r = Local_State_Vector[e_r];

    const double& Power = Mag_field_params.Mag_field_power;
    const double& B_0 = Mag_field_params.Mag_field_B_0;
    const double& r_0 = Mag_field_params.Mag_field_r_0;

    /* ====================================================================================== */

    B_eulerian[e_t] = 0.0;

    switch (Mag_field_params.e_Mag_field_geometry) {

    case Magnetic_field_geometry_enums::Constant:

        B_eulerian[e_r]     = Mag_field_params.Mag_field_geometry[e_r - 1];
        B_eulerian[e_theta] = Mag_field_params.Mag_field_geometry[e_theta - 1];
        B_eulerian[e_phi]   = Mag_field_params.Mag_field_geometry[e_phi - 1];

        break;

    case Magnetic_field_geometry_enums::Vertical:

        B_eulerian[e_r]     = cos(Local_State_Vector[e_theta]);
        B_eulerian[e_theta] = -sin(Local_State_Vector[e_theta]);
        B_eulerian[e_phi]   = 0;

        break;

    case Magnetic_field_geometry_enums::Toroidal:

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

        B_plasma[e_t] += Emission_medium_state->Plasma_Velocity[left_idx] * B_eulerian[left_idx] / p_Metric->Lapse_function;
        
    }

    for (int index = 0; index < 4; index++) {

        B_plasma[index] = (B_eulerian[index] + p_Metric->Lapse_function * B_plasma[e_t] * Emission_medium_state->Plasma_Velocity[index]) / Lorentz_factor;

    }

    switch (Mag_field_params.e_Mag_field_magnitude_profile) {

    case Magnetic_field_magnitude_enums::Magnetization_based:

        B_plasma_norm = sqrt(Mag_field_params.Magnetization * C_LIGHT_CGS * C_LIGHT_CGS * Emission_medium_state->Density * M_PROTON_CGS * 4 * std::numbers::pi);
        break;

    case Magnetic_field_magnitude_enums::Power_law_based:

        B_plasma_norm = B_0 * pow(r_0 / r, Power);
        break;

    default:

        throw std::runtime_error("Unsupported magnetic field magnitude profile! \n");

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
        if (fabs(r_source) < this->s_Disk_params.r_ISCO) { throw std::runtime_error("Disk with a Keplarian velocity profile extends below the ISCO orbit!"); }

        Omega = (-dr_g[e_t][e_phi] + sqrt(dr_g[e_t][e_phi] * dr_g[e_t][e_phi] - dr_g[e_t][e_t] * dr_g[e_phi][e_phi])) / dr_g[e_phi][e_phi];
        Normalization = g[e_t][e_t] + 2 * g[e_t][e_phi] * Omega + g[e_phi][e_phi] * Omega * Omega;

        this->Disk_Velocity[e_t] = sqrt(-1.0 / Normalization);
        this->Disk_Velocity[e_r] = 0.0;
        this->Disk_Velocity[e_theta] = 0.0;
        this->Disk_Velocity[e_phi] = this->Disk_Velocity[e_t] * Omega;

        break;

    case e_von_Zeipel_cylinder:

        ell = this->get_disk_ang_momentum_profile(Local_State_Vector);

        u_t = -1.0 / sqrt(-(inv_metric[e_t][e_t] - 2 * inv_metric[e_t][e_phi] * ell + inv_metric[e_phi][e_phi] * ell * ell));
        u_phi = -u_t * ell;

        /* Convert U_source to contravariant components to compute the circular velocity profile */
        this->Disk_Velocity[e_t] = inv_metric[e_t][e_t] * u_t + inv_metric[e_t][e_phi] * u_phi;
        this->Disk_Velocity[e_r] = 0.0;
        this->Disk_Velocity[e_theta] = 0.0;
        this->Disk_Velocity[e_phi] = inv_metric[e_phi][e_phi] * u_phi + inv_metric[e_phi][e_t] * u_t;

        break;

    default:

        rho = r_source * fabs(sin(theta_source));
        ell = sqrt(rho * rho * rho) / (1 + rho);

        /* I have noticed that this velocity profile becomes ill-defined in some places for the metrics in the below "if" clause.
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