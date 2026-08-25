#pragma once
#include <iostream>
#include <format>
#include "Structs.h"
#include "Constants.h"
#include "General_math_functions.h"
#include "General_GR_functions.h"

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>

struct Disk_model_type {

private:

    gsl_spline2d* Spline_instance_density;

    gsl_interp_accel* Radial_interp_accelerator;
    gsl_interp_accel* Theta_interp_accelerator;

    gsl_function Function_to_solve;
    gsl_root_fsolver* Root_finder;

    von_Zeipel_cylinder_condition_wrapper_struct von_Zeipel_cylinder_condition_wrapper_params;

    Simulation_Context_type* p_Sim_Context{};

    double Disk_Velocity[4]{};

    double Geometric_to_cgs_density_convertor{};

    ~Disk_model_type();

    double get_disk_gas_pressure(const double density) const;

    double get_disk_mag_pressure(const double density, const double* const Local_State_Vector) const;

    double get_disk_internal_energy(const double density) const;

    double get_disk_profile(const Disk_profile_parameters_type* const p_Profile_parameters,
                            Profile_enums e_Profile_type) const;

    void get_phenomenological_mag_field(const double* const Local_State_Vector,
                                          const Metric_type* const p_Metric,
                                          Emission_medium_state_type* const Emission_medium_state) const;

    void get_numerical_mag_field(const double* const Local_State_Vector,
                                   const Metric_type* const p_Metric,
                                   Emission_medium_state_type* const Emission_medium_state);

    double get_Keplarian_ang_momentum_profile(double r_0) const;

    double get_disk_eq_ang_momentum_profile(double r_0) const;

    double get_disk_ang_momentum_profile(const double* const Local_State_Vector);

public:

    double get_von_Zeipel_cylinder_condition(const Metric_type Metric, double r_0) const;

    /* Holds all the model parameteres for the background accretion disk. */
    Disk_model_parameters_type s_Disk_params{};

    /*! Copies over the initial data from the Simulation Context struct to internal variables for the sake of convenience
     *
     *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct.
     *   \return Nothing.
     */
    Disk_model_type(Simulation_Context_type* p_Sim_Context);
    
    void get_magnetic_field(const double* const Local_State_Vector,
                            const Metric_type* const p_Metric,
                            Emission_medium_state_type* const Emission_medium_state);

    /*! Computes the accretion disk at the current photon position.
     *
     *   \param [in] State_Vector - Pointer to the current photon state vector.
     *   \param [out] p_Emission_medium_state - Pointer to the struct that holds the temperature and density of the accretion disk.
     *   \return Nothing.
     */
    void get_density_and_temperature(const double* const State_Vector,
                                     Emission_medium_state_type* const p_Emission_medium_state) const;

    bool is_inside_disk(const double* const State_Vector, Emission_medium_state_type* const Disk_State) const;

    const double* const get_disk_velocity(const double* const Local_State_Vector);

};
