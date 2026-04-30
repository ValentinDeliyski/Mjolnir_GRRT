#pragma once
#include <iostream>
#include <format>
#include "Structs.h"
#include "Constants.h"
#include "General_math_functions.h"
#include "General_GR_functions.h"

#include"gsl/gsl_interp2d.h"
#include"gsl/gsl_spline.h"
#include"gsl/gsl_spline2d.h"
#include"gsl/gsl_errno.h"

struct Disk_model_type {

    gsl_spline2d* Spline_instance_density;

    gsl_interp_accel* Radial_interp_accelerator;
    gsl_interp_accel* Theta_interp_accelerator;

    /* Holds all the model parameteres for the background accretion disk. */
    Disk_model_parameters_type s_Disk_params{};
    Simulation_Context_type* p_Sim_Context{};

    double Disk_Velocity[4]{};

    //! Copies over the initial data from the Simulation Context struct to internal variables for the sake of convenience
    /*! Copies over the initial data from the Simulation Context struct to internal variables for the sake of convenience
     *
     *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct.
     *   \return Nothing.
     */
    Disk_model_type(Simulation_Context_type* p_Sim_Context);
    ~Disk_model_type();

    void get_magnetic_field(const double* const Local_State_Vector,
                            const Metric_type* const p_Metric,
                            Emission_medium_state_type* const Emission_medium_state) const ;

    double get_disk_internal_energy(double density) const;

    double get_disk_profile(const Disk_profile_parameters_type* const p_Profile_parameters,
                            Profile_enums e_Profile_type) const;

    //! Computes the accretion disk density
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
