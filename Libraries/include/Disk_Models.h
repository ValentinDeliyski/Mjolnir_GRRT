#pragma once
#include <iostream>
#include <format>
#include "Structs.h"
#include "General_math_functions.h"

#include"gsl/gsl_interp2d.h"
#include"gsl/gsl_spline.h"
#include"gsl/gsl_spline2d.h"
#include"gsl/gsl_errno.h"

struct Disk_model_type {

    gsl_spline2d* Spline_instance_density;
    gsl_spline2d* Spline_instance_mag_field_r;
    gsl_spline2d* Spline_instance_mag_field_theta;
    gsl_spline2d* Spline_instance_mag_field_phi;

    gsl_interp_accel* Radial_interp_accelerator;
    gsl_interp_accel* Theta_interp_accelerator;

    /* Holds all the model parameteres for the background accretion disk. */
    Disk_model_parameters_type s_Disk_params{};

    //! Copies over the initial data from the Simulation Context struct to internal variables for the sake of convenience
    /*! Copies over the initial data from the Simulation Context struct to internal variables for the sake of convenience
     *
     *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct.
     *   \return Nothing.
     */
    Disk_model_type(Simulation_Context_type* p_Sim_Context);
    ~Disk_model_type();

    double get_disk_profile(const Disk_profile_parameters_type* const p_Profile_parameters,
                            Profile_enums e_Profile_type) const;

    //! Computes the accretion disk density
    /*! Computes the accretion disk at the current photon position.
     *
     *   \param [in] State_Vector - Pointer to the current photon state vector.
     *   \param [in] e_Disk_model - Enum that specifies which model to use for the disk.
     *   \param [out] p_Emission_medium_state - Pointer to the struct that holds the temperature and density of the accretion disk.
     *   \return Nothing.
     */
    void get_density_and_temperature(const double* const State_Vector,
                                     Disk_model_enums e_Disk_model,
                                     Emission_medium_state_type* const p_Emission_medium_state) const;


    bool is_inside_disk(const double* const State_Vector, Disk_model_enums e_Disk_model,Emission_medium_state_type* const Disk_State) const;

};
