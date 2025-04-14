#pragma once
#include <iostream>
#include "Structs.h"
#include "General_math_functions.h"

struct Disk_model_type {

    /* Holds all the model parameteres for the background accretion disk. */
    Disk_model_parameters_type s_Disk_params{};

    //! Copies over the initial data from the Simulation Context struct to internal variables for the sake of convenience
    /*! Copies over the initial data from the Simulation Context struct to internal variables for the sake of convenience
     *
     *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct.
     *   \return Nothing.
     */
    Disk_model_type(Simulation_Context_type* p_Sim_Context);

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


};
