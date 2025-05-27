#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include "Structs.h"
#include "General_math_functions.h"

struct Hotspot_model_type {

    /* Holds all the model parameteres for hotspot. */
    Hotspot_model_parameters_type s_Hotspot_params;

    //! Copies over the initial data from the Simulation Context struct to internal variables for the sake of convenience
    /*! Copies over the initial data from the Simulation Context struct to internal variables for the sake of convenience
     *
     *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct.
     *   \return Nothing.
     */
    Hotspot_model_type(Simulation_Context_type* p_Sim_Context);


    Hotspot_position_type get_hotspot_position(const double* const State_Vector,
                                               const double* const Hotspot_Velocit) const;

    double get_hotspot_profile(const Hotspot_profile_parameters_type* const p_Profile_parameters,
                               Profile_enums e_Profile_type) const;

    //! Computes the hotspot density
    /*! Computes the hotspot density at the current photon position.
     *
     *   \param [in] State_Vector - Pointer to the current photon state vector.
     *   \param [in] Hotspot_Velocity - Pointer to the hotspot four-velocity.
     *   \param [out] p_Emission_medium_state - Pointer to the struct that holds the temperature and density of the hotspot.
     *   \return Nothing.
     */
    void get_density_and_temperature(const double* const State_Vector,
                                     const double* const Hotspot_Velocity,
                                     Emission_medium_state_type* const p_Emission_medium_state) const;

    bool is_inside_hotspot(const double* const State_Vector, const double* const Hotspot_Velocity, Emission_medium_state_type* const Hotspot_State) const;


};

