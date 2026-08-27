#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <format>
#include "Structs.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "General_math_functions.h"
#include "General_GR_functions.h"

struct Hotspot_model_type {

    Hotspot_position_type Current_Position;

    double Current_Velocity[4];

    /* Holds all the model parameteres for hotspot. */
    Hotspot_model_parameters_type s_Hotspot_params;

    std::shared_ptr<Spacetime_Base_Class> p_Spacetime;

    //! Copies over the initial data from the Simulation Context struct to internal variables for the sake of convenience
    /*! Copies over the initial data from the Simulation Context struct to internal variables for the sake of convenience
     *
     *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct.
     *   \return Nothing.
     */
    Hotspot_model_type(Simulation_Context_type* p_Sim_Context);


    Hotspot_position_type get_hotspot_position(const double* const State_Vector);

    double get_hotspot_profile(const Hotspot_profile_parameters_type* const p_Profile_parameters,
                               Profile_enums e_Profile_type) const;

    double* get_hotspot_velocity(bool Eval_at_hotspot_center, const double* const Local_State_Vector);

    void get_magnetic_field(const double* const Local_State_Vector,
                            const Metric_type* const p_Metric,
                            Emission_medium_state_type* const Emission_medium_state) const;

    //! Computes the hotspot density
    /*! Computes the hotspot density at the current photon position.
     *
     *   \param [in] State_Vector - Pointer to the current photon state vector.
     *   \param [out] p_Emission_medium_state - Pointer to the struct that holds the current state of the hotspot.
     *   \return Nothing.
     */
    void get_density_and_temperature(const double* const State_Vector,
                                     Emission_medium_state_type* const p_Emission_medium_state);

    bool is_inside_hotspot(const double* const State_Vector, Emission_medium_state_type* const Hotspot_State);


};

