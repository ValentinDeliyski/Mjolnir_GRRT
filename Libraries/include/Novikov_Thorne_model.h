#pragma once
#define _USE_MATH_DEFINES
#include "Structs.h"
#include "Spacetimes.h"
#include "Enumerations.h"
#include "General_GR_functions.h"
#include "General_math_functions.h"

class Novikov_Thorne_Model_class {

    private:

        /* The inner accretion disk radius in [M]. */
        double r_in;

        /* The outer accretion disk radius in [M]. */
        double r_out;

        /* The threshold value for the flux integral error estimate. */
        double flux_integral_accuracy;

        /* Pointer to the spacetime class. Stored in here so one does not have to pass it in as arguments to the functions.*/
        Spacetime_Base_Class* p_Spacetime;

        /* The spacetime enum. Certain calculations are spacetime specific. Stored in here so one does not have 
           to pass it in as arguments to the functions.*/
        Spacetime_enums e_Spacetime;

    public:

        //! Copies over initial conditions from the Simulation Context struct to internal class variables for the sake of convenicence.
        /*! Copies over initial conditions from the Simulation Context struct to internal class variables for the sake of convenicence.
         *
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct.
         *   \return Nothing.
         */
        Novikov_Thorne_Model_class(Simulation_Context_type* p_Sim_Context);

        //! Evaluates the Keplarian angular velocity of the Novikov-Thorne disk model.
        /*! Evaluates the Keplarian angular velocity of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The Keplarian angular velocity.
         */
        double Keplerian_angular_velocity(const double* const State_vector);

        //! Evaluates the radial derivative of the Keplarian angular velocity of the Novikov-Thorne disk model.
        /*! Evaluates the radial derivative of the Keplarian angular velocity of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The radial derivative of the Keplarian angular velocity.
         */
        double dr_Keplerian_angular_velocity(const double* const State_vector);

        //! Evaluates the redshift of the Novikov-Thorne disk model.
        /*! Evaluates the redshift of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The redshift.
         */
        double Redshift(const double* const State_vector, double r_obs, double theta_obs);

        //! Evaluates the energy of the Novikov-Thorne disk model.
        /*! Evaluates the energy of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The energy.
         */
        double disk_Energy(const double* const State_vector);

        //! Evaluates the angular momentum magnitude of the Novikov-Thorne disk model.
        /*! Evaluates the angular momentum magnitude of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The magnitude of the angular momentum.
         */
        double disk_Angular_Momentum(const double* const State_vector);

        //! Evaluates the integrand of the integral that appears in the flux expression of the Novikov-Thorne disk model.
        /*! Evaluates the integrand of the integral that appears in the flux expression of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The value of the integrand.
         */
        double Flux_integrand(const double* const State_vector);

        //! Evaluates the integral that appears in the flux expression of the Novikov-Thorne disk model.
        /*! Evaluates the integral that appears in the flux expression of the Novikov-Thorne disk model, using the adaptive Simpson method.
         *
         *   \param [in] r_in - The lower bound for the integral in units [M]
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The value of the integral term.
         */
        double solve_Flux_integral(double r_in, const double* const State_Vector, double tolerance);

        //! Evaluates the flux of the Novikov-Thorne disk model
        /*! Evaluates the flux of the Novikov-Thorne disk model
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The Novikov-Thorne flux in units [M_dot / M^2].
         */
        double get_flux(const double* const State_vector);

};
