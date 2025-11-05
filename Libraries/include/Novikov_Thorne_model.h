#pragma once
#define _USE_MATH_DEFINES
#include "Structs.h"
#include "Enumerations.h"
#include "General_GR_functions.h"
#include "General_math_functions.h"

class Novikov_Thorne_Model_class {

    private:

        /* @brief The inner accretion disk radius in [M]. */
        double r_in;

        /* @brief The outer accretion disk radius in [M]. */
        double r_out;

        /* @brief The threshold value for the flux integral error estimate. */
        double flux_integral_accuracy;

        /* @brief The maximum number of iteration the adaptive Simpson integrator can recursively make. */
        int max_flux_integration_teps;

        /* @brief The current iteration number of the adaptive Simpson integrator. */
        int current_flux_integration_step;

        /* @brief Holds the contravariant disk velocity vector. This exists so I can pass the disk velocity to the generic redshift functions easier. */
        double Disk_veclovity_vector[4];

        double Source_polarization_vector[4];

        /*! @brief Specifies the direction of the magnetic field. */
        Magnetic_field_geometry_enums e_Mag_field_geometry{};

        /*! @brief Array that holds the magnetic field direction vector in the plasma frame
            (for the case e_Mag_field_geometry == Constant). */
        double Mag_field_geometry[3]{};

        /* @brief Pointer to the spacetime class. Stored in here so one does not have to pass it in as arguments to the functions.*/
        Spacetime_Base_Class* p_Spacetime;

        /* @brief The spacetime enum. Certain calculations are spacetime specific. Stored in here so one does not have 
           to pass it in as arguments to the functions.*/
        Spacetime_enums e_Spacetime;

    public:
      
        /*! @brief Copies over initial conditions from the Simulation Context struct to internal class variables for the sake of convenicence.
         *
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct.
         *   \return Nothing.
         */
        Novikov_Thorne_Model_class(Simulation_Context_type* p_Sim_Context);

        /*! @brief Evaluates the Keplarian angular velocity of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The Keplarian angular velocity.
         */
        double Keplerian_angular_velocity(const double* const State_vector);

        /*! @brief Evaluates the radial derivative of the Keplarian angular velocity of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The radial derivative of the Keplarian angular velocity.
         */
        double dr_Keplerian_angular_velocity(const double* const State_vector);

        /*! @brief Constructs the contravariant angular velocity vector, and returns a pointer to it.
         *
         *   \param [in] State_Vector - Current photon state vector..
         *   \return Pointer to the velocity vector.
         * 
         */
        double* get_disk_velocity_vector(const double* const State_Vector);

        /*! @brief Evaluates the energy of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The energy.
         */
        double disk_Energy(const double* const State_vector);

        /*! @brief Evaluates the angular momentum magnitude of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The magnitude of the angular momentum.
         */
        double disk_Angular_Momentum(const double* const State_vector);

        /*! @brief Evaluates the integrand of the integral that appears in the flux expression of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The value of the integrand.
         */
        double Flux_integrand(const double* const State_vector);

        /*! @brief Evaluates the integral that appears in the flux expression of the Novikov-Thorne disk model, using the adaptive Simpson method.
         *
         *   \param [in] r_in - The lower bound for the integral in units [M].
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \param [in] tolerance - the error threshold value for the adaptive Simpson integrator that this function uses.
         *   \return The value of the integral term.
         */
        double solve_Flux_integral(double r_in, const double* const State_Vector, double tolerance);

        /*! @brief Evaluates the flux of the Novikov-Thorne disk model
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The Novikov-Thorne flux in units [M_dot / M^2].
         */
        double get_flux(const double* const State_vector);

        double* Construct_coord_polarization_vector(const double* const State_Vector);

};
