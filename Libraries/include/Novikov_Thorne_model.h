#pragma once
#define _USE_MATH_DEFINES
#include "Structs.h"
#include "Enumerations.h"
#include "General_GR_functions.h"
#include "General_math_functions.h"
#include "gsl/gsl_integration.h"

class Novikov_Thorne_Model_class {

    private:

        gsl_spline* Flux_integral_spline_instance;
        gsl_interp_accel* Flux_integral_accelerator;

        gsl_integration_cquad_workspace* Flux_integral_workspace;

        gsl_function Flux_integral_fucntion_struct;

        /* @brief The inner accretion disk radius in [M]. */
        double r_in;

        /* @brief The outer accretion disk radius in [M]. */
        double r_out;

        /* @brief The threshold value for the flux integral error estimate. */
        double flux_integral_accuracy;

        /* @brief Holds the contravariant disk velocity vector. This exists so I can pass the disk velocity to the generic redshift functions easier. */
        double Disk_veclovity_vector[4];

        double Source_polarization_vector[4];

        /*! @brief Array that holds the magnetic field direction vector in the plasma frame
            (for the case e_Mag_field_geometry == Constant). */
        double Mag_field_geometry[3]{};

        /* @brief The spacetime enum. Certain calculations are spacetime specific. Stored in here so one does not have 
           to pass it in as arguments to the functions.*/
        Spacetime_enums e_Spacetime;

        /*! @brief Evaluates the flux of the Novikov-Thorne disk model
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The Novikov-Thorne flux in units [M_dot / M^2].
         */
        double get_Flux(double* State_vector);

    public:
     
        std::unique_ptr<double[]> Flux_integral_array;
        std::unique_ptr<double[]> Flux_r_coords;

        /* @brief Pointer to the spacetime class. Stored in here so one does not have to pass it in as arguments to the functions.*/
        std::shared_ptr<Spacetime_Base_Class> p_Spacetime;

        /*! @brief Copies over initial conditions from the Simulation Context struct to internal class variables for the sake of convenicence.
         *
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct.
         *   \return Nothing.
         */
        Novikov_Thorne_Model_class(Simulation_Context_type* p_Sim_Context);

        ~Novikov_Thorne_Model_class();

        /*! @brief Evaluates the Keplarian angular velocity of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The Keplarian angular velocity.
         */
        double get_Disk_Angular_Velocity(const double* const State_vector) const;

        /*! @brief Evaluates the radial derivative of the Keplarian angular velocity of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The radial derivative of the Keplarian angular velocity.
         */
        double get_dr_Disk_Angular_Velocity(const double* const State_vector) const;

        /*! @brief Constructs the contravariant angular velocity vector, and returns a pointer to it.
         *
         *   \param [in] State_Vector - Current photon state vector..
         *   \return Pointer to the velocity vector.
         * 
         */
        double* get_Disk_Velocity_Vector(const double* const State_Vector);

        /*! @brief Evaluates the energy of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The energy.
         */
        double get_Disk_Energy(const double* const State_vector) const;

        /*! @brief Evaluates the angular momentum magnitude of the Novikov-Thorne disk model.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The magnitude of the angular momentum.
         */
        double get_Disk_Angular_Momentum(const double* const State_vector) const;

        double get_dr_Disk_Angular_Momentum(const double* const Local_State_Vector) const;

        double get_Interpolated_Flux(const double* const Local_State_Vector) const;

        double* Construct_coord_polarization_vector(const double* const State_Vector);

};
