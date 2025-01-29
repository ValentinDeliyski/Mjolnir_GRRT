#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include "Structs.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "General_math_functions.h"
#include "General_GR_functions.h"
#include "gsl/gsl_sf_hyperg.h"

class Novikov_Thorne_Model {

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
        Novikov_Thorne_Model(Simulation_Context_type* p_Sim_Context);

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

class Generic_Optically_Thin_Model {

    private:

        /* Holds all the model parameteres for the background accretion disk. */
        Disk_model_parameters_type s_Disk_params{};

        /* Holds all the model parameteres for hotspot. */
        Hotspot_model_parameters_type s_Hotspot_params{};

        /* Holds all the parameters for the emission models. */
        Emission_model_parameters_type s_Emission_params{};

        /* Holds all the precomupted electron pitch angles for use in averaging. */
        Precomputed_e_pitch_angles_type s_Precomputed_e_pitch_angles{};

        /* The number of electron pitch angle values (in the range [0, pi]) to average over. */
        int Num_samples_to_avg{};

        /* Flag that controls weather to include the polarization calculations. */
        bool Include_polarization{};

        // ====================== Thermally Distributed synchrotron Fit Functions ====================== //

        //! Evaluates the thermal sychrotron emission fit functions.
        /*! Evaluates the thermal sychrotron emission fit functions, based on the source: https://iopscience.iop.org/article/10.3847/1538-4357/ac1b28/pdf
        *
        *   \param [in] Transfer_arags - Pointer to a struct containing the dimensionless parameter that the fit functions depend on, its weird fractional powers, 
        *                               the dimentionless electron temperature, its weird powers, and the pitch angle.
        *   \param [out] Emission_functions - Pointer to the array that holds the fit functions.
        *   \return Nothing
        */
        void get_thermal_synchrotron_emission_fit_functions(const Thermal_transfer_f_arguments_type* const Transfer_arags,
                                                            double* const Emission_functions) const;


        void get_thermal_synchrotron_absorbtion_fit_functions(const Thermal_transfer_f_arguments_type* const Transfer_arags,
                                                              const Emission_medium_state_type* p_Emission_medium_state,
                                                              const double* const Emission_function,
                                                              double* const Absorbtion_function) const;

        //! Evaluates the thermal sychrotron Faradey fit functions.
        /*! Evaluates the thermal sychrotron Faradey fit functions, based on this source: https://arxiv.org/pdf/1602.03184.pdf
        *
        *   \param [in] Transfer_arags - Pointer to the struct containing the dimensionless parameter that the fit functions depend on, its weird fractional powers, 
        *                                the dimentionless electron temperature, its weird powers, and the pitch angle.
        *   \param [out] Faradey_fucntions - Pointer to the array that holds the evaluated Faradey functions.
        *   \return Nothing
        */
        void get_thermal_synchrotron_faradey_fit_functions(const Thermal_transfer_f_arguments_type* const Transfer_arags,
                                                           double* const Faradey_fucntions) const;

        // ====================== Kappa Distributed synchrotron Fit Functions ====================== //

        //! Evaluates the kappa sychrotron emission fit functions.
        /*! Evaluates the kappa sychrotron emission fit functions, based on the source: https://arxiv.org/pdf/1602.08749
        *
        *   
        *   \param [in] p_Transfer_args - Pointer to the struct containing the dimensionless parameter that the fit functions depend on, its weird fractional powers, 
        *                                 the dimentionless electron temperature, its weird powers, and the pitch angle.
        *    \param [out] Emission_functions - Pointer to the array that holds the fit functions.
        *   \return Nothing
        */
        void get_kappa_synchrotron_emission_fit_functions(const Kappa_transfer_f_arguments_type* const p_Transfer_args,
                                                          double* const Emission_functions);

        //! Evaluates the kappa sychrotron absorbtion fit functions.
        /*! Evaluates the kappa sychrotron absorbtion fit functions, based on the source: https://arxiv.org/pdf/1602.08749
        *
        *   
        *   \param [in] p_Transfer_args - Pointer to the struct containing the dimensionless parameter that the fit functions depend on, its weird fractional powers, 
        *                                 the dimentionless electron temperature, its weird powers, and the pitch angle.
        *   \param [out] Absorbtion_functions - Pointer to the array that holds the fit functions.
        *   \return Nothing
        */
        void get_kappa_synchrotron_absorbtion_fit_functions(const Kappa_transfer_f_arguments_type* const p_Transfer_args,
                                                            double* const Absorbtion_functions) const;

        // ====================== Phenomenologically Distributed synchrotron Fit Functions ====================== //

        void get_phenomenological_synchrotron_fit_functions(const Phenomenological_transfer_f_arguments_type* const p_Transfer_args,
                                                            Transfer_functions_type* const p_Transfer_functions) const;


        // ====================== Power Distributed synchrotron Fit Functions ====================== //


        //! Evaluates the kappa ensamble polarized synchrotron emission, absorbtion and Faradey functions.
        /*! Evaluates the kappa ensamble polarized synchrotron emission, absorbtion and Faradey functions.
         *
         *   \param [in] e_Ensamble_type - Enum that specifies the ensamble type.
         *   \param [in] p_Emission_medium_state - Pointer to the struct that holds the emmission medium state.
         *   \param [in] p_Transfer_args - Pointer to the sturct that holds the arguments for evaluating the emission fit functions.
         *   \param [out] p_Transfer_functions - Pointer to the struct that holds the transfer functions.
         *   \return Nothing.
         */
        void evaluate_synchrotron_transfer_functions(const Ensamble_enums e_Ensamble_type,
                                                     const Emission_medium_state_type* const p_Emission_medium_state,
                                                     const void* const p_Transfer_args,
                                                     const Simulation_Context_type* const p_Sim_Context,
                                                     Transfer_functions_type* const p_Transfer_functions);

    public:

        // ====================== Accretion Disk State Functions ====================== //

        Hotspot_position_type get_hotspot_position(const double* const State_Vector,
                                                   const Simulation_Context_type* const p_Sim_Context);

        //! Computes the background accretion disk temperature
        /*! Computes the background accretion disk at temperature the current photon position.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The temperature in [K].
         */
        double get_disk_temperature(const double* const State_Vector) const;

        //! Computes the hotspot temperature
        /*! Computes the hotspot at temperature the current photon position.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The temperature in [K].
         */
        double get_hotspot_temperature(const double* const State_Vector, const Simulation_Context_type* const p_Sim_Context);

        //! Computes the emission medium's plasma 4-velocity
        /*! Computes the emission medium's plasma 4-velocity
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to call the metric function
         *   \param [in] Velocity_profile - Enum for the type of velocity profile
         *   \return Pointer to the 4-velocity vector
         */
        double* get_plasma_velocity(const double* const State_Vector, 
                                    const Simulation_Context_type* const p_Sim_Context, 
                                    Velocity_enums const Velocity_profile,
                                    double const Radial_velocity_fraction);

        //! Computes the background accretion disk density
        /*! Computes the background accretion disk density at the current photon position.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The density in [g/cm^3].
         */
        double get_disk_density(const double* const State_Vector);

        //! Computes the hotspot density
        /*! Computes the hotspot density at the current photon position.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The density in [g/cm^3].
         */
        double get_hotspot_density(const double* const State_Vector, const Simulation_Context_type* const p_Sim_Context);

        //! Computes the magnetic field 4-vector in the coordinate and plasma frames.
        /*! Computes the magnetic field 4-vector, measured by a comoving obverver (with 4-velocity Plasma_velocity) in the following frames:
         *      1) That of a static observer (with 4-velocity n_mu = {1, 0, 0, 0} ) - a.e. the coordinate frame.
         *      2) The plasma rest frame.
         *
         *    NOTE: The magnitude of the magnetic field in these frames is different, because its not concerved under Lorentz boosts.
         *          In the plasma frame I set the geometry of the field, then scale it by B_Plasma_norm_CGS.
         *
         *    NOTE: The magnitudes of the magnetic fields in these frames are given in [G].
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to call the metric function for dot products.
         *   \param [out] Emission_medium_state - Pointer to the struct that holds the emission medium state.
         *   \return Nothing.
         */
        void get_magnetic_field(const double* const State_Vector,
                                const Simulation_Context_type* const p_Sim_Context,
                                Emission_medium_state_type* const p_Emission_medium_state) ;

        //! Main "Selector" For The Transfer Functions.
        /*! Calculates the density, temperature, magnetic field and 4-velocity of the chosen emission medium and calls the respective transfer functions evaluation.
         *
         *   \param [in] State_Vector - The current photon state vector.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to access the initial conditions.
         *   \param [out] p_Transfer_functions - Pointer to the struct that holds the radiative transfer functions.
         *   \param [in] Emission_medium - Enum that specifies which emission medium to evaluate.
         *   \return Nothing.
         */
        void get_radiative_transfer_functions(const double* const State_Vector,
                                              const Simulation_Context_type* const p_Sim_Context,
                                              const Emission_medium_enums Emission_medium,
                                              Transfer_functions_type* const p_Transfer_functions);

        //! Computes the necessary variables for evaluating the thermal ensamble polarized synchrotron transfer functions.
        /*! Computes the necessary variables(cyclotron frequency, emission angles and so on) for evaluating the thermal ensamble polarized
         *   synchrotron transfer functions, based on the current photon position.
         *
         *   \param [in] State_Vector - The current photon state vector in geometric units.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to access the initial conditions.
         *   \param [in] p_Emission_medium_state - Pointer to the struct that holds the emission medium state.
         *   \param [out] p_Transfer_functions - Pointer to the struct that holds the transfer functions.
         *   \return Nothing.
         */
        void get_thermal_synchrotron_transfer_functions(const double* const State_Vector,
                                                        const Simulation_Context_type* const p_Sim_Context,
                                                        const Emission_medium_state_type* const p_Emission_medium_state,
                                                        Transfer_functions_type* const p_Transfer_functions);

        //! Computes the necessary variables for evaluating the kappa ensamble polarized synchrotron transfer functions.
        /*! Computes the necessary variables (cyclotron frequency, emission angles and so on) for evaluating the kappa ensamble polarized
         *   synchrotron transfer functions, based on the current photon position.
         *
         *   \param [in] State_Vector - The current photon state vector in geometric units.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to access the initial conditions.
         *   \param [in] p_Emission_medium_state - Pointer to the struct that holds the emission medium state.
         *   \param [out] p_Transfer_functions - Pointer to the struct that holds the transfer functions.
         *   \return Nothing.
         */
        void get_kappa_synchrotron_transfer_functions(const double* const State_Vector,
                                                      const Simulation_Context_type* const p_Sim_Context,
                                                      const Emission_medium_state_type* const p_Emission_medium_state,
                                                      Transfer_functions_type* const p_Transfer_functions);

        //! Evaluates the phonomenological synchrotron transfer functions.
        /*! Evaluates the phonomenological synchrotron transfer functions.
         *
         *   \param [in] State_Vector - The current photon state vector in geometric units.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to access the initial conditions.
         *   \param [in] p_Emission_medium_state - Pointer to the struct that holds the emission medium state.
         *   \param [out] p_Transfer_functions - Pointer to the struct that holds the transfer functions.
         *   \return Nothing.
         */
        void get_phenomenological_synchrotron_functions(const double* const State_Vector,
                                                        const Simulation_Context_type* const p_Sim_Context, 
                                                        const Emission_medium_state_type* const Emission_medium_state,
                                                        Transfer_functions_type* const p_Transfer_functions);

        //! Computes the angle between the magnetic field and photon momentum 3-vectors in the plasma frame.
        /*! Computes the angle between the magnetic field and photon momentum 3-vectors in the plasma frame. There is a neat invariant way
         *   to do this by just operating on coordinate basis 4-vector using the projection tensor for an observer with 4-velocity = Plasma_velocity  .
         *
         *   \param [in] B_field_coord_frame - The magnetic field in the coordinate frame.
         *   \param [in] Plasma_velocity - The plasma velocity 4-vector.
         *   \param [in] State_Vector - Current photon state vector - used to get the photon momentum 4-vector.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to call the metric function for dot products.
         *   \return Cosine of the angle between the magnetic field and photon momentum 3-vectors in the plasma frame.
         */
        double get_electron_pitch_angle(const double* const B_field_coord_frame,
                                        const double* const Plasma_velocity,
                                        const double* const State_Vector, 
                                        const Simulation_Context_type* const p_Sim_Context);

        //! Precomputes the electron pitch angles and their weird powers to use in averaging.
        /*! Precomputes the electron pitch angles and their weird powers to use in averaging.
         *
         *   \param [in] p_Init_Conditions - Pointer to the struct that holds the initial conditions - used to determine how much memory to allocate.
         *   \return Nothing.
         */
        void precompute_electron_pitch_angles(Initial_conditions_type* p_Init_Conditions);

        //! Copies over the initial data from the Simulation Context struct to internal class variables for the sake of convenience
        /*! Copies over the initial data from the Simulation Context struct to internal class variables for the sake of convenience
         *
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct.
         *   \return Nothing.
         */
        int load_parameters(Simulation_Context_type* p_Sim_Context);

};



