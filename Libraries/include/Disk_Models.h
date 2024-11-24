#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include "Structs.h"
#include "Constants.h"
#include "Spacetimes.h"
#include "General_math_functions.h"
#include "General_GR_functions.h"
#include "gsl/gsl_sf_hyperg.h"

//! Evaluates the Planck function in the frequency domain in CGS units
/*! Evaluates the Planck function in the frequency domain in CGS units
 *
 *   \param [in] Frequency - Emission frequency in units [Hz].
 *   \param [in] Temperature - Emission medium temperature in units [K].
 *   \return The value of the Planck function in CGS.
 */
double get_planck_function_CGS(double Frequency, double Temperature);

class Novikov_Thorne_Model {

    private:

        double r_in;
        double r_out;
        double flux_integral_accuracy;
        Spacetime_Base_Class* p_Spacetime;
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
        Precomputed_e_pitch_angles s_Precomputed_e_pitch_angles{};

        /* The number of electron pitch angle values (in the range [0, pi]) to average over. */
        int Num_samples_to_avg{};

        /* Flag that controls weather to include the polarization calculations. */
        bool Include_polarization{};

        // ====================== Thermally Distributed synchrotron Fit Functions ====================== //

        //! Evaluates the thermal sychrotron emission fit functions.
        /*! Evaluates the thermal sychrotron emission fit functions, based on two sources:
        *       1) https://www.aanda.org/articles/aa/pdf/2022/11/aa44339-22.pdf for the unpolarized case.
        *       2) https://arxiv.org/pdf/1602.03184.pdf for the polarized case
        *
        *   \param [in] e_fit_functions - Enum that selects which fit functions to evaluate.
        *   \param [out] Emission_functions - Pointer to the array that holds the fit functions.
        *   \param [in] X - Dimensionless parameter that the fit functions depend on.
        *   \param [in] sqrt_X - Square root of X - computed ouside this function because it saves clock cycles when averaging over pitch angles.
        *   \param [in] cbrt_X - Cube root of X - computed ouside this function because it saves clock cycles when averaging over pitch angles.
        *   \return Nothing
        */
        void get_thermal_synchrotron_emission_fit_functions(Thermal_Syncotron_fit_selector e_fit_functions,
                                                          double Emission_functions[4],
                                                          double X,
                                                          double sqrt_X,
                                                          double cbrt_X);

        //! Evaluates the thermal sychrotron Faradey fit functions.
        /*! Evaluates the thermal sychrotron Faradey fit functions, based on this source: https://arxiv.org/pdf/1602.03184.pdf
        *
        *   \param [in] e_fit_functions - Enum that selects which fit functions to evaluate.
        *   \param [in] X - Dimensionless parameter that the fit functions depend on.
        *   \param [in] X_to_1_point_2 - pow(X,1.2) - computed ouside this function because it saves clock cycles when averaging over pitch angles.
        *   \param [in] X_frac - pow(X,1.035) - computed ouside this function because it saves clock cycles when averaging over pitch angles.
        *   \param [out] Faradey_fucntions - Pointer to the array that holds the evaluated Faradey functions.
        *   \return Nothing
        */
        void get_thermal_synchrotron_faradey_fit_functions(double const X,
                                                           double const X_to_1_point_2,
                                                           double const X_frac,
                                                           double* const Faradey_fucntions);
        //! Selects the approprite thermal synchrotron fit functions to evaluate, based on the "Include_polarization" flag.
        /*! Selects the approprite thermal synchrotron fit functions to evaluate, based on the "Include_polarization" flag.
        *
        *   \param [out] Emission_fucntions - ointer to the array that holds the emission fit functions
        *   \param [out] Faradey_functions - Pointer to the array that holds the Faradey fit functions.
        *   \param [in] Emission_args - Pointer to the struct that holds the arguments, necessary to evaluate the thermal synchrotron emission fit functions.
        *   \param [in] Faradey_args - Pointer to the struct that holds the arguments, necessary to evaluate the thermal synchrotron Faradey fit functions.
        *   \return Nothing
        */
        void get_thermal_synchrotron_fit_functions(double Emission_fucntions[STOKES_PARAM_NUM],
                                                  double Faradey_functions[STOKES_PARAM_NUM],
                                                  Thermal_emission_f_arguments* Emission_args,
                                                  Thermal_faradey_f_arguments* Faradey_args);

        //! Evaluates the thermal ensamble polarized synchrotron emission and Faradey functions.
        /*! Evaluates the thermal ensamble polarized synchrotron emission and Faradey functions.
         *
         *   \param [in] Density - The current emission medium density in [g/cm^3].
         *   \param [in] T_electron_dim - The current emission medium dimentionless temperature.
         *   \param [in] f_cyclo - The current cyclotron frequency in [Hz].
         *   \param [in] sin_pitch_angle - The sine of the angle between the magnetic field and the photon momentum 3-vector in the plasma frame.
         *   \param [in] cos_pitch_angle - The cosine of the angle between the magnetic field and the photon momentum 3-vector in the plasma frame.
         *   \param [out] Emission_functions - Vector to hold the emission functions.
         *   \param [out] Faradey_functions - Vector to hold the Faradey functions.
         *   \param [in] Emission_args - Sturct to hold the arguments for evaluating the emission fit functions @see get_thermal_synchrotron_fit_functions.
         *   \param [in] Faradey_args - Struct to hold the aruments for evaluating the Faradey fit function @see get_thermal_synchrotron_fit_functions.
         *   \return Nothing.
         */
        void evaluate_thermal_synchrotron_transfer_functions(double Density,
                                                            double T_electron_dim,
                                                            double f_cyclo,
                                                            double sin_pitch_angle,
                                                            double cos_pitch_angle,
                                                            double Emission_functions[STOKES_PARAM_NUM],
                                                            double Faradey_functions[STOKES_PARAM_NUM],
                                                            Thermal_emission_f_arguments Emission_args,
                                                            Thermal_faradey_f_arguments Faradey_args);


        // ====================== Kappa Distributed synchrotron Fit Functions ====================== //

        //! Evaluates the kappa sychrotron emission fit functions.
        /*! Evaluates the kappa sychrotron emission fit functions, based on the source: https://arxiv.org/pdf/1602.08749
        *
        *   \param [out] Emission_functions - Pointer to the array that holds the fit functions.
        *   \param [in] X - Dimensionless parameter that the fit functions depend on.
        *   \param [in] sqrt_X - Square root of X - computed ouside this function because it saves clock cycles when averaging over pitch angles.
        *   \param [in] cbrt_X - Cube root of X - computed ouside this function because it saves clock cycles when averaging over pitch angles.
        *   \param [in] X_to_7_over_20 - pow(X, 7./20) - computed ouside this function because it saves clock cycles when averaging over pitch angles.
        *   \param [in] kappa - The value of the dimensionless kappa parameter.
        *   \param [in] sin_emission_angle - The sin of the angle between the photon 3-momentum vector and the magnetic field vector, in the plasma frame.
        *   \param [in] T_electron_dim - The dimensionless electron temperature.
        *   \return Nothing
        */
        void get_kappa_synchrotron_emission_fit_functions(double* const Emission_functions,
                                                          const double X,
                                                          const double sqrt_X,
                                                          const double cbrt_X,
                                                          const double X_to_7_over_20,
                                                          const double kappa,
                                                          const double sin_emission_angle,
                                                          const double T_electron_dim);

        //! Evaluates the kappa sychrotron absorbtion fit functions.
        /*! Evaluates the kappa sychrotron absorbtion fit functions, based on the source: https://arxiv.org/pdf/1602.08749
        *
        *   \param [out] Absorbtion_functions - Pointer to the array that holds the fit functions.
        *   \param [in] X - Dimensionless parameter that the fit functions depend on.
        *   \param [in] sqrt_X - Square root of X - computed ouside this function because it saves clock cycles when averaging over pitch angles.
        *   \param [in] cbrt_X - Cube root of X - computed ouside this function because it saves clock cycles when averaging over pitch angles.
        *   \param [in] X_to_7_over_20 - pow(X, 7./20) - computed ouside this function because it saves clock cycles when averaging over pitch angles.
        *   \param [in] kappa - The value of the dimensionless kappa parameter.
        *   \param [in] sin_emission_angle - The sin of the angle between the photon 3-momentum vector and the magnetic field vector, in the plasma frame.
        *   \param [in] T_electron_dim - The dimensionless electron temperature.
        *   \return Nothing
        */
        void get_kappa_synchrotron_absorbtion_fit_functions(double* const Absorbtion_functions,
                                                            double const X,
                                                            double const sqrt_X,
                                                            double const cbrt_X,
                                                            double const X_to_7_over_20,
                                                            double const kappa,
                                                            double const sin_emission_angle,
                                                            double const T_electron_dim);

        //! Evaluates the kappa ensamble polarized synchrotron emission and Faradey functions.
        /*! Evaluates the kappa ensamble polarized synchrotron emission and Faradey functions.
         *
         *   \param [in] Density - The current emission medium density in [g/cm^3].
         *   \param [in] f_cyclo - The current cyclotron frequency in [Hz].
         *   \param [out] Emission_functions - Vector to hold the emission functions.
         *   \param [out] Faradey_functions - Vector to hold the Faradey functions.
         *   \param [out] Absorbtion_functions - Vector to hold the absorbtion functions.
         *   \param [in] Transfer_args - Sturct to hold the arguments for evaluating the emission fit functions @see get_kappa_synchrotron_fit_functions.
         *   \return Nothing.
         */
        void evaluate_kappa_synchrotron_transfer_functions(double Density,
                                                          double f_cyclo,
                                                          double Emission_functions[STOKES_PARAM_NUM],
                                                          double Faradey_functions[STOKES_PARAM_NUM],
                                                          double Absorbtion_functions[STOKES_PARAM_NUM],
                                                          Kappa_transfer_f_arguments Transfer_args);

        //! Selects the approprite kappa synchrotron fit functions to evaluate, based on the "Include_polarization" flag.
        /*! Selects the approprite kappa synchrotron fit functions to evaluate, based on the "Include_polarization" flag.
        *   NOTE: Currently Faradey functions are not implemented, so this functions just acts as a wrapper.
        *
        *   \param [out] Emission_fucntions - Pointer to the array that holds the emission fit functions
        *   \param [out] Faradey_functions - Pointer to the array that holds the Faradey fit functions.
        *   \param [out] Absorbtion_fucntions - Pointer to the array that holds the absorbtion fit functions.
        *   \param [in] args - Pointer to the struct that holds the arguments, necessary to evaluate the kappa synchrotron fit functions.
        *   \return Nothing
        */
        void get_kappa_synchrotron_fit_functions(double* const Emission_fucntions,
                                                 double* const Faradey_functions, 
                                                 double* const Absorbtion_fucntions,
                                                 const Kappa_transfer_f_arguments* const args);

        // ====================== Power Distributed synchrotron Fit Functions ====================== //


    public:

        // ====================== Accretion Disk State Functions ====================== //

        //! Computes the background accretion disk temperature
        /*! Computes the background accretion disk at temperature the current photon position.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The temperature in [K].
         */
        double get_disk_temperature(const double* const State_Vector);

        //! Computes the hotspot temperature
        /*! Computes the hotspot at temperature the current photon position.
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The temperature in [K].
         */
        double get_hotspot_temperature(const double* const State_Vector);

        //! Computes the emission medium's plasma 4-velocity
        /*! Computes the emission medium's plasma 4-velocity
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to call the metric function
         *   \param [in] Velocity_profile - Enum for the type of velocity profile
         *   \return Pointer to the 4-velocity vector
         */
        double* get_plasma_velocity(const double* const State_Vector, const Simulation_Context_type* const p_Sim_Context, Velocity_enums const Velocity_profile);

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
        double get_hotspot_density(const double* const State_Vector);

        //! Computes the magnetic field 4-vector in the coordinate and plasma frames.
        /*! Computes the magnetic field 4-vector, measured by a comoving obverver (with 4-velocity Plasma_velocity) in the following frames:
         *      1) That of a static observer (with 4-velocity n_mu = {1, 0, 0, 0} ) - a.e. the coordinate frame.
         *      2) The plasma rest frame.
         *
         *    NOTE: The magnitude of the magnetic field in these frames is different, because its not concerved under Lorentz boosts.
         *          In the plasma frame I set the geometry of the field, then scale it by B_Plasma_norm_CGS.
         *
         *    NOTE: The magnitudes of the magnetic fields in these frames are given in Gauss.
         *
         *   \param [out] Magnetic_fields - Struct that holds the magnetic field 4-vector in the two frames.
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to call the metric function for dot products.
         *   \param [in] Density - The current emission medium density - used to compute the field magnitude in the plasma frame.
         *   \param [in] Magnetization - The current emission medium magnetization - used to compute the field magnitude in the plasma frame.
         *   \return Nothing.
         */
        void get_magnetic_field(Magnetic_fields_type* Magnetic_fields,
                                const double* const State_Vector,
                                const Simulation_Context_type* const p_Sim_Context,
                                const double* const Plasma_Velocity,
                                double const Density,
                                double const Magnetization) ;

        //! Main "Selector" For The Transfer Functions.
        /*! Calculates the density, temperature, magnetic field and 4-velocity of the chosen emission medium and calls the respective transfer functions evaluation.
         *
         *   \param [in] State_Vector - The current photon state vector.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to access the initial conditions.
         *   \param [out] Emission_functions - Vector to hold the emission functions.
         *   \param [out] Faradey_functions - Vector to hold the Faradey functions.
         *   \param [out] Absorbtion_functions - Vector to hold the absorbtion functions.
         *   \param [in] Emission_medium - Enum that specifies which emission medium to evaluate.
         *   \return Nothing.
         */
        void get_radiative_transfer_functions(const double* const State_Vector,
                                              const Simulation_Context_type* const p_Sim_Context, 
                                              double* const Emission_functions,
                                              double* const Faradey_functions,
                                              double* const Absorbtion_functions,
                                              Emission_medium_enums Emission_medium);

        //! Computes the necessary variables for evaluating the thermal ensamble polarized synchrotron transfer functions.
        /*! Computes the necessary variables(cyclotron frequency, emission angles and so on) for evaluating the thermal ensamble polarized
         *   synchrotron transfer functions, based on the current photon position.
         *
         *   \param [in] State_Vector - The current photon state vector.
         *   \param [in] Plasma_velocity - The current emission medium plasma velocity.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to access the initial conditions.
         *   \param [out] Emission_functions - Vector to hold the emission functions.
         *   \param [out] Faradey_functions - Vector to hold the Faradey functions.
         *   \param [out] Absorbtion_functions - Vector to hold the absorbtion functions.
         *   \param [in] Density - The current density of the emission medium.
         *   \param [in] Temperature - The current temperrature of the emission medium.
         *   \param [in] B_field_coord_frame - The magnetic field 4-vector in the coordinate frame.
         *   \param [in] B_field_plasma_frame_norm - The norm of the magnetic field in the plasma frame.
         *   \return Nothing.
         */
        void get_thermal_synchrotron_transfer_functions(const double* const State_Vector,
                                                        const double* const Plasma_velocity,
                                                        const Simulation_Context_type* const p_Sim_Context,
                                                        double* const Emission_functions,
                                                        double* const Faradey_functions,
                                                        double* const Absorbtion_functions,
                                                        double  const Density,
                                                        double  const Temperature,
                                                        const double* const B_field_coord_frame,
                                                        double  const B_field_plasma_frame_norm);

        //! Computes the necessary variables for evaluating the kappa ensamble polarized synchrotron transfer functions.
        /*! Computes the necessary variables (cyclotron frequency, emission angles and so on) for evaluating the kappa ensamble polarized
         *   synchrotron transfer functions, based on the current photon position.
         *
         *   \param [in] State_Vector - The current photon state vector in geometric units.
         *   \param [in] Plasma_velocity - The current emission medium plasma velocity in geometric units.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to access the initial conditions.
         *   \param [out] Emission_functions - Vector to hold the emission functions in CGS units.
         *   \param [out] Faradey_functions - Vector to hold the Faradey functions in CGS units.
         *   \param [out] Absorbtion_functions - Vector to hold the absorbtion functions in CGS units.
         *   \param [in] Density - The current density of the emission medium in units [g/cm^3].
         *   \param [in] Temperature - The current temperrature of the emission medium in units [K].
         *   \param [in] B_field_coord_frame - The magnetic field 4-vector in the coordinate frame in units [G].
         *   \param [in] B_field_plasma_frame_norm - The norm of the magnetic field in the plasma frame in units [G].
         *   \return Nothing.
         */
        void get_kappa_synchrotron_transfer_functions(const double* const State_Vector,
                                                      const double* const Plasma_velocity,
                                                      const Simulation_Context_type* const p_Sim_Context,
                                                      double* const Emission_functions,
                                                      double* const Faradey_functions,
                                                      double* const Absorbtion_functions,
                                                      double  const Density,
                                                      double  const Temperature,
                                                      double* const B_field_coord_frame,
                                                      double  const B_field_plasma_frame_norm);

        //! Evaluates the phonomenological synchrotron transfer functions.
        /*! Evaluates the phonomenological synchrotron transfer functions.
         *
         *   \param [in] State_Vector - The current photon state vector
         *   \param [in] Plasma_velocity - The current emission medium plasma velocity.
         *   \param [in] p_Sim_Context - Pointer to the Simulation Context struct - used to access the initial conditions.
         *   \param [out] Emission_functions - Vector to hold the emission functions.
         *   \param [out] Faradey_functions - Vector to hold the Faradey functions.
         *   \param [out] Absorbtion_functions - Vector to hold the absorbtion functions.
         *   \param [in] Density - The current density of the emission medium.
         *   \return Nothing.
         */
        void get_phenomenological_synchrotron_functions(const double* const State_Vector,
                                                        const double* const Plasma_Veclocity,
                                                        const Simulation_Context_type* const p_Sim_Context, 
                                                        double* const Emission_functions,
                                                        double* const Faradey_functions,
                                                        double* const Absorbtion_functions,
                                                        const double Density);

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
        double get_electron_pitch_angle(const double* const B_field_local, 
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



