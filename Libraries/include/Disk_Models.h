#pragma once
#include "Structs.h"

class Spacetime_Base_Class;

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
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The value of the integral term.
         */
        double solve_Flux_integral(double r_in, const double* const State_Vector, double tolerance);

        //! Evaluates the flux of the Novikov-Thorne disk model
        /*! Evaluates the flux of the Novikov-Thorne disk model
         *
         *   \param [in] State_Vector - Current photon state vector - used to get the photon position.
         *   \return The Novikov-Thorne flux in [M_dot / M^2].
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

        void get_thermal_synchrotron_emission_fit_functions(Thermal_Syncotron_fit_selector e_fit_functions,
                                                          double Emission_functions[4],
                                                          double X,
                                                          double sqrt_X,
                                                          double cbrt_X);

        void get_thermal_synchrotron_faradey_fit_functions(double X,
                                                          double X_to_1_point_2,
                                                          double X_frac, 
                                                          double faradey_fucntions[STOKES_PARAM_NUM]);

        void get_thermal_synchrotron_fit_functions(double Emission_fucntions[STOKES_PARAM_NUM],
                                                  double Faradey_functions[STOKES_PARAM_NUM],
                                                  Thermal_emission_f_arguments* Emission_args,
                                                  Thermal_faradey_f_arguments* Faradey_args);

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

        void get_kappa_synchrotron_emission_fit_functions(double Emission_functions[STOKES_PARAM_NUM],
                                                         double X,
                                                         double sqrt_X,
                                                         double cbrt_X,
                                                         double X_to_7_over_20,
                                                         double kappa,
                                                         double sin_emission_angle,
                                                         double T_electron_dim);

        void get_kappa_synchrotron_absorbtion_fit_functions(double Absorbtion_functions[STOKES_PARAM_NUM],
                                                           double X,
                                                           double sqrt_X,
                                                           double cbrt_X,
                                                           double X_to_7_over_20,
                                                           double kappa,
                                                           double sin_emission_angle,
                                                           double T_electron_dim);

        void evaluate_kappa_synchrotron_transfer_functions(double Density,
                                                          double f_cyclo,
                                                          double Emission_functions[STOKES_PARAM_NUM],
                                                          double Faradey_functions[STOKES_PARAM_NUM],
                                                          double Absorbtion_functions[STOKES_PARAM_NUM],
                                                          Kappa_transfer_f_arguments Transfer_args);

        void get_kappa_synchrotron_fit_functions(double Emission_fucntions[STOKES_PARAM_NUM],
                                                double Faradey_functions[STOKES_PARAM_NUM], 
                                                double Absorbtion_fucntions[STOKES_PARAM_NUM],
                                                Kappa_transfer_f_arguments* args);

        // ====================== Power Distributed synchrotron Fit Functions ====================== //


    public:

        // ====================== Accretion Disk State Functions ====================== //

        double get_disk_temperature(const double* const State_Vector);

        double get_hotspot_temperature(const double* const State_Vector);

        double* get_plasma_velocity(const double* const State_Vector, const Simulation_Context_type* const p_Sim_Context, Velocity_enums const Velocity_profile);

        double get_disk_density(const double* const State_Vector);

        double get_hotspot_density(const double* const State_Vector);

        void get_magnetic_field(Magnetic_fields_type* Magnetic_fields,
                                const double* const State_Vector,
                                const Simulation_Context_type* const p_Sim_Context,
                                const double* const Plasma_Velocity,
                                double const Density,
                               double const Magnetization) ;

        // ======================= Radiative Fransfer Functions ======================= //

        void get_radiative_transfer_functions(const double* const State_Vector,
                                              const Simulation_Context_type* const p_Sim_Context, 
                                              double* const Emission_functions,
                                              double* const Faradey_functions,
                                              double* const Absorbtion_functions,
                                              Emission_medium_enums Emission_medium);

        void get_thermal_synchrotron_transfer_functions(const double* const State_Vector,
                                                       const double* const Plasma_velocity,
                                                       const Simulation_Context_type* const p_Sim_Context,
                                                       double* const Emission_functions,
                                                       double* const Faradey_functions,
                                                       double* const Absorbtion_functions,
                                                       double  const Density,
                                                       double  const Temperature,
                                                       double* const B_field_coord_frame,
                                                       double  const B_field_plasma_frame_norm);

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

        void get_phenomenological_synchrotron_functions(const double* const State_Vector,
                                                       const double* const Plasma_Veclocity,
                                                       const Simulation_Context_type* const p_Sim_Context, 
                                                       double* const Emission_functions,
                                                       double* const Faradey_functions,
                                                       double* const Absorbtion_functions,
                                                       const double Density);

        // ======================= Electron Pitch Angle Functions ======================= //

        double get_electron_pitch_angle(const double* const B_field_local, 
                                        const double* const Plasma_velocity,
                                        const double* const State_Vector, 
                                        const Simulation_Context_type* const p_Sim_Context);

        // =============================== Misc Functions =============================== //

        void precompute_electron_pitch_angles(Initial_conditions_type* p_Init_Conditions);

        int load_parameters(Simulation_Context_type* p_Sim_Context);

};



