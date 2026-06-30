#pragma once
#define _USE_MATH_DEFINES
#include "IO_files.h"
#include "Constants.h"
#include "Disk_models.h"
#include "Spacetimes.h"
#include <filesystem>
#include <iostream>

File_manager_class::File_manager_class(Initial_conditions_type *p_Initial_Conditions) {

    this->p_Initial_Conditions  = p_Initial_Conditions;
    this->Truncate_files        = p_Initial_Conditions->File_manager_params.Truncate_files;

    this->sim_mode_1_ray_number = 0;

    this->Output_File_Path = "";
    this->Output_File = new std::ofstream();
}

void File_manager_class::get_geodesic_data(double J_data[], double p_theta_data[]) {

    std::ifstream geodesic_data;
    std::string line;

    double J_input{};
    double P_input{};

    geodesic_data.open(this->p_Initial_Conditions->File_manager_params.Sim_mode_1_imput_path, std::ios::in);

    while (true) {

        if (geodesic_data >> J_input >> P_input) {

            J_data[this->sim_mode_1_ray_number] = J_input;
            p_theta_data[this->sim_mode_1_ray_number] = P_input;
            
            this->sim_mode_1_ray_number += 1;
        }

        if (geodesic_data.eof()) {

            break;
        }
    }

    geodesic_data.close();
}

void File_manager_class::write_observer_metadata(std::ofstream* Output_file) {

    *Output_file << "------------------------------------------------------- Observer Metadata -------------------------------------------------------" << "\n"
                 << "Observervation Time [M]: " << p_Initial_Conditions->Observer_params.init_time << '\n'
                 << "Observer Distance [M]: " << p_Initial_Conditions->Observer_params.distance << '\n'
                 << "Observer Inclination [Deg]: " << p_Initial_Conditions->Observer_params.inclination * 180.0 / M_PI << '\n'
                 << "Observer Azimuth [Deg]: " << p_Initial_Conditions->Observer_params.azimuth * 180.0 / M_PI << '\n'
                 << "Observation Frequency [Hz]: " << p_Initial_Conditions->Observer_params.obs_frequency << '\n';

    switch (p_Initial_Conditions->Simulation_mode) {

    default:

        *Output_file << "Observation Window Dimentions (-X,+X,-Y,+Y) [M]: "
                     << this->p_Initial_Conditions->Observer_params.x_min << ","
                     << this->p_Initial_Conditions->Observer_params.x_max << ","
                     << this->p_Initial_Conditions->Observer_params.y_min << ","
                     << this->p_Initial_Conditions->Observer_params.y_max
                     << '\n'
                     << "Simulation Resolution: "
                     << this->p_Initial_Conditions->Observer_params.resolution_x
                     << " x "
                     << this->p_Initial_Conditions->Observer_params.resolution_y
                     << '\n';
        break;

    case 2:

        *Output_file << "Number Of Photons Per Param Value: " << this->sim_mode_1_ray_number << '\n'
                     << "Number Of Param Values: " << p_Initial_Conditions->Sim_mode_1_param_value_number << '\n';

        break;

    case 3:

        *Output_file << "X_init [M] = " << this->p_Initial_Conditions->Sim_mode_2_X_init << "\n" 
                     << "Y_init [M] = " << this->p_Initial_Conditions->Sim_mode_2_Y_init << "\n";

        break;

    }
}

void File_manager_class::write_metric_metadata(std::ofstream* Output_file) {

    *Output_file << "------------------------------------------------------- Metric Metadata -------------------------------------------------------" << "\n";

    *Output_file << "Spacetime [-]: " << this->Base_File_Names[this->p_Initial_Conditions->Metric_parameters.e_Spacetime] << "\n";

    switch (this->p_Initial_Conditions->Metric_parameters.e_Spacetime) {

    case Kerr:

        *Output_file << "Mass [M]: " << this->p_Initial_Conditions->Metric_parameters.Mass << '\n';
        *Output_file << "Spin Parameter [M]: " << this->p_Initial_Conditions->Metric_parameters.Spin << '\n';
        break;

    case Wormhole:

        *Output_file << "Mass [M]: " << this->p_Initial_Conditions->Metric_parameters.Mass << '\n';
        *Output_file << "Spin Parameter [M]: " << this->p_Initial_Conditions->Metric_parameters.Spin << "\n"
            << "Redshift Parameter [-]: " << this->p_Initial_Conditions->Metric_parameters.Redshift_Parameter << '\n';
        break;

    case Reg_Black_Hole:

        *Output_file << "Mass [M]: " << this->p_Initial_Conditions->Metric_parameters.Mass << '\n';
        *Output_file << "Parameter [M]: " << this->p_Initial_Conditions->Metric_parameters.RBH_Parameter << '\n';
        break;

    case Janis_Newman_Winicour:

        *Output_file << "Mass [M]: " << this->p_Initial_Conditions->Metric_parameters.Mass << '\n';
        *Output_file << "Gamma [-]: " << this->p_Initial_Conditions->Metric_parameters.JNW_Gamma_Parameter << '\n';
        break;

    case Einstein_Gauss_Bonnet:

        *Output_file << "Mass [M]: " << this->p_Initial_Conditions->Metric_parameters.Mass << '\n';
        *Output_file << "Gamma [M^2]: " << this->p_Initial_Conditions->Metric_parameters.GB_Gamma_Parameter << '\n';
        break;

    case BH_w_Dark_Matter:

        *Output_file << "Mass [M]: " << this->p_Initial_Conditions->Metric_parameters.Mass << '\n';
        *Output_file << "Halo Mass [M]: " << this->p_Initial_Conditions->Metric_parameters.Halo_Mass << '\n'
                     << "Halo Compactness [-]: " << this->p_Initial_Conditions->Metric_parameters.Compactness << '\n';

        break;

    case Numerical:

        *Output_file << "Numerical metric file path [-]: " << this->p_Initial_Conditions->Metric_parameters.Numerical_metric_params.Metric_file_path << "\n"
                     << "ADM Mass [M]: " << this->p_Initial_Conditions->Metric_parameters.Numerical_metric_params.M_ADM << "\n"
                     << "ADM spin parameter [M]: " << this->p_Initial_Conditions->Metric_parameters.Numerical_metric_params.a_ADM << "\n"
                     << "Metric Anzatz [-]: " << int(this->p_Initial_Conditions->Metric_parameters.Numerical_metric_params.e_Anzatz) << "\n";

        break;
    }

}

void File_manager_class::write_accretion_disk_metadata(std::ofstream* Output_file) {


    *Output_file << "------------------------------------------------------- Accretion Disk Metadata -------------------------------------------------------\n";

    switch (this->p_Initial_Conditions->Disk_params.e_Disk_model) {

    case e_Numerical:

        *Output_file << "Active disk model: Phenomenological_RIAF_1\n";
        break;

    case e_Phenom_RIAF_1:

        *Output_file << "Active disk model: Phenomenological_RIAF_1\n";
        break;

    case e_Phenom_RIAF_2:

        *Output_file << "Active disk model: Phenomenological_RIAF_2\n";
        break;

    case e_Phenom_RIAF_3:

        *Output_file << "Active disk model: Phenomenological_RIAF_3\n";
        break;

    case e_Colab_test_1:

        *Output_file << "Active disk model: Colaboration_test_1\n";
        break;

    case e_Novikov_Thorne:

        *Output_file << "Active disk model: Novikov-Thorne\n";
        break;
    }

    if (e_Phenom_RIAF_1 == this->p_Initial_Conditions->Disk_params.e_Disk_model or
        e_Phenom_RIAF_2 == this->p_Initial_Conditions->Disk_params.e_Disk_model or
        e_Phenom_RIAF_3 == this->p_Initial_Conditions->Disk_params.e_Disk_model) {

        
        *Output_file << "Maximum Density [g / cm^3]: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Electron_density_scale << "\n"
                     << "Maximum Temperature [K]: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Electron_temperature_scale << "\n";

        *Output_file << "--------------------------- Model Parameters\n"
                     << "Disk Opening Angle Parameter: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Disk_opening_angle << "\n"
                     << "Disk Density power law scale: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Density_power_law_scale << "\n"
                     << "Disk Density power law power: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Density_power_law_power << "\n"
                     << "Disk Density cutoff radius: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Density_cutoff_radius << "\n"
                     << "Disk Density cutoff scale: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Density_cutoff_scale << "\n"
                     << "Disk Temperature power law scale: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Temperature_power_law_scale << "\n"
                     << "Disk Temperature power law power: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Temperature_power_law_power << "\n"
                     << "Disk Temperature cutoff radius: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Temperature_cutoff_radius << "\n"
                     << "Disk Temperature cutoff scale: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Temperature_cutoff_scale << "\n";

        *Output_file << "--------------------------- Magnetic Field Parameters\n";

        switch (this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Mag_field_params.e_Mag_field_geometry) {

        case Toroidal:

            *Output_file << "Disk Magnetic field geometry: Toroidal\n";
            break;

        case Vertical:

            *Output_file << "Disk Magnetic field geometry: Vertical\n";
            break;

        case Constant:

            *Output_file << "Disk Magnetic field geometry: [" << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Mag_field_params.Mag_field_geometry[e_r - 1] << " "
                << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Mag_field_params.Mag_field_geometry[e_theta - 1] << " "
                << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Mag_field_params.Mag_field_geometry[e_phi - 1] << "]"
                << "\n";
            break;

        }

        switch (this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Mag_field_params.e_Mag_field_magnitude_profile) {

        case Magnetization_based:

            *Output_file << "Magnetic field magnitude profile: Magnetization based\n";
            *Output_file << "Disk Magnetization [-]: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Mag_field_params.Magnetization << "\n";

            break;

        case Power_law_based:

            *Output_file << "Magnetic field magnitude profile: Power law based\n"
                         << "Magnetic field scale: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Mag_field_params.Mag_field_B_0 << "\n"
                         << "Magnetic field radial scale: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Mag_field_params.Mag_field_r_0 << "\n"
                         << "Magnetic field power law: " << this->p_Initial_Conditions->Disk_params.Common_RIAF_params.Mag_field_params.Mag_field_power << "\n";

            break;
        }
   
    }
    else if (e_Novikov_Thorne == this->p_Initial_Conditions->Disk_params.e_Disk_model) {


        *Output_file << "--------------------------- Model Parameters\n"
                     << "Inner Disk Radius [M]: "
                     << this->p_Initial_Conditions->Disk_params.Novikov_Thorne_params.r_in
                     << "\n"
                     << "Outer Disk Radius [M]: "
                     << this->p_Initial_Conditions->Disk_params.Novikov_Thorne_params.r_out
                     << "\n";

        *Output_file << "--------------------------- Magnetic Field Parameters\n";

        *Output_file << "Disk Magnetic field geometry: [" << this->p_Initial_Conditions->Disk_params.Novikov_Thorne_params.Mag_field_geometry[e_r - 1] << " "
                                                          << this->p_Initial_Conditions->Disk_params.Novikov_Thorne_params.Mag_field_geometry[e_theta - 1] << " "
                                                          << this->p_Initial_Conditions->Disk_params.Novikov_Thorne_params.Mag_field_geometry[e_phi - 1] << "]"
                                                          << "\n";

    }
    else {

        *Output_file << "--------------------------- Model Parameters\n"
                     << "Disk Radial Scale: " << this->p_Initial_Conditions->Disk_params.Colab_test_1_params.Radial_scale << "\n"
                     << "Disk Vertical Scale: " << this->p_Initial_Conditions->Disk_params.Colab_test_1_params.Vertical_scale << "\n";

        *Output_file << "--------------------------- Magnetic Field Parameters\n";

        switch (this->p_Initial_Conditions->Disk_params.Colab_test_1_params.Mag_field_params.e_Mag_field_geometry) {

        case Toroidal:

            *Output_file << "Disk Magnetic field geometry: Toroidal\n";
            break;

        case Vertical:

            *Output_file << "Disk Magnetic field geometry: Vertical\n";
            break;

        case Constant:

            *Output_file << "Disk Magnetic field geometry: [" << this->p_Initial_Conditions->Disk_params.Colab_test_1_params.Mag_field_params.Mag_field_geometry[e_r - 1] << " "
                                                              << this->p_Initial_Conditions->Disk_params.Colab_test_1_params.Mag_field_params.Mag_field_geometry[e_theta - 1] << " "
                                                              << this->p_Initial_Conditions->Disk_params.Colab_test_1_params.Mag_field_params.Mag_field_geometry[e_phi - 1] << "]"
                                                              << "\n";
            break;

        }

        switch (this->p_Initial_Conditions->Disk_params.Colab_test_1_params.Mag_field_params.e_Mag_field_magnitude_profile) {

        case Magnetization_based:

            *Output_file << "Magnetic field magnitude profile: Magnetization based\n";
            *Output_file << "Disk Magnetization [-]: " << this->p_Initial_Conditions->Disk_params.Colab_test_1_params.Mag_field_params.Magnetization << "\n";

            break;

        case Power_law_based:

            *Output_file << "Magnetic field magnitude profile: Power law based\n"
                         << "Magnetic field scale: " << this->p_Initial_Conditions->Disk_params.Colab_test_1_params.Mag_field_params.Mag_field_B_0 << "\n"
                         << "Magnetic field radial scale: " << this->p_Initial_Conditions->Disk_params.Colab_test_1_params.Mag_field_params.Mag_field_r_0 << "\n"
                         << "Magnetic field power law: " << this->p_Initial_Conditions->Disk_params.Colab_test_1_params.Mag_field_params.Mag_field_power << "\n";

            break;
        }
    }
}

void File_manager_class::write_hotspot_metadata(std::ofstream* Output_file) {

    *Output_file << "------------------------------------------------------- Hotspot Parameters -------------------------------------------------------"
                 << "\n"
                 << "--------------------------- Density Model Parameters"
                 << "\n";

    switch (this->p_Initial_Conditions->Hotspot_params.Density_profile_type) {

    case e_Gaussian:

        *Output_file << "Density Profile: Gaussian"
                     << "\n"
                     << "Spread [M]: "
                     << this->p_Initial_Conditions->Hotspot_params.Profile_params.Density_gaussian_spread
                     << "\n";

        break;

    case e_Spherical:

        *Output_file << "Density Profile: Spherical"
                     << "\n"
                     << "Radius [M]: "
                     << this->p_Initial_Conditions->Hotspot_params.Profile_params.Radius
                     << "\n";

        break;

    case e_Hybrid_power_gaussian:

        *Output_file << "Density Profile: Hybrid power law gaussian"
                     << "\n"
                     << "Spread [M]: "
                     << this->p_Initial_Conditions->Hotspot_params.Profile_params.Density_gaussian_spread
                     << "\n"
                     << "Power law power [-]: "
                     << this->p_Initial_Conditions->Hotspot_params.Profile_params.Density_power_law_power
                     << "\n"
                     << "Power law scale [-]: "
                     << this->p_Initial_Conditions->Hotspot_params.Profile_params.Density_power_law_scale
                     << "\n";

        break;

    default:

        *Output_file << "Unsupported Density Profile!" << "\n";

        break;
    }

    *Output_file << "Maximum Density [g / cm^3]: "
                 << this->p_Initial_Conditions->Hotspot_params.Electron_density_scale
                 << "\n";

    *Output_file << "--------------------------- Temperature Model Parameters" << "\n";

    switch (this->p_Initial_Conditions->Hotspot_params.Temperature_profile_type) {

    case e_Spherical:

        *Output_file << "Temperature Profile: Spherical"
                     << "\n"
                     << "Radius [M]: "
                     << this->p_Initial_Conditions->Hotspot_params.Profile_params.Radius
                     << "\n";

        break;

    case e_Hybrid_power_gaussian:

        *Output_file << "Temperature Profile: Hybrid power law gaussian"
                     << "\n"
                     << "Spread [M]: "
                     << this->p_Initial_Conditions->Hotspot_params.Profile_params.Temperature_gaussian_spread
                     << "\n"
                     << "Power law power [-]: "
                     << this->p_Initial_Conditions->Hotspot_params.Profile_params.Temperature_power_law_power
                     << "\n"
                     << "Power law scale [-]: "
                     << this->p_Initial_Conditions->Hotspot_params.Profile_params.Temperature_power_law_scale
                     << "\n";

        break;

    default:

        *Output_file << "Temperature Profile : Gaussian"
            << "\n"
            << "Spread [M]: "
            << this->p_Initial_Conditions->Hotspot_params.Profile_params.Temperature_gaussian_spread
            << "\n";

        break;
    }

    *Output_file << "Maximum Temperature [K]: "
                 << this->p_Initial_Conditions->Hotspot_params.Electron_temperature_scale
                 << "\n";

    *Output_file << "--------------------------- Initial Hotspot Position" << "\n";

    *Output_file << "Distance [M]: "
                 << this->p_Initial_Conditions->Hotspot_params.Init_Position[e_r]
                 << "\n"
                 << "Inclination [Deg]: "
                 << this->p_Initial_Conditions->Hotspot_params.Init_Position[e_theta] * 180.0 / M_PI
                 << "\n"
                 << "Azimuth [Deg]: "
                 << this->p_Initial_Conditions->Hotspot_params.Init_Position[e_phi] * 180.0 / M_PI
                 << "\n";

    *Output_file << "Coordinate time offset [M]: " << this->p_Initial_Conditions->Hotspot_params.Profile_params.Coord_time_offset - this->p_Initial_Conditions->Observer_params.distance << "\n";


    *Output_file << "--------------------------- Magnetic field parameters " << "\n";

    *Output_file << "Hotspot Magnetization [-]: " << this->p_Initial_Conditions->Hotspot_params.Magnetization << "\n";

    switch (this->p_Initial_Conditions->Hotspot_params.e_Mag_field_geometry) {

    case Toroidal:

        *Output_file << "Hotspot Magnetic field geometry: Toroidal\n";
        break;

    case Vertical:

        *Output_file << "Hotspot Magnetic field geometry: Vertical\n";
        break;

    case Constant:

        *Output_file << "Hotspot Magnetic field geometry: [" << this->p_Initial_Conditions->Hotspot_params.Mag_field_geometry[e_r - 1] << " "
                                                     << this->p_Initial_Conditions->Hotspot_params.Mag_field_geometry[e_theta - 1] << " "
                                                     << this->p_Initial_Conditions->Hotspot_params.Mag_field_geometry[e_phi - 1] << "]"
                                                     << "\n";
        break;

    }

    *Output_file << "--------------------------- Misc Parameters\n";


    switch (this->p_Initial_Conditions->Hotspot_params.Ensamble_type) {

    case e_Phenomenological_ensamble:

        *Output_file << "Hotspot Ensamble: Phenomenological" << "\n";

        break;

    case e_Kappa_ensamble:
        *Output_file << "Hotspot Ensamble: Kappa" << "\n";

        break;

    default:

        *Output_file << "Hotspot Ensamble: Thermal" << "\n";

        break;
    }
}

void File_manager_class::write_emission_models_metadata(std::ofstream* Output_file) {

    *Output_file << "------------------------------------------------------- Emission models metadata -------------------------------------------------------" << "\n";

    if ((e_Phenomenological_ensamble == this->p_Initial_Conditions->Disk_params.Ensamble_type) or
        (e_Phenomenological_ensamble == this->p_Initial_Conditions->Hotspot_params.Ensamble_type)) {

        *Output_file << "--------------------------- Phenomenological ensamble parameters " << "\n";

        *Output_file << "Emission Power Law Exponent [-]: "
                     << this->p_Initial_Conditions->Emission_params.Phenomenological_emission_power_law
                     << "\n"
                     << "Absorbtion Coefficient [?]: "
                     << this->p_Initial_Conditions->Emission_params.Phenomenological_absorbtion_coeff
                     << "\n"
                     << "Source Function Power Law Exponent [-]: "
                     << this->p_Initial_Conditions->Emission_params.Phenomenological_source_f_power_law
                     << "\n"
                     << "Emission Scale [erg / (cm^3 s sr Hz)]: "
                     << this->p_Initial_Conditions->Emission_params.Phenomenological_emission_coeff
                     << "\n";

    }

    if ((e_Kappa_ensamble == this->p_Initial_Conditions->Disk_params.Ensamble_type) or
        (e_Kappa_ensamble == this->p_Initial_Conditions->Hotspot_params.Ensamble_type)) {

        *Output_file << "--------------------------- Kappa ensamble parameters " << "\n";

        *Output_file << "Kappa value [-]: " << this->p_Initial_Conditions->Emission_params.Kappa << "\n";
    }

    if ((e_Debug_constant_functions == this->p_Initial_Conditions->Disk_params.Ensamble_type) or
        (e_Debug_constant_functions == this->p_Initial_Conditions->Hotspot_params.Ensamble_type)) {

        *Output_file << "--------------------------- Debug ensamble parameters " << "\n";

        *Output_file << "Disk Ensamble: Debug"
                     << "\n"
                     << "j_I value [-]: "
                     << this->p_Initial_Conditions->Emission_params.Debug_j_I_value
                     << "\n"
                     << "j_Q value [-]: "
                     << this->p_Initial_Conditions->Emission_params.Debug_j_Q_value
                     << "\n"
                     << "j_U value [-]: "
                     << this->p_Initial_Conditions->Emission_params.Debug_j_U_value
                     << "\n"
                     << "j_B value [-]: "
                     << this->p_Initial_Conditions->Emission_params.Debug_j_V_value
                     << "\n"
                     << "alpha_I value [-]: "
                     << this->p_Initial_Conditions->Emission_params.Debug_alpha_I_value
                     << "\n"
                     << "alpha_Q value [-]: "
                     << this->p_Initial_Conditions->Emission_params.Debug_alpha_Q_value
                     << "\n"
                     << "alpha_U value [-]: "
                     << this->p_Initial_Conditions->Emission_params.Debug_alpha_U_value
                     << "\n"
                     << "alpha_V value [-]: "
                     << this->p_Initial_Conditions->Emission_params.Debug_alpha_V_value
                     << "\n"
                     << "rho_I value [-]: "
                     << this->p_Initial_Conditions->Emission_params.Debug_rho_I_value
                     << "\n"
                     << "rho_Q value [-]: "
                     << this->p_Initial_Conditions->Emission_params.Debug_rho_Q_value
                     << "\n"
                     << "rho_U value [-]: "
                     << this->p_Initial_Conditions->Emission_params.Debug_rho_U_value
                     << "\n"
                     << "rho_V value [-]: "
                     << this->p_Initial_Conditions->Emission_params.Debug_rho_V_value
                     << "\n";

    }

}

void File_manager_class::write_integrator_metadata(std::ofstream* Output_file) {

    *Output_file << "------------------------------------------------------- Geodesic Integrator metadata -------------------------------------------------------" << "\n";

    if (this->p_Initial_Conditions->Integrator_params.Geodesic_Step_Controller_Params.Use_adaptive_step) {

        *Output_file << "Step type [-]: Adaptive \n";
        *Output_file << "Max relative step increase [-]: " << this->p_Initial_Conditions->Integrator_params.Geodesic_Step_Controller_Params.Max_rel_step_increase << "\n";
        *Output_file << "Max stepsize [M]: " << this->p_Initial_Conditions->Integrator_params.Geodesic_Step_Controller_Params.Max_upper_stepsize << "\n";


    }
    else {

        *Output_file << "Step type [-]: Fixed \n";
        *Output_file << "Stepsize [M]: " << this->p_Initial_Conditions->Integrator_params.Geodesic_Step_Controller_Params.Init_stepzie;
    }

    *Output_file << "Geodesic integrator abs accuracy parameter [M]: " << this->p_Initial_Conditions->Integrator_params.Geodesic_Step_Controller_Params.RK_abs_accuracy << "\n";
    *Output_file << "Geodesic integrator rel accuracy parameter [M]: " << this->p_Initial_Conditions->Integrator_params.Geodesic_Step_Controller_Params.RK_rel_accuracy << "\n";
    *Output_file << "Max affine parameter [M]: " << this->p_Initial_Conditions->Integrator_params.Max_affine_param << "\n";

}

void File_manager_class::write_simulation_metadata() {

    *this->Output_File << "============================================================ SIMULATION METADATA ============================================================" << "\n";

    *this->Output_File << "Active Simulation Mode: " << static_cast<int>(p_Initial_Conditions->Simulation_mode) << '\n';  
    
    this->write_metric_metadata(this->Output_File);
    this->write_observer_metadata(this->Output_File); 
    this->write_accretion_disk_metadata(this->Output_File);
    this->write_hotspot_metadata(this->Output_File);
    this->write_emission_models_metadata(this->Output_File);
    this->write_integrator_metadata(this->Output_File);

    *this->Output_File << "============================================================ Simulation Results ============================================================"
                       << "\n";


    if (p_Initial_Conditions->Simulation_mode != Make_geodesic_log) {

        *this->Output_File << "Image X Coord [M],"
                           << "Image Y Coord [M],";

        if (e_Novikov_Thorne != this->p_Initial_Conditions->Disk_params.e_Disk_model){

           *this->Output_File << "Synchotron Intensity I [Jy/sRad],"
                              << "Synchotron Intensity Q [Jy/sRad],"
                              << "Synchotron Intensity U [Jy/sRad],"
                              << "Synchotron Intensity V [Jy/sRad],"
                              << "Final t Coordinate [M],"
                              << "Total Optical Depth [-],"
                              << "Total Faraday Q Depth [-],"
                              << "Total Faraday V Depth [-],";
        }
        else {

            *this->Output_File << "Disk Redshift [-],"
                               << "Disk Flux [M_dot/M^2]," 
                               << "Polarization vector X [-],"
                               << "Polarization vector Y [-],"
                               << "Source t Coord [M],"
                               << "Source r Coord [M],"
                               << "Source phi Coord [Rad],"
                               << "Radial Momentum (covariant),"
                               << "Theta Momentum (covariant),"
                               << "Phi Momentum (covariant),";
        }

        *this->Output_File << "Celestial Sphere Crossing Theta [Rad],"
                           << "Celestial Sphere Crossing Phi [Rad],";

    } else {

       *this->Output_File << "t_coord [M],"
                          << "r_coord [M],"
                          << "theta_coord [rad],"
                          << "phi_coord [rad],"
                          << "p_t [-],"
                          << "p_r [-],"
                          << "p_theta [rad/M],"
                          << "p_phi [rad/M],"
                          << "Integration Step [M],"
                          << "Affine Parameter [M],"
                          << "Geodesic State Error [-],"
                          << "Number of rejected steps [-],"
                          << "Synchotron Intensity I [Jy/sRad],"
                          << "Synchotron Intensity Q [Jy/sRad],"
                          << "Synchotron Intensity U [Jy/sRad],"
                          << "Synchotron Intensity V [Jy/sRad],"
                          << "ZAMO Polarization vector x [-],"
                          << "ZAMO Polarization vector y [-],";

       if (Spacetime_enums::Kerr == this->p_Initial_Conditions->Metric_parameters.e_Spacetime) {

            *this->Output_File << "PW Constant 1 [-]," << "PW Constant 2 [-],";

       }

    }

    *this->Output_File << '\n';
    
}

void File_manager_class::write_debug_metadata() {

    *this->Output_File << "============================================================ SIMULATION METADATA ============================================================" << "\n";

    *this->Output_File << "Active Simulation Mode: " << static_cast<int>(p_Initial_Conditions->Simulation_mode) << '\n';

    this->write_metric_metadata(this->Output_File);
    this->write_observer_metadata(this->Output_File);
    this->write_accretion_disk_metadata(this->Output_File);
    this->write_hotspot_metadata(this->Output_File);
    this->write_emission_models_metadata(this->Output_File);
    this->write_integrator_metadata(this->Output_File);

    *this->Output_File << "============================================================ Simulation Results ============================================================"
                       << "\n";

    *this->Output_File << "NT Flux r coord [M],"
                       << "NT Flux [M_dot/M^2]"
                       << "\n";

}

void File_manager_class::create_output_file() {

    // Create the path to the main results directory

    std::string Output_directory_path = this->p_Initial_Conditions->File_manager_params.Output_file_directory + "/" + this->p_Initial_Conditions->File_manager_params.Simulation_name;
    std::error_code error_code;

    if (!std::filesystem::exists(Output_directory_path)) {

        std::filesystem::create_directories(Output_directory_path, error_code);

    }

    if (0 != error_code.value()) {

        std::cout << "Could not create output directory!" << "\n";

        exit(ERROR);

    }

    std::filesystem::path dir(Output_directory_path);

    // Init the std::path variables where we will store the names of the output files for sim modes 1 and 2
    std::filesystem::path Image_file_name;

    // Init the std::path variables where we will store the names of the output files for sim mode 3
    std::filesystem::path Photon_log_name;

    // Init the std::path variables where we will store the names of the output file the debug sim mode
    std::filesystem::path Debug_file_name = "Debug_log";

    // Specify the output file extention
    std::filesystem::path file_extention(".txt");

    // Set weather we truncate the file upon opening or not
    auto open_type = std::ios::app;

    if (this->Truncate_files) {

        open_type = std::ios::trunc;

    }

    // Loop over all the files and populate the (so far empty) 
    
    if (this->p_Initial_Conditions->Simulation_mode == Make_geodesic_log) {

        if (0 == strcmp(static_cast<const char*>(this->p_Initial_Conditions->File_manager_params.Common_file_names.c_str()), "")) {

            Photon_log_name = this->Base_File_Names[this->p_Initial_Conditions->Metric_parameters.e_Spacetime] + "_photon_log";

        }
        else {

            Photon_log_name = this->p_Initial_Conditions->File_manager_params.Common_file_names + "_photon_log";

        }

        Photon_log_name.replace_extension(file_extention);
        this->Output_File_Path = dir / Photon_log_name;

        this->Output_File->open(this->Output_File_Path, open_type);

    }
    else if (this->p_Initial_Conditions->Simulation_mode == Image_generation or this->p_Initial_Conditions->Simulation_mode == Make_geodesic_sweep) {

        if (0 == strcmp(static_cast<const char*>(this->p_Initial_Conditions->File_manager_params.Common_file_names.c_str()), "")) {

            Image_file_name = this->Base_File_Names[this->p_Initial_Conditions->Metric_parameters.e_Spacetime];

        }
        else {

            Image_file_name = this->p_Initial_Conditions->File_manager_params.Common_file_names;

        }

        Image_file_name.replace_extension(file_extention);
        this->Output_File_Path = dir / Image_file_name;

        this->Output_File->open(this->Output_File_Path, open_type);

    }
    else {

        this->Output_File_Path = dir / Debug_file_name.replace_extension(file_extention);
        this->Output_File->open(this->Output_File_Path, open_type);

    }

    if (this->Truncate_files and !(this->p_Initial_Conditions->Simulation_mode == Debug_mode)) {

        // If we are truncating the file (and we are not in the debug sim mode), we should write the metadata to it.
        // The debug mode has its own header type
        this->write_simulation_metadata();

    }

    if (this->p_Initial_Conditions->Simulation_mode == Debug_mode) {

        this->write_debug_metadata();

    }

    this->Output_File->close();

}

void File_manager_class::write_image_data_to_file(Results_type* s_Ray_results) {

    *this->Output_File << s_Ray_results->Image_Coords[e_x]
                            << ","
                            << s_Ray_results->Image_Coords[e_y]
                            << "," 
                            << std::setprecision(15);

    if (e_Novikov_Thorne == this->p_Initial_Conditions->Disk_params.e_Disk_model) {

        *this->Output_File << s_Ray_results->Redshift_NT
                                << ","
                                << s_Ray_results->Flux_NT
                                << ","
                                << s_Ray_results->Projected_polarization_vector[e_x]
                                << ","
                                << s_Ray_results->Projected_polarization_vector[e_y]
                                << ","
                                << s_Ray_results->Thin_Disk_State_Vector[e_t]
                                << ","
                                << s_Ray_results->Thin_Disk_State_Vector[e_r]
                                << ","
                                << s_Ray_results->Thin_Disk_State_Vector[e_phi]
                                << ","
                                << s_Ray_results->Thin_Disk_State_Vector[e_p_r]
                                << ","
                                << s_Ray_results->Thin_Disk_State_Vector[e_p_theta]
                                << ","
                                << s_Ray_results->Thin_Disk_State_Vector[e_p_phi]
                                << ",";
    }
    else {

        *this->Output_File << s_Ray_results->Intensity[I] * CGS_TO_JANSKY
                                << ","
                                << s_Ray_results->Intensity[Q] * CGS_TO_JANSKY
                                << ","
                                << s_Ray_results->Intensity[U] * CGS_TO_JANSKY
                                << ","
                                << s_Ray_results->Intensity[V] * CGS_TO_JANSKY
                                << ","
                                << s_Ray_results->Final_State_Vector[e_t]
                                << ","
                                << s_Ray_results->Optical_Depth
                                << ","
                                << s_Ray_results->Faraday_Q_Depth
                                << ","
                                << s_Ray_results->Faraday_V_Depth
                                << ",";
    }

    *this->Output_File << s_Ray_results->Celestial_sphere_crossing_coords[e_theta]
                            << ","
                            << s_Ray_results->Celestial_sphere_crossing_coords[e_phi];

    *this->Output_File << '\n';

}

void File_manager_class::write_debug_data_to_file(Debug_mode_struct* Debug_results) {

    for (int idx = 0; idx < Debug_results->Array_length; idx++){

        *this->Output_File << Debug_results->NT_Flux_r_coord_array[idx]
                           << ","
                           << Debug_results->NT_Flux_integral_array[idx]
                           << '\n';
    }

}

void File_manager_class::log_photon_path(Results_type* p_Ray_results) {

    *this->Output_File << std::setprecision(15);

    for (int log_index = 1; log_index < p_Ray_results->Ray_log_struct.Log_length - 1; log_index++) {

        for (int state_index = 0; state_index < e_Full_state_size; state_index++) {

            *this->Output_File << p_Ray_results->Ray_log_struct.Ray_path_log_global[state_index + log_index * e_Full_state_size] << ",";

        }

        *this->Output_File << p_Ray_results->RK_integrator_debug_log.State_error_history[log_index] << ",";
        *this->Output_File << p_Ray_results->RK_integrator_debug_log.N_steps_rejected[log_index] << ",";

        if (1) {

            for (int stokes_index = I; stokes_index < e_Stokes_param_num; stokes_index++) {

                *this->Output_File << p_Ray_results->Ray_log_struct.Ray_emission_log[stokes_index][log_index - 1] << ",";

            }

            for (int Pol_component = 0; Pol_component < 2; Pol_component++) {

                *this->Output_File << 0.0 << ",";

            }


            if (Spacetime_enums::Kerr == this->p_Initial_Conditions->Metric_parameters.e_Spacetime) {

                for (int Pol_component = 0; Pol_component < 2; Pol_component++) {

                    *this->Output_File << 0.0 << ",";

                }

            }

        }
        else {

            size_t idx = log_index - (p_Ray_results->Ray_log_struct.Log_length - p_Ray_results->Ray_log_struct.Log_offet_at_disk_edge);

            for (int stokes_index = I; stokes_index < e_Stokes_param_num; stokes_index++) {

                *this->Output_File << p_Ray_results->Ray_log_struct.Ray_emission_log[stokes_index][idx] << ",";

            }

            for (int Pol_component = 0; Pol_component < 2; Pol_component++) {

                *this->Output_File << p_Ray_results->Ray_log_struct.Ray_polarization_log[Pol_component][idx] << ",";

            }

            if (Spacetime_enums::Kerr == this->p_Initial_Conditions->Metric_parameters.e_Spacetime) {

                for (int Pol_component = 0; Pol_component < 2; Pol_component++) {

                    *this->Output_File << p_Ray_results->Polarization_debug_log.PW_constant[Pol_component][idx] << ",";

                }

            }

        }

        *this->Output_File << '\n';
    }
   
};

void File_manager_class::open_output_file() {

    this->Output_File->open(this->Output_File_Path, std::ios::app);

}

void File_manager_class::close_output_file() {

    this->Output_File->close();
    
}