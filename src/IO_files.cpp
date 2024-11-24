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
    this->sim_mode_2_ray_number = 0;
}

void File_manager_class::get_geodesic_data(double J_data[], double p_theta_data[]) {

    std::ifstream geodesic_data;
    std::string line;


    double J_input{};
    double P_input{};

    geodesic_data.open(this->p_Initial_Conditions->File_manager_params.Sim_mode_2_imput_path, std::ios::in);

    while (true) {

        if (geodesic_data >> J_input >> P_input) {

            J_data[this->sim_mode_2_ray_number] = J_input;
            p_theta_data[this->sim_mode_2_ray_number] = P_input;
            
            this->sim_mode_2_ray_number += 1;
        }

        if (geodesic_data.eof()) {

            break;
        }
    }

    geodesic_data.close();
}

void File_manager_class::write_simulation_metadata() {

    std::ofstream* Output_file;
    int Output_file_number = ORDER_NUM;

    switch (this->p_Initial_Conditions->Simulation_mode){

    case 3:

        Output_file = &this->Log_Output_File;
        Output_file_number = 1;

        break;

    default:

        Output_file = this->Image_Output_files;

        break;
    }

    for (int Image_order = direct; Image_order <= Output_file_number - 1; Image_order += 1) {

        *(Output_file + Image_order) << "============================================================ SIMULATION METADATA ============================================================"
                                        << "\n"
                                        << "Spacetime: "
                                        << this->Base_File_Names[this->p_Initial_Conditions->Metric_params.e_Spacetime]
                                        << "\n";

        Metric_parameters_type& Parameters = this->p_Initial_Conditions->Metric_params;

        switch (this->p_Initial_Conditions->Metric_params.e_Spacetime) {

        case Kerr:

            *(Output_file + Image_order) << "Spin Parameter [M]: " << Parameters.Spin << '\n';
            break;

        case Wormhole:

            *(Output_file + Image_order) << "Spin Parameter [M]: " << Parameters.Spin << "\n"
                                            << "Redshift Parameter [-]: " << Parameters.Redshift_Parameter << '\n';
            break;

        case Reg_Black_Hole:

            *(Output_file + Image_order) << "Parameter [M]: " << Parameters.RBH_Parameter << '\n';
            break;

        case Janis_Newman_Winicour:

            *(Output_file + Image_order) << "Gamma [-]: " << Parameters.JNW_Gamma_Parameter << '\n';
            break;

        case Einstein_Gauss_Bonnet:

            *(Output_file + Image_order) << "Gamma [M^2]: " << Parameters.GB_Gamma_Parameter << '\n';
            break;

        case BH_w_Dark_Matter:

            *(Output_file + Image_order) << "Halo Mass [M]: " << Parameters.Halo_Mass << '\n'
                                            << "Halo Compactness [-]: " << Parameters.Compactness << '\n';
            break;
        }

        

        *(Output_file + Image_order) << "Active Simulation Mode: " << p_Initial_Conditions->Simulation_mode << '\n';

        if (this->p_Initial_Conditions->Simulation_mode != 3) {

            *(Output_file + Image_order) << "Image Order [-]: " << Image_order << "\n";

        }             
     
        *(Output_file + Image_order) << "------------------------------------------------------- Observer Parameters -------------------------------------------------------" << "\n"
                                        << "Observer Distance [M]: " << p_Initial_Conditions->Observer_params.distance << '\n'
                                        << "Observer Inclination [Deg]: " << p_Initial_Conditions->Observer_params.inclination * 180.0 / M_PI << '\n'
                                        << "Observer Azimuth [Deg]: " << p_Initial_Conditions->Observer_params.azimuth * 180.0 / M_PI << '\n'
                                        << "Observation Frequency [Hz]: " << p_Initial_Conditions->Observer_params.obs_frequency << '\n';

        switch (p_Initial_Conditions->Simulation_mode) {

        case 1:

            *(Output_file + Image_order) << "Observation Window Dimentions (-X,+X,-Y,+Y) [M]: "
                                            << this->p_Initial_Conditions->Observer_params.x_min << ","
                                            << this->p_Initial_Conditions->Observer_params.x_max << ","
                                            << this->p_Initial_Conditions->Observer_params.y_min << ","
                                            << this->p_Initial_Conditions->Observer_params.y_max
                                            << '\n'
                                            << "Simulation Resolutoin: "
                                            << this->p_Initial_Conditions->Observer_params.resolution_x
                                            << " x "
                                            << this->p_Initial_Conditions->Observer_params.resolution_y
                                            << '\n';
            break;

        case 2:

            *(Output_file + Image_order) << "Number Of Photons Per Param Value: " << this->sim_mode_2_ray_number << '\n'
                                            << "Number Of Param Values: " << p_Initial_Conditions->Sim_mode_2_param_value_number << '\n';

            break;

        case 3:

            *(Output_file + Image_order) << "X_init [M] = " << this->p_Initial_Conditions->Sim_mode_3_X_init << "\n" 
                                         << "Y_init [M] = " << this->p_Initial_Conditions->Sim_mode_3_Y_init << "\n";

            break;


        default:

            *(Output_file + Image_order) << "Unsupported simulation mode!" << "\n";

        }


        *(Output_file + Image_order) << "------------------------------------------------------- Accretion Disk Parameters -------------------------------------------------------"
                                        << "\n"
                                        << "--------------------------- Density Model Parameters"
                                        << "\n";

        /*
        
        --------------------------------------- Print the accretion disk density parameters ---------------------------------------
        
        */

        switch (this->p_Initial_Conditions->Disk_params.Density_profile_type) {

        case e_Power_law_profile:

            *(Output_file + Image_order) << "Density Profile: Power law"
                                            << "\n"
                                            << "Disk Opening Angle [tan(angle)]: "
                                            << this->p_Initial_Conditions->Disk_params.Power_law_disk_opening_angle
                                            << "\n"
                                            << "Density R_0 [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Power_law_density_R_0
                                            << "\n"
                                            << "Density R_Cutoff [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Power_law_density_R_cutoff 
                                            << "\n"
                                            << "Density Cutoff Scale [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Power_law_density_cutoff_scale
                                            << "\n";
                                            

            break;

        case e_Exponential_law_profile:

            *(Output_file + Image_order) << "Density Profile: Exponential law"
                                            << "\n"
                                            << "Density Height Scale [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Exp_law_density_height_scale
                                            << "\n"
                                            << "Density Radial Scale [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Exp_law_density_radial_scale
                                            << "\n";

            break;

        default:

            *(Output_file + Image_order) << "Unsupported Density Profile!" << "\n";

            break;

        }

        *(Output_file + Image_order) << "Maximum Density [g / cm^3]: "
                                        << this->p_Initial_Conditions->Disk_params.Electron_density_scale
                                        << "\n";
   
        /*
        
        --------------------------------------- Print the accretion disk temperature parameters ---------------------------------------
        
        */

        *(Output_file + Image_order) << "--------------------------- Temperature Model Parameters"
                                        << "\n";

        switch (this->p_Initial_Conditions->Disk_params.Temperature_profile_type) {

        case e_Power_law_profile:

            *(Output_file + Image_order) << "Temperature Profile: Power law"
                                            << "\n"
                                            << "Temperature R_0 [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Power_law_temperature_R_0
                                            << "\n"
                                            << "Temperature R_Cutoff [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Power_law_temperature_R_cutoff << "\n"
                                            << "Temperature Cutoff Scale [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Power_law_temperature_cutoff_scale << "\n";

            break;

        case e_Exponential_law_profile:

            *(Output_file + Image_order) << "Temperature Profile: Exponential law"
                                            << "\n"
                                            << "Temperature Height Scale [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Exp_law_temperature_height_scale
                                            << "\n"
                                            << "Temperature Radial Scale [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Exp_law_temperature_radial_scale
                                            << "\n";

            break;

        default:

            *(Output_file + Image_order) << "Unsupported Temperature Profile!" << "\n";

            break;
        }

        *(Output_file + Image_order) << "Maximum Temperature [K]: "
                                        << this->p_Initial_Conditions->Disk_params.Electron_temperature_scale
                                        << "\n";

        /*
        
        --------------------------------------- Print the disk ensamble parameters ---------------------------------------
        
        */

        *(Output_file + Image_order) << "--------------------------- Disk Synchrotron Emission Model Parameters"
                                        << "\n";

        switch (this->p_Initial_Conditions->Disk_params.Ensamble_type) {

        case e_Phenomenological_ensamble:

            *(Output_file + Image_order) << "Disk Ensamble: Phenomenological"
                                            << "\n"
                                            << "Emission Power Law Exponent [-]: "
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

            break;

        case e_Thermal_ensamble:

            *(Output_file + Image_order) << "Disk Ensamble: Thermal"
                                            << "\n";

            break;

        case e_Kappa_ensamble:
            *(Output_file + Image_order) << "Disk Ensamble: Kappa"
                                            << "\n"
                                            << "Kappa value [-]: "
                                            << this->p_Initial_Conditions->Emission_params.Kappa
                                            << "\n";

            break;

        default:

            *(Output_file + Image_order) << "Unsupported Ensamble!" << "\n";

            break;
        }

        *(Output_file + Image_order) << "Disk Magnetization [-]: "
                                        << this->p_Initial_Conditions->Disk_params.Magnetization << "\n"
                                        << "Disk Magnetic Field Geometry [-]: "
                                        << "[" << this->p_Initial_Conditions->Disk_params.Mag_field_geometry[0] << " "
                                        << this->p_Initial_Conditions->Disk_params.Mag_field_geometry[1] << " "
                                        << this->p_Initial_Conditions->Disk_params.Mag_field_geometry[2] << "]"
                                        << "\n";

        /*
        
        ================================================ Print the hotspot density parameters ================================================
        
        */

        *(Output_file + Image_order) << "------------------------------------------------------- Hotspot Parameters -------------------------------------------------------"
                                        << "\n" 
                                        << "--------------------------- Density Model Parameters"
                                        << "\n";

        switch (this->p_Initial_Conditions->Hotspot_params.Density_profile_type) {

        case e_Gaussian_profile:

            *(Output_file + Image_order) << "Density Profile: Gaussian"
                                            << "\n"
                                            << "Spread [M]: "
                                            << this->p_Initial_Conditions->Hotspot_params.Density_spread
                                            << "\n";
                                            
            break;

        case e_Spherical_profile:

            *(Output_file + Image_order) << "Density Profile: Spherical"
                                            << "\n"
                                            << "Radius [M]: "
                                            << this->p_Initial_Conditions->Hotspot_params.Radius
                                            << "\n";            

            break;

        default:

            *(Output_file + Image_order) << "Unsupported Density Profile!" << "\n";

            break;
        }

        *(Output_file + Image_order) << "Maximum Density [g / cm^3]: "
                                        << this->p_Initial_Conditions->Hotspot_params.Electron_density_scale
                                        << "\n";

        /*

        --------------------------------------- Print the hotspot temperature parameters ---------------------------------------

        */

        *(Output_file + Image_order) << "--------------------------- Temperature Model Parameters"
                                        << "\n";

        switch (this->p_Initial_Conditions->Hotspot_params.Temperature_profile_type) {

        case e_Gaussian_profile:

            *(Output_file + Image_order) << "Temperature Profile : Gaussian"
                                            << "\n"
                                            << "Spread [M]: "
                                            << this->p_Initial_Conditions->Hotspot_params.Temperature_spread
                                            << "\n";
                                            
            break;

        case e_Spherical_profile:

            *(Output_file + Image_order) << "Temperature Profile: Spherical"
                                            << "\n"
                                            << "Radius [M]: "
                                            << this->p_Initial_Conditions->Hotspot_params.Radius
                                            << "\n";

            break;

        default:

            *(Output_file + Image_order) << "Unsupported Temperature Profile!" << "\n";

            break;
        }

        *(Output_file + Image_order) << "Maximum Temperature [K]: "
                                        << this->p_Initial_Conditions->Hotspot_params.Electron_temperature_scale
                                        << "\n";

        /*
        
        --------------------------------------- Print the hotspot ensamble parameters ---------------------------------------
        
        */

        *(Output_file + Image_order) << "--------------------------- Hotspot Synchrotron Emission Model Parameters"
                                        << "\n";

        switch (this->p_Initial_Conditions->Hotspot_params.Ensamble_type) {

        case e_Phenomenological_ensamble:

            *(Output_file + Image_order) << "Hotspot Ensamble: Phenomenological"
                                            << "\n"
                                            << "Emission Power Law Exponent [-]: "
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

            break;

        case e_Thermal_ensamble:

            *(Output_file + Image_order) << "Hotspot Ensamble: Thermal"
                                            << "\n";

            break;

        case e_Kappa_ensamble:
            *(Output_file + Image_order) << "Hotspot Ensamble: Kappa"
                                            << "\n"
                                            << "Kappa value [-]: "
                                            << this->p_Initial_Conditions->Emission_params.Kappa
                                            << "\n";

            break;

        default:

            *(Output_file + Image_order) << "Unsupported Ensamble!" << "\n";

            break;
        }

        *(Output_file + Image_order) << "Hotspot Magnetization [-]: "
                                        << this->p_Initial_Conditions->Disk_params.Magnetization << "\n"
                                        << "Hotspot Magnetic Field Geometry [-]: "
                                        << "[" << this->p_Initial_Conditions->Disk_params.Mag_field_geometry[0] << " "
                                        << this->p_Initial_Conditions->Disk_params.Mag_field_geometry[1] << " "
                                        << this->p_Initial_Conditions->Disk_params.Mag_field_geometry[2] << "]"
                                        << "\n";

        *(Output_file + Image_order) << "--------------------------- Hotspot Position"
                                        << "\n";

        *(Output_file + Image_order) << "Distance [M]: "
                                        << this->p_Initial_Conditions->Hotspot_params.Position[e_r - 1]
                                        << "\n"
                                        << "Inclination [Deg]: "
                                        << this->p_Initial_Conditions->Hotspot_params.Position[e_theta - 1] * 180.0 / M_PI
                                        << "\n"
                                        << "Azimuth [Deg]: "
                                        << this->p_Initial_Conditions->Hotspot_params.Position[e_phi - 1] * 180.0 / M_PI
                                        << "\n";


        *(Output_file + Image_order) << "------------------------------------------------------- Novikov - Thorne Model Parameters -------------------------------------------------------"
                                        << "\n";

        /*
        
        --------------------------------------- Novikov - Thorne Disk Parameters ---------------------------------------
        
        */

        if (this->p_Initial_Conditions->NT_params.evaluate_NT_disk){

            *(Output_file + Image_order) << "Inner Disk Radius [M]: "
                                            << this->p_Initial_Conditions->NT_params.r_in
                                            << "\n"
                                            << "Outer Disk Radius [M]: "
                                            << this->p_Initial_Conditions->NT_params.r_out
                                            << "\n";
        }
        else {

            *(Output_file + Image_order) << "Novikov - Thorne Disk Evaluation Is Disabled." << "\n";

        }

        *(Output_file + Image_order) << "------------------------------------------------------- Simulation Results -------------------------------------------------------"
                                        << "\n";


        if (p_Initial_Conditions->Simulation_mode != 3) {

            *(Output_file + Image_order) << "Image X Coord [M],"
                                            << " "
                                            << "Image Y Coord [M],"
                                            << " "
                                            << "Novikov-Thorne Disk Redshift [-],"
                                            << " "
                                            << "Novikov-Thorne Flux [M_dot/M^2],"
                                            << " "
                                            << "Synchotron Intensity I [Jy/sRad],"
                                            << " "
                                            << "Synchotron Intensity Q [Jy/sRad],"
                                            << " "
                                            << "Synchotron Intensity U [Jy/sRad],"
                                            << " "
                                            << "Synchotron Intensity V [Jy/sRad]";

            if (p_Initial_Conditions->Simulation_mode == 2) {

                *(Output_file + Image_order) << ", Source r Coord [M],"
                    << " "
                    << "Source phi Coord [Rad],"
                    << " "
                    << "Radial Momentum (covariant),"
                    << " "
                    << "Theta Momentum (covariant),"
                    << " "
                    << "Phi Momentum (covariant),"
                    << " ";

                switch (this->p_Initial_Conditions->Metric_params.e_Spacetime) {

                case Kerr:

                    *(Output_file + Image_order) << "Spin Parameter";
                    break;

                case Wormhole:

                    *(Output_file + Image_order) << "Spin Parameter,"
                        << " "
                        << "Redshift Parameter";
                    break;

                case Reg_Black_Hole:

                    *(Output_file + Image_order) << "Parameter";
                    break;

                case Janis_Newman_Winicour:

                    *(Output_file + Image_order) << "Gamma";
                    break;

                case Einstein_Gauss_Bonnet:

                    *(Output_file + Image_order) << "Gamma";
                    break;

                case BH_w_Dark_Matter:

                    *(Output_file + Image_order) << "Halo Mass,"
                        << " "
                        << "Halo Compactness";
                    break;
                }
            }

        }else{

           *(Output_file + Image_order) << "t_coord [M],"
                                                << " "
                                                << "r_coord [M],"
                                                << " "
                                                << "theta_coord [rad],"
                                                << " "
                                                << "phi_coord [rad],"
                                                << " "
                                                << "p_t [-],"
                                                << " "
                                                << "p_r [-],"
                                                << " "
                                                << "p_theta [rad/M],"
                                                << " "
                                                << "p_phi [rad/M],"
                                                << " "
                                                << "Integration Step [M]"
                                                << " "
                                                << "Synchotron Intensity I [Jy/sRad],"
                                                << " "
                                                << "Synchotron Intensity Q [Jy/sRad],"
                                                << " "
                                                << "Synchotron Intensity U [Jy/sRad],"
                                                << " "
                                                << "Synchotron Intensity V [Jy/sRad]";

        }

        *(Output_file + Image_order) << '\n';
    }
}

void File_manager_class::open_image_output_files() {

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

    // Init the std::path variables where we will store the names of the output files
    std::filesystem::path Image_file_names[ORDER_NUM];

    // Init the std::path variables of the full file paths
    std::filesystem::path Image_full_path[ORDER_NUM]{};

    // Specify the output file extention
    std::filesystem::path file_extention(".txt");

    // Set weather we truncate the file upon opening or not
    auto open_type = std::ios::app;

    if (this->Truncate_files) {

        open_type = std::ios::trunc;

    }

    // Loop over all the files and populate the (so far empty) 
    
    for (int File_Index = 0; File_Index <= ORDER_NUM - 1; File_Index += 1) {

        if (0 == strcmp(static_cast<const char*>(this->p_Initial_Conditions->File_manager_params.Common_file_names.c_str()), "")) {

            Image_file_names[File_Index] = this->Base_File_Names[this->p_Initial_Conditions->Metric_params.e_Spacetime]
                                         + "_n"
                                         + std::to_string(File_Index);

        }
        else {

            Image_file_names[File_Index] = this->p_Initial_Conditions->File_manager_params.Common_file_names
                                         + "_n"
                                         + std::to_string(File_Index);

        }

        Image_file_names[File_Index].replace_extension(file_extention);
        Image_full_path[File_Index] = dir / Image_file_names[File_Index];
       

        this->Image_Output_files[File_Index].open(Image_full_path[File_Index], open_type);

    }


    if (this->Truncate_files) {

        // If we are truncating the file, we should write the metadata to it
        this->write_simulation_metadata();

    }

}

void File_manager_class::open_log_output_file() {

    std::filesystem::path dir(std::filesystem::current_path() / "Sim_Results"); // Main results directory

    std::filesystem::path Log_File_Name;
    std::filesystem::path Log_File_full_path;
    std::filesystem::path file_extention(".txt");

    auto open_type = std::ios::trunc;

    Log_File_Name = "Photon_Log";
    Log_File_Name.replace_extension(file_extention);
    Log_File_full_path = dir / Log_File_Name;
     
    Log_Output_File.open(Log_File_full_path, open_type);

}

void File_manager_class::write_image_data_to_file(Results_type* s_Ray_results) {

    for (int Image_order = direct; Image_order <= ORDER_NUM - 1; Image_order += 1) {

        Image_Output_files[Image_order] << s_Ray_results->Image_Coords[x]
                                        << " "
                                        << s_Ray_results->Image_Coords[y]
                                        << " "
                                        << s_Ray_results->Redshift_NT[Image_order]
                                        << " "
                                        << s_Ray_results->Flux_NT[Image_order]
                                        << " "
                                        << s_Ray_results->Intensity[Image_order][I] * CGS_TO_JANSKY
                                        << " "
                                        << s_Ray_results->Intensity[Image_order][Q] * CGS_TO_JANSKY
                                        << " "
                                        << s_Ray_results->Intensity[Image_order][U] * CGS_TO_JANSKY
                                        << " "
                                        << s_Ray_results->Intensity[Image_order][V] * CGS_TO_JANSKY;

        if (p_Initial_Conditions->Simulation_mode == 2) {

            Image_Output_files[Image_order] << " "
                                            << s_Ray_results->Source_Coords[e_r][Image_order]
                                            << " "
                                            << s_Ray_results->Source_Coords[e_phi][Image_order]
                                            << " "
                                            << s_Ray_results->Photon_Momentum[e_r][Image_order]
                                            << " "
                                            << s_Ray_results->Photon_Momentum[e_theta][Image_order]
                                            << " "
                                            << s_Ray_results->Photon_Momentum[e_phi][Image_order]
                                            << " ";

                switch (this->p_Initial_Conditions->Metric_params.e_Spacetime) {

                case Kerr:
                    Image_Output_files[Image_order] << s_Ray_results->Parameters.Spin;
                    break;

                case Wormhole:
                    Image_Output_files[Image_order] << s_Ray_results->Parameters.Spin
                                                    << " " 
                                                    << s_Ray_results->Parameters.Redshift_Parameter;
                    break;

                case Reg_Black_Hole:
                    Image_Output_files[Image_order] << s_Ray_results->Parameters.RBH_Parameter;
         
                    break;

                case Janis_Newman_Winicour:
                    Image_Output_files[Image_order] << s_Ray_results->Parameters.JNW_Gamma_Parameter;
                    break;

                case Einstein_Gauss_Bonnet:
                    Image_Output_files[Image_order] << s_Ray_results->Parameters.GB_Gamma_Parameter;
                    break;

                case BH_w_Dark_Matter:
                    Image_Output_files[Image_order] << s_Ray_results->Parameters.Halo_Mass << " " << s_Ray_results->Parameters.Compactness;
                    break;
                }
        }

        Image_Output_files[Image_order] << '\n';

    }
}

void File_manager_class::log_photon_path(Results_type* s_Ray_results) {

    for (int log_index = 0; log_index <= s_Ray_results->Ray_log_struct.Log_length - 1; log_index++) {

        for (int state_index = 0; state_index <= e_State_Number - 1; state_index++) {

            Log_Output_File << s_Ray_results->Ray_log_struct.Ray_path_log[state_index + log_index * e_State_Number] << " ";
          
        }

        for (int stokes_index = I; stokes_index <= STOKES_PARAM_NUM - 1; stokes_index++) {

            Log_Output_File << s_Ray_results->Ray_log_struct.Ray_emission_log[stokes_index][0 + 2 * log_index] << " ";

        }

        Log_Output_File << '\n';
    }
   
};

void File_manager_class::close_image_output_files() {

    for (int File_Index = 0; File_Index <= ORDER_NUM - 1; File_Index++) {

        Image_Output_files[File_Index].close();
    }
}

void File_manager_class::close_log_output_file() {

    Log_Output_File.close();

}