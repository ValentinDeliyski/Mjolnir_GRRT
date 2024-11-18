#pragma once

#define _USE_MATH_DEFINES

#include "IO_files.h"
#include "Constants.h"
#include "Disk_models.h"
#include "Spacetimes.h"
#include <filesystem>
#include <iostream>

File_manager_class::File_manager_class(Initial_conditions_type *p_Initial_Conditions, bool truncate) {

    this->p_Initial_Conditions       = p_Initial_Conditions;
    this->Input_file_path_sim_mode_2 = p_Initial_Conditions->File_paths.Sim_mode_2_imput_path;
    this->Truncate_files             = truncate;
}

void File_manager_class::get_geodesic_data(double J_data[], double p_theta_data[], int* Data_number) {

    std::ifstream geodesic_data;
    std::string line;

    double J_input{};
    double P_input{};

    geodesic_data.open(Input_file_path_sim_mode_2, std::ios::in);

    while (true) {

        if (geodesic_data >> J_input >> P_input) {

            J_data[*Data_number] = J_input;
            p_theta_data[*Data_number] = P_input;
            
            *Data_number += 1;
        }

        if (geodesic_data.eof()) {

            break;
        }
    }

    geodesic_data.close();
}

void File_manager_class::write_simulation_metadata(int Sim_mode_2_number) {

    for (int Image_order = direct; Image_order <= ORDER_NUM - 1; Image_order += 1) {

        Image_Output_files[Image_order] << "============================================================ SIMULATION METADATA ============================================================"
                                        << "\n"
                                        << "Spacetime: "
                                        << this->Base_File_Names[this->p_Initial_Conditions->Metric_params.e_Spacetime]
                                        << "\n";

        Metric_parameters_type& Parameters = this->p_Initial_Conditions->Metric_params;

        switch (this->p_Initial_Conditions->Metric_params.e_Spacetime) {

        case Kerr:

            Image_Output_files[Image_order] << "Spin Parameter [M]: " << Parameters.Spin << '\n';
            break;

        case Wormhole:

            Image_Output_files[Image_order] << "Spin Parameter [M]: " << Parameters.Spin << "\n"
                                            << "Redshift Parameter [-]: " << Parameters.Redshift_Parameter << '\n';
            break;

        case Reg_Black_Hole:

            Image_Output_files[Image_order] << "Parameter [M]: " << Parameters.RBH_Parameter << '\n';
            break;

        case Janis_Newman_Winicour:

            Image_Output_files[Image_order] << "Gamma [-]: " << Parameters.JNW_Gamma_Parameter << '\n';
            break;

        case Einstein_Gauss_Bonnet:

            Image_Output_files[Image_order] << "Gamma [M^2]: " << Parameters.GB_Gamma_Parameter << '\n';
            break;

        case BH_w_Dark_Matter:

            Image_Output_files[Image_order] << "Halo Mass [M]: " << Parameters.Halo_Mass << '\n'
                                            << "Halo Compactness [-]: " << Parameters.Compactness << '\n';
            break;
        }

        

        Image_Output_files[Image_order] << "Active Simulation Mode: " << Active_Sim_Mode << '\n'
                                        << "Image Order [-]: " << Image_order << "\n";

     
        Image_Output_files[Image_order] << "------------------------------------------------------- Observer Parameters -------------------------------------------------------" << "\n"
                                        << "Observer Distance [M]: " << p_Initial_Conditions->Observer_params.distance << '\n'
                                        << "Observer Inclination [Deg]: " << p_Initial_Conditions->Observer_params.inclination * 180.0 / M_PI << '\n'
                                        << "Observer Azimuth [Deg]: " << p_Initial_Conditions->Observer_params.azimuth * 180.0 / M_PI << '\n'
                                        << "Observation Frequency [Hz]: " << p_Initial_Conditions->Observer_params.obs_frequency << '\n';

        switch (Active_Sim_Mode) {

        case 1:

            Image_Output_files[Image_order] << "Observation Window Dimentions (-X,+X,-Y,+Y) [M]: "
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

            Image_Output_files[Image_order] << "Number Of Photons Per Param Value: " << Sim_mode_2_number << '\n'
                                            << "Number Of Param Values: " << PARAM_SWEEP_NUMBER << '\n';

            break;

        default:

            Image_Output_files[Image_order] << "Unsupported simulation mode!" << "\n";

        }


        Image_Output_files[Image_order] << "------------------------------------------------------- Accretion Disk Parameters -------------------------------------------------------"
                                        << "\n"
                                        << "--------------------------- Density Model Parameters"
                                        << "\n";

        /*
        
        --------------------------------------- Print the accretion disk density parameters ---------------------------------------
        
        */

        switch (this->p_Initial_Conditions->Disk_params.Density_profile_type) {

        case e_Power_law_profile:

            Image_Output_files[Image_order] << "Density Profile: Power law"
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

            Image_Output_files[Image_order] << "Density Profile: Exponential law"
                                            << "\n"
                                            << "Density Height Scale [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Exp_law_density_height_scale
                                            << "\n"
                                            << "Density Radial Scale [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Exp_law_density_radial_scale
                                            << "\n";

            break;

        default:

            Image_Output_files[Image_order] << "Unsupported Density Profile!" << "\n";

            break;

        }

        Image_Output_files[Image_order] << "Maximum Density [g / cm^3]: "
                                        << this->p_Initial_Conditions->Disk_params.Electron_density_scale
                                        << "\n";
   
        /*
        
        --------------------------------------- Print the accretion disk temperature parameters ---------------------------------------
        
        */

        Image_Output_files[Image_order] << "--------------------------- Temperature Model Parameters"
                                        << "\n";

        switch (this->p_Initial_Conditions->Disk_params.Temperature_profile_type) {

        case e_Power_law_profile:

            Image_Output_files[Image_order] << "Temperature Profile: Power law"
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

            Image_Output_files[Image_order] << "Temperature Profile: Exponential law"
                                            << "\n"
                                            << "Temperature Height Scale [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Exp_law_temperature_height_scale
                                            << "\n"
                                            << "Temperature Radial Scale [M]: "
                                            << this->p_Initial_Conditions->Disk_params.Exp_law_temperature_radial_scale
                                            << "\n";

            break;

        default:

            Image_Output_files[Image_order] << "Unsupported Temperature Profile!" << "\n";

            break;
        }

        Image_Output_files[Image_order] << "Maximum Temperature [K]: "
                                        << this->p_Initial_Conditions->Disk_params.Electron_temperature_scale
                                        << "\n";

        /*
        
        --------------------------------------- Print the disk ensamble parameters ---------------------------------------
        
        */

        Image_Output_files[Image_order] << "--------------------------- Disk Synchrotron Emission Model Parameters"
                                        << "\n";

        switch (this->p_Initial_Conditions->Disk_params.Ensamble_type) {

        case e_Phenomenological_ensamble:

            Image_Output_files[Image_order] << "Disk Ensamble: Phenomenological"
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

            Image_Output_files[Image_order] << "Disk Ensamble: Thermal"
                                            << "\n";

            break;

        case e_Kappa_ensamble:
            Image_Output_files[Image_order] << "Disk Ensamble: Kappa"
                                            << "\n"
                                            << "Kappa value [-]: "
                                            << this->p_Initial_Conditions->Emission_params.Kappa
                                            << "\n";

            break;

        default:

            Image_Output_files[Image_order] << "Unsupported Ensamble!" << "\n";

            break;
        }

        Image_Output_files[Image_order] << "Disk Magnetization [-]: "
                                        << this->p_Initial_Conditions->Disk_params.Magnetization << "\n"
                                        << "Disk Magnetic Field Geometry [-]: "
                                        << "[" << this->p_Initial_Conditions->Disk_params.Mag_field_geometry[0] << " "
                                        << this->p_Initial_Conditions->Disk_params.Mag_field_geometry[1] << " "
                                        << this->p_Initial_Conditions->Disk_params.Mag_field_geometry[2] << "]"
                                        << "\n";

        /*
        
        ================================================ Print the hotspot density parameters ================================================
        
        */

        Image_Output_files[Image_order] << "------------------------------------------------------- Hotspot Parameters -------------------------------------------------------"
                                        << "\n" 
                                        << "--------------------------- Density Model Parameters"
                                        << "\n";

        switch (this->p_Initial_Conditions->Hotspot_params.Density_profile_type) {

        case e_Gaussian_profile:

            Image_Output_files[Image_order] << "Density Profile: Gaussian"
                                            << "\n"
                                            << "Spread [M]: "
                                            << this->p_Initial_Conditions->Hotspot_params.Density_spread
                                            << "\n";
                                            
            break;

        case e_Spherical_profile:

            Image_Output_files[Image_order] << "Density Profile: Spherical"
                                            << "\n"
                                            << "Radius [M]: "
                                            << this->p_Initial_Conditions->Hotspot_params.Radius
                                            << "\n";            

            break;

        default:

            Image_Output_files[Image_order] << "Unsupported Density Profile!" << "\n";

            break;
        }

        Image_Output_files[Image_order] << "Maximum Density [g / cm^3]: "
                                        << this->p_Initial_Conditions->Hotspot_params.Electron_density_scale
                                        << "\n";

        /*

        --------------------------------------- Print the hotspot temperature parameters ---------------------------------------

        */

        Image_Output_files[Image_order] << "--------------------------- Temperature Model Parameters"
                                        << "\n";

        switch (this->p_Initial_Conditions->Hotspot_params.Temperature_profile_type) {

        case e_Gaussian_profile:

            Image_Output_files[Image_order] << "Temperature Profile : Gaussian"
                                            << "\n"
                                            << "Spread [M]: "
                                            << this->p_Initial_Conditions->Hotspot_params.Temperature_spread
                                            << "\n";
                                            
            break;

        case e_Spherical_profile:

            Image_Output_files[Image_order] << "Temperature Profile: Spherical"
                                            << "\n"
                                            << "Radius [M]: "
                                            << this->p_Initial_Conditions->Hotspot_params.Radius
                                            << "\n";

            break;

        default:

            Image_Output_files[Image_order] << "Unsupported Temperature Profile!" << "\n";

            break;
        }

        Image_Output_files[Image_order] << "Maximum Temperature [K]: "
                                        << this->p_Initial_Conditions->Hotspot_params.Electron_temperature_scale
                                        << "\n";

        /*
        
        --------------------------------------- Print the hotspot ensamble parameters ---------------------------------------
        
        */

        Image_Output_files[Image_order] << "--------------------------- Hotspot Synchrotron Emission Model Parameters"
                                        << "\n";

        switch (this->p_Initial_Conditions->Hotspot_params.Ensamble_type) {

        case e_Phenomenological_ensamble:

            Image_Output_files[Image_order] << "Hotspot Ensamble: Phenomenological"
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

            Image_Output_files[Image_order] << "Hotspot Ensamble: Thermal"
                                            << "\n";

            break;

        case e_Kappa_ensamble:
            Image_Output_files[Image_order] << "Hotspot Ensamble: Kappa"
                                            << "\n"
                                            << "Kappa value [-]: "
                                            << this->p_Initial_Conditions->Emission_params.Kappa
                                            << "\n";

            break;

        default:

            Image_Output_files[Image_order] << "Unsupported Ensamble!" << "\n";

            break;
        }

        Image_Output_files[Image_order] << "Hotspot Magnetization [-]: "
                                        << this->p_Initial_Conditions->Disk_params.Magnetization << "\n"
                                        << "Hotspot Magnetic Field Geometry [-]: "
                                        << "[" << this->p_Initial_Conditions->Disk_params.Mag_field_geometry[0] << " "
                                        << this->p_Initial_Conditions->Disk_params.Mag_field_geometry[1] << " "
                                        << this->p_Initial_Conditions->Disk_params.Mag_field_geometry[2] << "]"
                                        << "\n";

        Image_Output_files[Image_order] << "--------------------------- Hotspot Position"
                                        << "\n";

        Image_Output_files[Image_order] << "Distance [M]: "
                                        << this->p_Initial_Conditions->Hotspot_params.Position[e_r - 1]
                                        << "\n"
                                        << "Inclination [Deg]: "
                                        << this->p_Initial_Conditions->Hotspot_params.Position[e_theta - 1] * 180.0 / M_PI
                                        << "\n"
                                        << "Azimuth [Deg]: "
                                        << this->p_Initial_Conditions->Hotspot_params.Position[e_phi - 1] * 180.0 / M_PI
                                        << "\n";


        Image_Output_files[Image_order] << "------------------------------------------------------- Novikov - Thorne Model Parameters -------------------------------------------------------"
                                        << "\n";

        /*
        
        --------------------------------------- Novikov - Thorne Disk Parameters ---------------------------------------
        
        */

        if (this->p_Initial_Conditions->NT_params.evaluate_NT_disk){

            Image_Output_files[Image_order] << "Inner Disk Radius [M]: "
                                            << this->p_Initial_Conditions->NT_params.r_in
                                            << "\n"
                                            << "Outer Disk Radius [M]: "
                                            << this->p_Initial_Conditions->NT_params.r_out
                                            << "\n";
        }
        else {

            Image_Output_files[Image_order] << "Novikov - Thorne Disk Evaluation Is Disabled." << "\n";

        }

        Image_Output_files[Image_order] << "------------------------------------------------------- Simulation Results -------------------------------------------------------"
                                        << "\n";


        Image_Output_files[Image_order] << "Image X Coord [M],"
                                        << " "
                                        << "Image Y Coord [M],"
                                        << " "
                                        << "Novikov-Thorne Disk Redshift [-],"
                                        << " "
                                        << "Novikov-Thorne Flux [M^-2],"
                                        << " "
                                        << "Synchotron Intensity I [Jy/sRad],"
                                        << " "
                                        << "Synchotron Intensity Q [Jy/sRad],"
                                        << " "
                                        << "Synchotron Intensity U [Jy/sRad],"
                                        << " "
                                        << "Synchotron Intensity V [Jy/sRad]";

            if (Active_Sim_Mode == 2) {

                Image_Output_files[Image_order] << ", Source r Coord [M],"
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

                    Image_Output_files[Image_order] << "Spin Parameter";
                    break;

                case Wormhole:

                    Image_Output_files[Image_order] << "Spin Parameter,"
                        << " "
                        << "Redshift Parameter";
                    break;

                case Reg_Black_Hole:

                    Image_Output_files[Image_order] << "Parameter";
                    break;

                case Janis_Newman_Winicour:

                    Image_Output_files[Image_order] << "Gamma";
                    break;

                case Einstein_Gauss_Bonnet:

                    Image_Output_files[Image_order] << "Gamma";
                    break;

                case BH_w_Dark_Matter:

                    Image_Output_files[Image_order] << "Halo Mass,"
                                                    << " " 
                                                    << "Halo Compactness";
                    break;
                }
            }

            Image_Output_files[Image_order] << '\n';
    }
}

void File_manager_class::open_image_output_files(int Sim_mode_2_number) {

    const int File_number = SPACETIME_NUMBER;

    // Create the path to the main results directory

    std::string Output_directory_path = this->p_Initial_Conditions->File_paths.Output_file_directory + "/" + this->p_Initial_Conditions->File_paths.Simulation_name;
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
    std::filesystem::path Image_file_names[File_number];
    std::filesystem::path Momenta_file_names[File_number];

    // Init the std::path variables of the full file paths
    std::filesystem::path Image_full_path[ORDER_NUM]{};
    std::filesystem::path Momenta_full_path[ORDER_NUM]{};

    // Specify the output file extention
    std::filesystem::path file_extention(".txt");

    // Set weather we truncate the file upon opening or not
    auto open_type = std::ios::app;

    if (Truncate_files) {

        open_type = std::ios::trunc;

    }

    // Loop over all the files and populate the (so far empty) 
    
    for (int File_Index = 0; File_Index <= ORDER_NUM - 1; File_Index += 1) {

        if (0 == strcmp(static_cast<const char*>(this->p_Initial_Conditions->File_paths.Common_file_names.c_str()), "")) {

            Image_file_names[File_Index] = this->Base_File_Names[this->p_Initial_Conditions->Metric_params.e_Spacetime]
                                         + "_n"
                                         + std::to_string(File_Index);

        }
        else {

            Image_file_names[File_Index] = this->p_Initial_Conditions->File_paths.Common_file_names
                                         + "_n"
                                         + std::to_string(File_Index);

        }

        Image_file_names[File_Index].replace_extension(file_extention);
        Image_full_path[File_Index] = dir / Image_file_names[File_Index];
       

        Image_Output_files[File_Index].open(Image_full_path[File_Index], open_type);

    }


    if (this->Truncate_files) {

        // If we are truncating the file, we should write the metadata to it
        File_manager_class::write_simulation_metadata(Sim_mode_2_number);

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

    File_manager_class::write_simulation_metadata(int(0));
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

        if (Active_Sim_Mode == 2) {

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

    //for (int log_index = 0; log_index <= s_Ray_results->Ray_log.size() - 1; log_index++) {

    //    for (int state_index = e_r; state_index <= e_State_Number - 1; state_index++) {

    //        Log_Output_File << s_Ray_results->Ray_log[log_index][state_index] << " ";
    //      
    //    }

    //    Log_Output_File << '\n';
    //}
   
};

void File_manager_class::close_image_output_files() {

    for (int File_Index = 0; File_Index <= ORDER_NUM - 1; File_Index++) {

        Image_Output_files[File_Index].close();
    }
}

void File_manager_class::close_log_output_file() {

    Log_Output_File.close();

}