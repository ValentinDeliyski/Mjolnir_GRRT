#include "Console_printing.h"

void Console_Printer_class::print_ASCII_art() {

        std::cout << " __       __                    __            __                   ______   _______   _______   ________   \n"
                  << "/  \\     /  |                  /  |          /  |                 /      \\ /       \\ /      \\ /        |  \n"
                  << "$$  \\   /$$ |     __   ______  $$ | _______  $$/   ______        /$$$$$$  |$$$$$$$  |$$$$$$$  |$$$$$$$$/   \n"
                  << "$$$  \\ /$$$ |    /  | /      \\ $$ |/       \\ /  | /      \\       $$ | _$$/ $$ |__$$ |$$ |__$$ |   $$ |     \n"
                  << "$$$$  /$$$$ |    $$/ /$$$$$$  |$$ |$$$$$$$  |$$ |/$$$$$$  |      $$ |/    |$$    $$< $$    $$<    $$ |     \n"
                  << "$$ $$ $$/$$ |    /  |$$ |  $$ |$$ |$$ |  $$ |$$ |$$ |  $$/       $$ |$$$$ |$$$$$$$  |$$$$$$$  |   $$ |     \n"
                  << "$$ |$$$/ $$ |    $$ |$$ \\__$$ |$$ |$$ |  $$ |$$ |$$ |            $$ \\__$$ |$$ |  $$ |$$ |  $$ |   $$ |     \n"
                  << "$$ | $/  $$ |    $$ |$$    $$/ $$ |$$ |  $$ |$$ |$$ |            $$    $$/ $$ |  $$ |$$ |  $$ |   $$ |     \n"
                  << "$$/      $$/__   $$ | $$$$$$/  $$/ $$/   $$/ $$/ $$/              $$$$$$/  $$/   $$/ $$/   $$/    $$/      \n"
                  << "           /  \\__$$ |                                                                                      \n"
                  << "           $$    $$/                                                                                       \n"
                  << "            $$$$$$/                                                                                        \n";

        std::cout << '\n';

}

void Console_Printer_class::print_sim_parameters(Initial_conditions_type* p_Initial_Conditions) {


    std::cout << "============================================================ SIMULATION METADATA ============================================================"
        << "\n"
        << "Spacetime: "
        << this->Metric_strings[p_Initial_Conditions->Metric_params.e_Spacetime]
        << "\n";

    Metric_parameters_type& Parameters = p_Initial_Conditions->Metric_params;

    switch (p_Initial_Conditions->Metric_params.e_Spacetime) {

    case Kerr:

        std::cout << "Spin Parameter [M]: " << Parameters.Spin << '\n';
        break;

    case Wormhole:

        std::cout << "Spin Parameter [M]: " << Parameters.Spin << '\n'
                  << "Redshift Parameter [-]: " << Parameters.Redshift_Parameter << '\n';
        break;

    case Reg_Black_Hole:

        std::cout << "Parameter [M]: " << Parameters.RBH_Parameter << '\n';
        break;

    case Janis_Newman_Winicour:

        std::cout << "Gamma [-]: " << Parameters.JNW_Gamma_Parameter << '\n';
        break;

    case Einstein_Gauss_Bonnet:

        std::cout << "Gamma [M^2]: " << Parameters.GB_Gamma_Parameter << '\n';
        break;

    case BH_w_Dark_Matter:

        std::cout << "Halo Mass [M]: " << Parameters.Halo_Mass << '\n'
                  << "Halo Compactness [-]: " << Parameters.Compactness << '\n';
        break;
    }

    std::cout << "Active Simulation Mode: " << p_Initial_Conditions->Simulation_mode << '\n';

    std::cout << "------------------------------------------------------- Observer Parameters -------------------------------------------------------" << "\n"
        << "Observer Distance [M]: " << p_Initial_Conditions->Observer_params.distance << '\n'
        << "Observer Inclination [Deg]: " << p_Initial_Conditions->Observer_params.inclination * 180.0 / M_PI << '\n'
        << "Observer Azimuth [Deg]: " << p_Initial_Conditions->Observer_params.azimuth * 180.0 / M_PI << '\n'
        << "Observation Frequency [Hz]: " << p_Initial_Conditions->Observer_params.obs_frequency << '\n';

    switch (p_Initial_Conditions->Simulation_mode) {

    case 1:

        std::cout << "Observation Window Dimentions (-X,+X,-Y,+Y) [M]: "
            <<  p_Initial_Conditions->Observer_params.x_min << ","
            <<  p_Initial_Conditions->Observer_params.x_max << ","
            <<  p_Initial_Conditions->Observer_params.y_min << ","
            <<  p_Initial_Conditions->Observer_params.y_max
            << '\n'
            << "Simulation Resolutoin: "
            <<  p_Initial_Conditions->Observer_params.resolution_x
            << " x "
            <<  p_Initial_Conditions->Observer_params.resolution_y
            << '\n';
        break;

    case 2:

        std::cout << "Number Of Photons Per Param Value: " << "" << '\n'
                  << "Number Of Param Values: " << p_Initial_Conditions->Sim_mode_2_param_value_number << '\n';

        break;

    default:

        std::cout << "Unsupported simulation mode!" << "\n";

    }


    std::cout << "------------------------------------------------------- Accretion Disk Parameters -------------------------------------------------------"
        << "\n"
        << "--------------------------- Density Model Parameters"
        << "\n";

    /*

    --------------------------------------- Print the accretion disk density parameters ---------------------------------------

    */

    switch ( p_Initial_Conditions->Disk_params.Density_profile_type) {

    case e_Power_law_profile:

        std::cout << "Density Profile: Power law"
            << "\n"
            << "Disk Opening Angle [tan(angle)]: "
            <<  p_Initial_Conditions->Disk_params.Power_law_disk_opening_angle
            << "\n"
            << "Density R_0 [M]: "
            <<  p_Initial_Conditions->Disk_params.Power_law_density_R_0
            << "\n"
            << "Density R_Cutoff [M]: "
            <<  p_Initial_Conditions->Disk_params.Power_law_density_R_cutoff
            << "\n"
            << "Density Cutoff Scale [M]: "
            <<  p_Initial_Conditions->Disk_params.Power_law_density_cutoff_scale
            << "\n";


        break;

    case e_Exponential_law_profile:

        std::cout << "Density Profile: Exponential law"
            << "\n"
            << "Density Height Scale [M]: "
            <<  p_Initial_Conditions->Disk_params.Exp_law_density_height_scale
            << "\n"
            << "Density Radial Scale [M]: "
            <<  p_Initial_Conditions->Disk_params.Exp_law_density_radial_scale
            << "\n";

        break;

    default:

        std::cout << "Unsupported Density Profile!" << "\n";

        break;

    }

    std::cout << "Maximum Density [g / cm^3]: "
        <<  p_Initial_Conditions->Disk_params.Electron_density_scale
        << "\n";

    /*

    --------------------------------------- Print the accretion disk temperature parameters ---------------------------------------

    */

    std::cout << "--------------------------- Temperature Model Parameters"
        << "\n";

    switch ( p_Initial_Conditions->Disk_params.Temperature_profile_type) {

    case e_Power_law_profile:

        std::cout << "Temperature Profile: Power law"
            << "\n"
            << "Temperature R_0 [M]: "
            <<  p_Initial_Conditions->Disk_params.Power_law_temperature_R_0
            << "\n"
            << "Temperature R_Cutoff [M]: "
            <<  p_Initial_Conditions->Disk_params.Power_law_temperature_R_cutoff << "\n"
            << "Temperature Cutoff Scale [M]: "
            <<  p_Initial_Conditions->Disk_params.Power_law_temperature_cutoff_scale << "\n";

        break;

    case e_Exponential_law_profile:

        std::cout << "Temperature Profile: Exponential law"
            << "\n"
            << "Temperature Height Scale [M]: "
            <<  p_Initial_Conditions->Disk_params.Exp_law_temperature_height_scale
            << "\n"
            << "Temperature Radial Scale [M]: "
            <<  p_Initial_Conditions->Disk_params.Exp_law_temperature_radial_scale
            << "\n";

        break;

    default:

        std::cout << "Unsupported Temperature Profile!" << "\n";

        break;
    }

    std::cout << "Maximum Temperature [K]: "
        <<  p_Initial_Conditions->Disk_params.Electron_temperature_scale
        << "\n";

    /*

    --------------------------------------- Print the disk ensamble parameters ---------------------------------------

    */

    std::cout << "--------------------------- Disk Synchrotron Emission Model Parameters"
        << "\n";

    switch ( p_Initial_Conditions->Disk_params.Ensamble_type) {

    case e_Phenomenological_ensamble:

        std::cout << "Disk Ensamble: Phenomenological"
            << "\n"
            << "Emission Power Law Exponent [-]: "
            <<  p_Initial_Conditions->Emission_params.Phenomenological_emission_power_law
            << "\n"
            << "Absorbtion Coefficient [?]: "
            <<  p_Initial_Conditions->Emission_params.Phenomenological_absorbtion_coeff
            << "\n"
            << "Source Function Power Law Exponent [-]: "
            <<  p_Initial_Conditions->Emission_params.Phenomenological_source_f_power_law
            << "\n"
            << "Emission Scale [erg / (cm^3 s sr Hz)]: "
            <<  p_Initial_Conditions->Emission_params.Phenomenological_emission_coeff
            << "\n";

        break;

    case e_Thermal_ensamble:

        std::cout << "Disk Ensamble: Thermal"
            << "\n";

        break;

    case e_Kappa_ensamble:
        std::cout << "Disk Ensamble: Kappa"
            << "\n"
            << "Kappa value [-]: "
            <<  p_Initial_Conditions->Emission_params.Kappa
            << "\n";

        break;

    default:

        std::cout << "Unsupported Ensamble!" << "\n";

        break;
    }

    std::cout << "Disk Magnetization [-]: "
        <<  p_Initial_Conditions->Disk_params.Magnetization << "\n"
        << "Disk Magnetic Field Geometry [-]: "
        << "[" <<  p_Initial_Conditions->Disk_params.Mag_field_geometry[0] << " "
        <<  p_Initial_Conditions->Disk_params.Mag_field_geometry[1] << " "
        <<  p_Initial_Conditions->Disk_params.Mag_field_geometry[2] << "]"
        << "\n";

    /*

    ================================================ Print the hotspot density parameters ================================================

    */

    std::cout << "------------------------------------------------------- Hotspot Parameters -------------------------------------------------------"
        << "\n"
        << "--------------------------- Density Model Parameters"
        << "\n";

    switch ( p_Initial_Conditions->Hotspot_params.Density_profile_type) {

    case e_Gaussian_profile:

        std::cout << "Density Profile: Gaussian"
            << "\n"
            << "Spread [M]: "
            <<  p_Initial_Conditions->Hotspot_params.Density_spread
            << "\n";

        break;

    case e_Spherical_profile:

        std::cout << "Density Profile: Spherical"
            << "\n"
            << "Radius [M]: "
            << p_Initial_Conditions->Hotspot_params.Radius
            << "\n";

        break;

    default:

        std::cout << "Unsupported Density Profile!" << "\n";

        break;
    }

    std::cout << "Maximum Density [g / cm^3]: "
        <<  p_Initial_Conditions->Hotspot_params.Electron_density_scale
        << "\n";

    /*

    --------------------------------------- Print the hotspot temperature parameters ---------------------------------------

    */

    std::cout << "--------------------------- Temperature Model Parameters"
        << "\n";

    switch ( p_Initial_Conditions->Hotspot_params.Temperature_profile_type) {

    case e_Gaussian_profile:

        std::cout << "Temperature Profile : Gaussian"
            << "\n"
            << "Spread [M]: "
            <<  p_Initial_Conditions->Hotspot_params.Temperature_spread
            << "\n";

        break;

    case e_Spherical_profile:

        std::cout << "Density Profile: Spherical"
            << "\n"
            << "Radius [M]: "
            << p_Initial_Conditions->Hotspot_params.Radius
            << "\n";
        break;

    default:

        std::cout << "Unsupported Temperature Profile!" << "\n";

        break;
    }

    std::cout << "Maximum Temperature [K]: "
        <<  p_Initial_Conditions->Hotspot_params.Electron_temperature_scale
        << "\n";

    /*

    --------------------------------------- Print the hotspot ensamble parameters ---------------------------------------

    */

    std::cout << "--------------------------- Hotspot Synchrotron Emission Model Parameters"
        << "\n";

    switch ( p_Initial_Conditions->Hotspot_params.Ensamble_type) {

    case e_Phenomenological_ensamble:

        std::cout << "Hotspot Ensamble: Phenomenological"
            << "\n"
            << "Emission Power Law Exponent [-]: "
            <<  p_Initial_Conditions->Emission_params.Phenomenological_emission_power_law
            << "\n"
            << "Absorbtion Coefficient [?]: "
            <<  p_Initial_Conditions->Emission_params.Phenomenological_absorbtion_coeff
            << "\n"
            << "Source Function Power Law Exponent [-]: "
            <<  p_Initial_Conditions->Emission_params.Phenomenological_source_f_power_law
            << "\n"
            << "Emission Scale [erg / (cm^3 s sr Hz)]: "
            <<  p_Initial_Conditions->Emission_params.Phenomenological_emission_coeff
            << "\n";

        break;

    case e_Thermal_ensamble:

        std::cout << "Hotspot Ensamble: Thermal"
            << "\n";

        break;

    case e_Kappa_ensamble:
        std::cout << "Hotspot Ensamble: Kappa"
            << "\n"
            << "Kappa value [-]: "
            <<  p_Initial_Conditions->Emission_params.Kappa
            << "\n";

        break;

    default:

        std::cout << "Unsupported Ensamble!" << "\n";

        break;
    }

    std::cout << "Hotspot Magnetization [-]: "
        <<  p_Initial_Conditions->Hotspot_params.Magnetization << "\n"
        << "Hotspot Magnetic Field Geometry [-]: "
        << "[" <<  p_Initial_Conditions->Hotspot_params.Mag_field_geometry[0] << " "
        <<  p_Initial_Conditions->Hotspot_params.Mag_field_geometry[1] << " "
        <<  p_Initial_Conditions->Hotspot_params.Mag_field_geometry[2] << "]"
        << "\n";


    std::cout << "--------------------------- Hotspot Position"
        << "\n";

    std::cout << "Distance [M]: "
        << p_Initial_Conditions->Hotspot_params.Position[e_r - 1]
        << "\n"
        << "Inclination [Deg]: "
        << p_Initial_Conditions->Hotspot_params.Position[e_theta - 1] * 180.0 / M_PI
        << "\n"
        << "Azimuth [Deg]: "
        << p_Initial_Conditions->Hotspot_params.Position[e_phi - 1] * 180.0 / M_PI
        << "\n";


    std::cout << "------------------------------------------------------- Novikov - Thorne Model Parameters -------------------------------------------------------"
        << "\n";

    /*

    --------------------------------------- Novikov - Thorne Disk Parameters ---------------------------------------

    */

    if ( p_Initial_Conditions->NT_params.evaluate_NT_disk) {

        std::cout << "Inner Disk Radius [M]: "
            <<  p_Initial_Conditions->NT_params.r_in
            << "\n"
            << "Outer Disk Radius [M]: "
            <<  p_Initial_Conditions->NT_params.r_out
            << "\n";
    }
    else {

        std::cout << "Novikov - Thorne Disk Evaluation Is Disabled.";

    }

    std::cout << '\n';

}



