#include "Input_parser.h"
#include "tinyxml2.h"

Return_Values static parse_hotspot_params(tinyxml2::XMLElement* Hotspot_element, Hotspot_model_parameters_type *Hotspot_params) {

    tinyxml2::XMLElement* temp_param_var;

    // -------------------- The ensamble type

    temp_param_var = Hotspot_element->FirstChildElement("Ensamble_type");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the ensamble type!" << "\n"; return ERROR; }
    std::string Ensamble_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Kappa")) {

        Hotspot_params->Ensamble_type = e_Kappa_ensamble;

    }
    else if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Thermal")) {

        Hotspot_params->Ensamble_type = e_Thermal_ensamble;

    }
    else if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Phenomenological")) {

        Hotspot_params->Ensamble_type = e_Phenomenological_ensamble;

    }
    else { std::cout << "Unsupported ensamble type for the hotspot!" << "\n"; return ERROR; }

    
    // -------------------- The density profile

    temp_param_var = Hotspot_element->FirstChildElement("Density_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the density profile type!" << "\n"; return ERROR; }
    std::string Density_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Density_type_string.c_str()), "Gaussian")) {

        Hotspot_params->Density_profile_type = e_Gaussian_profile;

    }
    else { std::cout << "Unsupported density profile type for the hotspot!" << "\n"; return ERROR; }

    // -------------------- The temperature profile

    temp_param_var = Hotspot_element->FirstChildElement("Temperature_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the temperature profile type!" << "\n"; return ERROR; }
    std::string Temperature_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Temperature_type_string.c_str()), "Gaussian")) {

        Hotspot_params->Temperature_profile_type = e_Gaussian_profile;

    }
    else { std::cout << "Unsupported temperature profile type for the hotspot!" << "\n"; return ERROR; }


    // -------------------- The density sclae factor
    temp_param_var = Hotspot_element->FirstChildElement("Density_scale_factor");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the density scale factor!" << "\n"; return ERROR; }
    Hotspot_params->Electron_density_scale = std::stod(temp_param_var->GetText());

    // -------------------- The density sclae factor
    temp_param_var = Hotspot_element->FirstChildElement("Temperature_scale_factor");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the temperature scale factor!" << "\n"; return ERROR; }
    Hotspot_params->Electron_temperature_scale = std::stod(temp_param_var->GetText());

    // -------------------- The magnetic field geometry X component
    temp_param_var = Hotspot_element->FirstChildElement("Mag_field_geometry_X");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the magnetic field geometry X component!" << "\n"; return ERROR; }
    Hotspot_params->Mag_field_geometry[0] = std::stod(temp_param_var->GetText());

    // -------------------- The magnetic field geometry Y component
    temp_param_var = Hotspot_element->FirstChildElement("Mag_field_geometry_Y");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the magnetic field geometry Y component!" << "\n"; return ERROR; }
    Hotspot_params->Mag_field_geometry[1] = std::stod(temp_param_var->GetText());

    // -------------------- The magnetic field geometry Z component
    temp_param_var = Hotspot_element->FirstChildElement("Mag_field_geometry_Z");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the magnetic field geometry Z component!" << "\n"; return ERROR; }
    Hotspot_params->Mag_field_geometry[2] = std::stod(temp_param_var->GetText());

    // -------------------- The distance to the hotspot center
    temp_param_var = Hotspot_element->FirstChildElement("Distance");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the distance to the hotspot center!" << "\n"; return ERROR; }
    Hotspot_params->Position[e_r] = std::stod(temp_param_var->GetText());

    // -------------------- The hotspot inclination
    temp_param_var = Hotspot_element->FirstChildElement("Inclination");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot inclination!" << "\n"; return ERROR; }
    Hotspot_params->Position[e_theta] = std::stod(temp_param_var->GetText());

    // -------------------- The hotspot azimuth
    temp_param_var = Hotspot_element->FirstChildElement("Azimuth");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot Azimuth!" << "\n"; return ERROR; }
    Hotspot_params->Position[e_phi] = std::stod(temp_param_var->GetText());

    // -------------------- The magnetization
    temp_param_var = Hotspot_element->FirstChildElement("Magnetization");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the magnetization!" << "\n"; return ERROR; }
    Hotspot_params->Magnetization = std::stod(temp_param_var->GetText());

    /* ================================== The gaussian profile parameters ================================== */

    tinyxml2::XMLElement* Gaussian_node = Hotspot_element->FirstChildElement("Gaussian_profile");
    if (Gaussian_node == nullptr) { std::cout << "Failed to parse the ensamble type!" << "\n"; return ERROR; }

    if (e_Gaussian_profile == Hotspot_params->Density_profile_type) {

        // -------------------- The magnetization
        temp_param_var = Gaussian_node->FirstChildElement("Density_spread");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the density Gaussian spread!" << "\n"; return ERROR; }
        Hotspot_params->Density_spread = std::stod(temp_param_var->GetText());


    }

    if (e_Gaussian_profile == Hotspot_params->Temperature_profile_type) {

        // -------------------- The magnetization
        temp_param_var = Gaussian_node->FirstChildElement("Temperature_spread");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the temperature Gaussian spread!" << "\n"; return ERROR; }
        Hotspot_params->Temperature_spread = std::stod(temp_param_var->GetText());


    }

    return OK;

}

Return_Values static parse_disk_params(tinyxml2::XMLElement* Accretion_disk_element, Disk_model_parameters_type* Disk_params) {

    tinyxml2::XMLElement* temp_param_var;

    /* ======================================== Common paramaters ======================================== */

    tinyxml2::XMLElement* Common_paramaters_element = Accretion_disk_element->FirstChildElement("Common_parameters");
    if (Common_paramaters_element == nullptr) { std::cout << "Failed to find the Common parameters element!" << "\n"; return ERROR; }

    // -------------------- The ensamble type

    temp_param_var = Common_paramaters_element->FirstChildElement("Ensamble_type");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the ensamble type!" << "\n"; return ERROR; }
    std::string Ensamble_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Kappa")) {

        Disk_params->Ensamble_type = e_Kappa_ensamble;

    }
    else if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Thermal")) {

        Disk_params->Ensamble_type = e_Thermal_ensamble;

    }
    else if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Phenomenological")) {

        Disk_params->Ensamble_type = e_Phenomenological_ensamble;

    }
    else { std::cout << "Unsupported ensamble type for the disk!" << "\n"; return ERROR; }

    // -------------------- The density profile

    temp_param_var = Common_paramaters_element->FirstChildElement("Density_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the density profile type!" << "\n"; return ERROR; }
    std::string Profile_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Profile_type_string.c_str()), "Power Law")) {

        Disk_params->Density_profile_type = e_Power_law_profile;

    }
    else if (0 == strcmp(static_cast<const char*>(Profile_type_string.c_str()), "Exponential Law")) {


        Disk_params->Density_profile_type = e_Exponential_law_profile;

    }
    else { std::cout << "Unsupported density profile type for the disk!" << "\n"; return ERROR; }

    // -------------------- The temperature profile

    temp_param_var = Common_paramaters_element->FirstChildElement("Temperature_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the temperature profile type!" << "\n"; return ERROR; }

    Profile_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Profile_type_string.c_str()), "Power Law")) {

        Disk_params->Temperature_profile_type = e_Power_law_profile;

    }
    else if (0 == strcmp(static_cast<const char*>(Profile_type_string.c_str()), "Exponential Law")) {


        Disk_params->Temperature_profile_type = e_Exponential_law_profile;

    }
    else { std::cout << "Unsupported temperature profile type for the disk!" << "\n"; return ERROR; }

    // -------------------- The disk electron density scale factor
    temp_param_var = Common_paramaters_element->FirstChildElement("Density_scale_factor");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk density scale factor!" << "\n"; return ERROR; }
    Disk_params->Electron_density_scale = std::stod(temp_param_var->GetText());

    // -------------------- The disk electron temperature scale factor
    temp_param_var = Common_paramaters_element->FirstChildElement("Temperature_scale_factor");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk temperature scale factor!" << "\n"; return ERROR; }
    Disk_params->Electron_temperature_scale = std::stod(temp_param_var->GetText());

    // -------------------- The disk magnetic field geometry X component
    temp_param_var = Common_paramaters_element->FirstChildElement("Mag_field_geometry_X");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field geometry!" << "\n"; return ERROR; }
    Disk_params->Mag_field_geometry[0] = std::stod(temp_param_var->GetText());

    // -------------------- The disk magnetic field geometry Y component
    temp_param_var = Common_paramaters_element->FirstChildElement("Mag_field_geometry_Y");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field geometry!" << "\n"; return ERROR; }
    Disk_params->Mag_field_geometry[1] = std::stod(temp_param_var->GetText());

    // -------------------- The disk magnetic field geometry Z component
    temp_param_var = Common_paramaters_element->FirstChildElement("Mag_field_geometry_Z");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field geometry!" << "\n"; return ERROR; }
    Disk_params->Mag_field_geometry[2] = std::stod(temp_param_var->GetText());

    // -------------------- The disk magnetization
    temp_param_var = Common_paramaters_element->FirstChildElement("Magnetization");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk Magnetization!" << "\n"; return ERROR; }
    Disk_params->Magnetization = std::stod(temp_param_var->GetText());

    /* ======================================== Disk profile paramaters ======================================== */

    tinyxml2::XMLElement* Power_law_profile_element = Accretion_disk_element->FirstChildElement("Power_law_profile");
    tinyxml2::XMLElement* Exponential_law_profile_element = Accretion_disk_element->FirstChildElement("Exponential_law_profile");

    if (e_Power_law_profile == Disk_params->Density_profile_type) {

        if (Power_law_profile_element == nullptr) { std::cout << "Failed to find the power law profile element!" << "\n"; return ERROR; }

        // -------------------- The density radial power law
        temp_param_var = Power_law_profile_element->FirstChildElement("Density_radial_power_law");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk density radial power law!" << "\n"; return ERROR; }
        Disk_params->Power_law_density_radial_power_law = std::stod(temp_param_var->GetText());

        // -------------------- The density radial cutoff scale
        temp_param_var = Power_law_profile_element->FirstChildElement("Density_cutoff_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk density radial cutoff scale!" << "\n"; return ERROR; }
        Disk_params->Power_law_density_cutoff_scale = std::stod(temp_param_var->GetText());

        // -------------------- The density r_cutoff
        temp_param_var = Power_law_profile_element->FirstChildElement("Density_r_cutoff");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk density r_cutoff!" << "\n"; return ERROR; }
        Disk_params->Power_law_density_R_cutoff = std::stod(temp_param_var->GetText());

        // -------------------- The density r_0
        temp_param_var = Power_law_profile_element->FirstChildElement("Density_r_0");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk density r_0!" << "\n"; return ERROR; }
        Disk_params->Power_law_density_R_0 = std::stod(temp_param_var->GetText());

        // -------------------- The disk opening angle
        temp_param_var = Power_law_profile_element->FirstChildElement("Opening_angle");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk opening angle!" << "\n"; return ERROR; }
        Disk_params->Power_law_disk_opening_angle = std::stod(temp_param_var->GetText());

    }
    else if (e_Exponential_law_profile == Disk_params->Density_profile_type) {

        if (Exponential_law_profile_element == nullptr) { std::cout << "Failed to find the exponential law profile element!" << "\n"; return ERROR; }

        // -------------------- The density height scale
        temp_param_var = Exponential_law_profile_element->FirstChildElement("Density_exp_height_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk density height scale!" << "\n"; return ERROR; }
        Disk_params->Exp_law_density_height_scale = std::stod(temp_param_var->GetText());

        // -------------------- The density radial scale
        temp_param_var = Exponential_law_profile_element->FirstChildElement("Density_exp_radial_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk density radial scale!" << "\n"; return ERROR; }
        Disk_params->Exp_law_density_radial_scale = std::stod(temp_param_var->GetText());

    }

    if (e_Power_law_profile == Disk_params->Temperature_profile_type) {

        // -------------------- The temperature radial power law
        temp_param_var = Power_law_profile_element->FirstChildElement("Temperature_radial_power_law");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk temperature radial power law!" << "\n"; return ERROR; }
        Disk_params->Power_law_temperature_radial_power_law = std::stod(temp_param_var->GetText());

        // -------------------- The temperature radial cutoff scale
        temp_param_var = Power_law_profile_element->FirstChildElement("Temperature_cutoff_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk temperature radial cutoff scale!" << "\n"; return ERROR; }
        Disk_params->Power_law_temperature_cutoff_scale = std::stod(temp_param_var->GetText());

        // -------------------- The temperature r_cutoff
        temp_param_var = Power_law_profile_element->FirstChildElement("Temperature_r_cutoff");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk temperature r_cutoff!" << "\n"; return ERROR; }
        Disk_params->Power_law_temperature_R_cutoff = std::stod(temp_param_var->GetText());

        // -------------------- The temperature r_0
        temp_param_var = Power_law_profile_element->FirstChildElement("Temperature_r_0");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk temperature r_0!" << "\n"; return ERROR; }
        Disk_params->Power_law_temperature_R_0 = std::stod(temp_param_var->GetText());

    }
    else if (e_Exponential_law_profile == Disk_params->Density_profile_type) {

        // -------------------- The temperature height scale
        temp_param_var = Exponential_law_profile_element->FirstChildElement("Temperature_exp_height_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk temperature height scale!" << "\n"; return ERROR; }
        Disk_params->Exp_law_temperature_height_scale = std::stod(temp_param_var->GetText());

        // -------------------- The temperature radial scale
        temp_param_var = Exponential_law_profile_element->FirstChildElement("Temperature_exp_radial_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk temperature radial scale!" << "\n"; return ERROR; }
        Disk_params->Exp_law_temperature_radial_scale = std::stod(temp_param_var->GetText());

    }
    
    return OK;

}

Return_Values static parse_integrator_params(tinyxml2::XMLElement* Integrator_element, Integrator_parameters_type* Integrator_params) {

    tinyxml2::XMLElement* temp_param_var;

    // -------------------- Init stepsize
    temp_param_var = Integrator_element->FirstChildElement("init_stepsize");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the initial step size!" << "\n"; return ERROR; }
    Integrator_params->Init_stepzie = std::stod(temp_param_var->GetText());

    // -------------------- RK45 accuracy
    temp_param_var = Integrator_element->FirstChildElement("RK45_accuracy");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the RK45 accuracy parameter!" << "\n"; return ERROR; }
    Integrator_params->RK_45_accuracy = std::stod(temp_param_var->GetText());

    // -------------------- Step controller safety 1
    temp_param_var = Integrator_element->FirstChildElement("step_controller_safety_factor_1");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the step controller safety parameter 1!" << "\n"; return ERROR; }
    Integrator_params->Safety_1 = std::stod(temp_param_var->GetText());

    // -------------------- Step controller safety 2
    temp_param_var = Integrator_element->FirstChildElement("step_controller_safety_factor_2");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the step controller safety parameter 2!" << "\n"; return ERROR; }
    Integrator_params->Safety_2 = std::stod(temp_param_var->GetText());

    // -------------------- Step controller I gain
    temp_param_var = Integrator_element->FirstChildElement("step_controller_I_gain");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the step controller I gain!" << "\n"; return ERROR; }
    Integrator_params->Gain_I = std::stod(temp_param_var->GetText());

    // -------------------- Step controller P gain
    temp_param_var = Integrator_element->FirstChildElement("step_controller_P_gain");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the step controller P gain!" << "\n"; return ERROR; }
    Integrator_params->Gain_P = std::stod(temp_param_var->GetText());

    // -------------------- Step controller D gain
    temp_param_var = Integrator_element->FirstChildElement("step_controller_D_gain");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the step controller D gain!" << "\n"; return ERROR; }
    Integrator_params->Gain_D = std::stod(temp_param_var->GetText());

    // -------------------- Max integration count 
    temp_param_var = Integrator_element->FirstChildElement("max_integration_count");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the maximum integration count!" << "\n"; return ERROR; }
    Integrator_params->Max_integration_count = std::stod(temp_param_var->GetText());

    return OK;

}

Return_Values static parse_emission_model_params(tinyxml2::XMLElement* Emission_model_element, Initial_conditions_type* p_Init_conditions) {

    std::ifstream input_params;
    input_params.open(input_file_path, std::ios::in);

    tinyxml2::XMLElement* temp_param_var;

    if ((e_Kappa_ensamble == p_Init_conditions->Disk_params.Ensamble_type &&    0 != p_Init_conditions->Disk_params.Electron_density_scale) ||
         e_Kappa_ensamble == p_Init_conditions->Hotspot_params.Ensamble_type && 0 != p_Init_conditions->Hotspot_params.Electron_density_scale) {

        // -------------------- Kappa
        temp_param_var = Emission_model_element->FirstChildElement("Kappa");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the kappa value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Kappa = std::stod(temp_param_var->GetText());

    }

    if (e_Phenomenological_ensamble == p_Init_conditions->Disk_params.Ensamble_type ||
        e_Phenomenological_ensamble == p_Init_conditions->Hotspot_params.Ensamble_type) {


        // -------------------- Emission power law
        temp_param_var = Emission_model_element->FirstChildElement("Emission_power_law");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the phenomenological emission power law!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Phenomenological_emission_power_law = std::stod(temp_param_var->GetText());

        // -------------------- Source function power law
        temp_param_var = Emission_model_element->FirstChildElement("Source_f_power_law");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the phenomenological source function power law!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Phenomenological_source_f_power_law = std::stod(temp_param_var->GetText());

        // -------------------- Absorbtion coefficient
        temp_param_var = Emission_model_element->FirstChildElement("Absorbtion_coeff");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the phenomenological absorbtion coefficient!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Phenomenological_absorbtion_coeff = std::stod(temp_param_var->GetText());

        // -------------------- Emission coefficient
        temp_param_var = Emission_model_element->FirstChildElement("Emission_coeff");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the phenomenological emission coefficient!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Phenomenological_emission_coeff = std::stod(temp_param_var->GetText());

    }

    return OK;

}

Return_Values static parse_observer_parameters(tinyxml2::XMLElement* Observer_element, Observer_parameters_type* Observer_params) {

    tinyxml2::XMLElement* temp_param_var;

    // -------------------- Distance
    temp_param_var = Observer_element->FirstChildElement("Distance");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observer distance!" << "\n"; return ERROR; }
    Observer_params->distance = std::stod(temp_param_var->GetText());

    // -------------------- Inclination
    temp_param_var = Observer_element->FirstChildElement("Inclination");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observer inclination!" << "\n"; return ERROR; }
    Observer_params->inclination = std::stod(temp_param_var->GetText());

    // -------------------- Azimuth
    temp_param_var = Observer_element->FirstChildElement("Azimuth");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observer inclination!" << "\n"; return ERROR; }
    Observer_params->azimuth = std::stod(temp_param_var->GetText());

    // -------------------- Camera Rotation Angle
    temp_param_var = Observer_element->FirstChildElement("Cam_rotation_angle");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the camera rotation angle!" << "\n"; return ERROR; }
    Observer_params->cam_rotation_angle = std::stod(temp_param_var->GetText());

    // -------------------- Image Y min
    temp_param_var = Observer_element->FirstChildElement("Image_y_min");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observation window Y min!" << "\n"; return ERROR; }
    Observer_params->y_min = std::stod(temp_param_var->GetText());

    // -------------------- Image Y max
    temp_param_var = Observer_element->FirstChildElement("Image_y_max");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observation window Y max!" << "\n"; return ERROR; }
    Observer_params->y_max = std::stod(temp_param_var->GetText());

    // -------------------- Image X min
    temp_param_var = Observer_element->FirstChildElement("Image_x_min");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observation window X min!" << "\n"; return ERROR; }
    Observer_params->x_min = std::stod(temp_param_var->GetText());

    // -------------------- Image X max
    temp_param_var = Observer_element->FirstChildElement("Image_x_max");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observation window X max!" << "\n"; return ERROR; }
    Observer_params->x_max = std::stod(temp_param_var->GetText());

    // -------------------- Image resolution Y
    temp_param_var = Observer_element->FirstChildElement("Resolution_y");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse Y resolution!" << "\n"; return ERROR; }
    Observer_params->resolution_y = std::stod(temp_param_var->GetText());

    // -------------------- Image resolution X
    temp_param_var = Observer_element->FirstChildElement("Resolution_x");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse X resolution!" << "\n"; return ERROR; }
    Observer_params->resolution_x = std::stod(temp_param_var->GetText());

    // -------------------- Observation frequency
    temp_param_var = Observer_element->FirstChildElement("Obs_frequency");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the observation frequency!" << "\n"; return ERROR; }
    Observer_params->obs_frequency = std::stod(temp_param_var->GetText());

    // -------------------- Polarization flag
    temp_param_var = Observer_element->FirstChildElement("Include_polarization");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the polarization flag!" << "\n"; return ERROR; }
    Observer_params->include_polarization = std::stoi(temp_param_var->GetText());

    return OK;

}

Return_Values static parse_NT_params(tinyxml2::XMLElement* NT_element, NT_parameters_type* NT_params) {

    tinyxml2::XMLElement* temp_param_var;

    // -------------------- Evaluation flag
    temp_param_var = NT_element->FirstChildElement("Evaluate_NT_disk");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the NT disk evaluation flag!" << "\n"; return ERROR; }
    NT_params->evaluate_NT_disk = std::stoi(temp_param_var->GetText());

    if (NT_params->evaluate_NT_disk) {

        // -------------------- Inner disk radius
        temp_param_var = NT_element->FirstChildElement("r_in");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the inner NT disk radius!" << "\n"; return ERROR; }
        NT_params->r_in = std::stoi(temp_param_var->GetText());

        // -------------------- Outer disk radius
        temp_param_var = NT_element->FirstChildElement("r_out");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the outer NT disk radius!" << "\n"; return ERROR; }
        NT_params->r_out = std::stoi(temp_param_var->GetText());

    }

    return OK;

}

Return_Values static parse_metric_parameters(tinyxml2::XMLElement* Metric_element, Metric_parameters_type* Metric_params) {

    tinyxml2::XMLElement* temp_param_var;

    std::string Metric_type = Metric_element->FirstChildElement("Metric_type")->GetText();

    if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Kerr")) {

        temp_param_var = Metric_element->FirstChildElement("Spin_parameter");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the metric spin parameter!" << "\n"; return ERROR; }
        Metric_params->Spin = std::stod(temp_param_var->GetText());

    }
    else if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Wormhole")) {

        temp_param_var = Metric_element->FirstChildElement("Spin_parameter");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the metric spin parameter!" << "\n"; return ERROR; }
        Metric_params->Spin = std::stod(temp_param_var->GetText());

        temp_param_var = Metric_element->FirstChildElement("WH_redshift");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the wormhole redshift parameter!" << "\n"; return ERROR; }
        Metric_params->Redshift_Parameter = std::stod(temp_param_var->GetText());

        temp_param_var = Metric_element->FirstChildElement("WH_r_throat");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the wormhole throat radius!" << "\n"; return ERROR; }
        Metric_params->R_throat = std::stod(temp_param_var->GetText());

        temp_param_var = Metric_element->FirstChildElement("WH_stop_at_throat");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the wormhole \"stop at the throat\" flag!" << "\n"; return ERROR; }
        Metric_params->Stop_At_Throat = std::stoi(temp_param_var->GetText());

    }
    else if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Janis-Newman-Winicour")) {

        temp_param_var = Metric_element->FirstChildElement("JNW_gamma");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Janis-Newman-Winicour metric parameter!" << "\n"; return ERROR; }
        Metric_params->JNW_Gamma_Parameter = std::stod(temp_param_var->GetText());

    }
    else if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Einstein-Gauss-Bonnet")) {

        temp_param_var = Metric_element->FirstChildElement("EGB_gamma");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Einstein-Gauss-Bonnet metric parameter!" << "\n"; return ERROR; }
        Metric_params->GB_Gamma_Parameter = std::stod(temp_param_var->GetText());

    }
    else if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Regular-Black-Hole")) {

        temp_param_var = Metric_element->FirstChildElement("RBH_param");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the regular black hole metric parameter!" << "\n"; return ERROR; }
        Metric_params->RBH_Parameter = std::stod(temp_param_var->GetText());

    }
    else if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Black-Hole-w-Dark-Matter")) {

        temp_param_var = Metric_element->FirstChildElement("Halo_compactness");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the dark matter halo compactness!" << "\n"; return ERROR; }
        Metric_params->Compactness = std::stod(temp_param_var->GetText());

        temp_param_var = Metric_element->FirstChildElement("Halo_mass");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the dark matter halo mass!" << "\n"; return ERROR; }
        Metric_params->Halo_Mass = std::stod(temp_param_var->GetText());

    }
    else { std::cout << "Unsupported metric type!" << "\n"; return ERROR; }

    return OK;

}

Return_Values parse_simulation_input_XML(std::string input_file_path, Initial_conditions_type* p_Initial_conditions) {

    tinyxml2::XMLDocument xml_doc;

    tinyxml2::XMLError e_parse_result = xml_doc.LoadFile(static_cast<const char*>(input_file_path.c_str()));
    if (e_parse_result != tinyxml2::XML_SUCCESS) { return ERROR; }

    tinyxml2::XMLNode* Root_node = xml_doc.FirstChildElement("Simulation_Input");
    if (Root_node == nullptr) { return ERROR; }

    /* ====================================== Parse the average pitch angle flag ====================================== */

    tinyxml2::XMLElement* Average_pitch_angle_flag = Root_node->FirstChildElement("Average_emission_pitch_angle");
    if (Average_pitch_angle_flag == nullptr) { std::cout << "Failed to find the Observer node!" << "\n"; return ERROR; }
    p_Initial_conditions->Average_electron_pitch_angle = std::stoi(Average_pitch_angle_flag->GetText());

    /* ====================================== Parse the observer parameters ====================================== */

    tinyxml2::XMLElement* Observer_element = Root_node->FirstChildElement("Observer");
    if (Observer_element == nullptr) { std::cout << "Failed to find the Observer node!" << "\n"; return ERROR; }
    if (OK != parse_observer_parameters(Observer_element, &p_Initial_conditions->Observer_params)) { return ERROR; }

    /* ====================================== Parse the metric parameters ====================================== */

    tinyxml2::XMLElement* Metric_element = Root_node->FirstChildElement("Metric");
    if (Metric_element == nullptr) { std::cout << "Failed to find the Metric node!" << "\n"; return ERROR; }
    if (OK != parse_metric_parameters(Metric_element, &p_Initial_conditions->Metric_params)) { return ERROR; };

    /* ====================================== Parse the integrator parameters ====================================== */

    tinyxml2::XMLElement* Integrator_element = Root_node->FirstChildElement("Integrator");
    if (Integrator_element == nullptr) { std::cout << "Failed to find the Metric node!" << "\n"; return ERROR; }
    if (OK != parse_integrator_params(Integrator_element, &p_Initial_conditions->Integrator_params)) { return ERROR; };

    /* ====================================== Parse the accretion disk parameters ====================================== */

    tinyxml2::XMLElement* Accretion_disk_element = Root_node->FirstChildElement("Accretion_Disk");
    if (Accretion_disk_element == nullptr) { std::cout << "Failed to find the Accretion Disk node!" << "\n"; return ERROR; }
    if (OK != parse_disk_params(Accretion_disk_element, &p_Initial_conditions->Disk_params)) { return ERROR; };

    /* ====================================== Parse the hotspot parameters ====================================== */

    tinyxml2::XMLElement* Hotspot_element = Root_node->FirstChildElement("Hotspot");
    if (Hotspot_element == nullptr) { std::cout << "Failed to find the Hotspot node!" << "\n"; return ERROR; }
    if (OK != parse_hotspot_params(Hotspot_element, &p_Initial_conditions->Hotspot_params)) { return ERROR; };

    /* ====================================== Parse the emission model parameters ====================================== */

    tinyxml2::XMLElement* Emission_model_element = Root_node->FirstChildElement("Emission_models");
    if (Emission_model_element == nullptr) { std::cout << "Failed to find the Hotspot node!" << "\n"; return ERROR; }
    if (OK != parse_emission_model_params(Emission_model_element, p_Initial_conditions)) { return ERROR; };

    /* ====================================== Parse the Novikov-Thorne model parameters ====================================== */

    tinyxml2::XMLElement* NT_element = Root_node->FirstChildElement("Novikov_Thorne_disk");
    if (NT_element == nullptr) { std::cout << "Failed to find the Hotspot node!" << "\n"; return ERROR; }
    if (OK != parse_NT_params(NT_element, &p_Initial_conditions->NT_params)) { return ERROR; };

    return OK;

}