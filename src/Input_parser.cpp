#include "Input_parser.h"
#include "tinyxml2.h"

Return_Values static parse_hotspot_params(tinyxml2::XMLElement* Hotspot_element, Hotspot_model_parameters_type *Hotspot_params) {

    tinyxml2::XMLElement* temp_param_var;

    // -------------------- The ensamble type

    temp_param_var = Hotspot_element->FirstChildElement("Ensamble_type");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot ensamble type!" << "\n"; return ERROR; }
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
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot density profile type!" << "\n"; return ERROR; }
    std::string Density_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Density_type_string.c_str()), "Gaussian")) {

        Hotspot_params->Density_profile_type = e_Gaussian_profile;

    }
    else if (0 == strcmp(static_cast<const char*>(Density_type_string.c_str()), "Sphere")) {

        Hotspot_params->Density_profile_type = e_Spherical_profile;

    }
    else { std::cout << "Unsupported density profile type for the hotspot!" << "\n"; return ERROR; }

    // -------------------- The temperature profile

    temp_param_var = Hotspot_element->FirstChildElement("Temperature_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot temperature profile type!" << "\n"; return ERROR; }
    std::string Temperature_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Temperature_type_string.c_str()), "Gaussian")) {

        Hotspot_params->Temperature_profile_type = e_Gaussian_profile;

    }
    else if (0 == strcmp(static_cast<const char*>(Temperature_type_string.c_str()), "Sphere")) {

        Hotspot_params->Temperature_profile_type = e_Spherical_profile;

    }
    else { std::cout << "Unsupported temperature profile type for the hotspot!" << "\n"; return ERROR; }

    // -------------------- The velocity profile

    temp_param_var = Hotspot_element->FirstChildElement("Velocity_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot velocity profile type!" << "\n"; return ERROR; }
    std::string Velocity_profile_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Velocity_profile_string.c_str()), "Theta Dependant")) {

        Hotspot_params->Velocity_profile_type = e_Theta_dependant;

    }
    else if (0 == strcmp(static_cast<const char*>(Velocity_profile_string.c_str()), "Keplarian")) {

        Hotspot_params->Velocity_profile_type = e_Keplarian;

    }else{std::cout << "Unsupported velocity profile type for the hotspot!" << "\n"; return ERROR; }

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
    Hotspot_params->Position[e_r - 1] = std::stod(temp_param_var->GetText());

    // -------------------- The hotspot inclination
    temp_param_var = Hotspot_element->FirstChildElement("Inclination");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot inclination!" << "\n"; return ERROR; }
    Hotspot_params->Position[e_theta - 1] = std::stod(temp_param_var->GetText());

    // -------------------- The hotspot azimuth
    temp_param_var = Hotspot_element->FirstChildElement("Azimuth");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot Azimuth!" << "\n"; return ERROR; }
    Hotspot_params->Position[e_phi - 1] = std::stod(temp_param_var->GetText());

    // -------------------- The magnetization
    temp_param_var = Hotspot_element->FirstChildElement("Magnetization");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the magnetization!" << "\n"; return ERROR; }
    Hotspot_params->Magnetization = std::stod(temp_param_var->GetText());

    /* ================================== The gaussian profile parameters ================================== */

    tinyxml2::XMLElement* Gaussian_node = Hotspot_element->FirstChildElement("Gaussian_profile");
    tinyxml2::XMLElement* Spherical_node = Hotspot_element->FirstChildElement("Spherical_profile");

    if (e_Gaussian_profile == Hotspot_params->Density_profile_type) {

        if (Gaussian_node == nullptr) { std::cout << "Failed to parse the ensamble type!" << "\n"; return ERROR; }

        // -------------------- The Gaussian density standard deviation
        temp_param_var = Gaussian_node->FirstChildElement("Density_spread");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the density Gaussian spread!" << "\n"; return ERROR; }
        Hotspot_params->Density_spread = std::stod(temp_param_var->GetText());


    }

    if (e_Gaussian_profile == Hotspot_params->Temperature_profile_type) {

        if (Gaussian_node == nullptr) { std::cout << "Failed to parse the ensamble type!" << "\n"; return ERROR; }

        // -------------------- The Gaussian temperature standard devivation
        temp_param_var = Gaussian_node->FirstChildElement("Temperature_spread");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the temperature Gaussian spread!" << "\n"; return ERROR; }
        Hotspot_params->Temperature_spread = std::stod(temp_param_var->GetText());


    }

    if (e_Spherical_profile == Hotspot_params->Temperature_profile_type || e_Spherical_profile == Hotspot_params->Density_profile_type) {

        if (Spherical_node == nullptr) { std::cout << "Failed to parse the ensamble type!" << "\n"; return ERROR; }

        // -------------------- The Spherical radius
        temp_param_var = Spherical_node->FirstChildElement("Radius");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Spherical hotspot radius!" << "\n"; return ERROR; }
        Hotspot_params->Radius = std::stod(temp_param_var->GetText());

    }

    // -------------------- The Temporal spread
    temp_param_var = Hotspot_element->FirstChildElement("Temporal_spread");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the temporal spread!" << "\n"; return ERROR; }
    Hotspot_params->Temporal_spread = std::stod(temp_param_var->GetText());

    // -------------------- The Coordiante time at max emission
    temp_param_var = Hotspot_element->FirstChildElement("Coord_time_at_max");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the coordiante time at max emission!" << "\n"; return ERROR; }
    Hotspot_params->Coord_time_at_max = std::stod(temp_param_var->GetText());

    return OK;

}

Return_Values static parse_disk_params(tinyxml2::XMLElement* Accretion_disk_element, Disk_model_parameters_type* Disk_params) {

    tinyxml2::XMLElement* temp_param_var;

    /* ======================================== Common paramaters ======================================== */

    tinyxml2::XMLElement* Common_paramaters_element = Accretion_disk_element->FirstChildElement("Common_parameters");
    if (Common_paramaters_element == nullptr) { std::cout << "Failed to find the Common parameters element!" << "\n"; return ERROR; }

    // -------------------- The ensamble type

    temp_param_var = Common_paramaters_element->FirstChildElement("Ensamble_type");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk ensamble type!" << "\n"; return ERROR; }
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
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk density profile type!" << "\n"; return ERROR; }
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
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk temperature profile type!" << "\n"; return ERROR; }

    Profile_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Profile_type_string.c_str()), "Power Law")) {

        Disk_params->Temperature_profile_type = e_Power_law_profile;

    }
    else if (0 == strcmp(static_cast<const char*>(Profile_type_string.c_str()), "Exponential Law")) {


        Disk_params->Temperature_profile_type = e_Exponential_law_profile;

    }
    else { std::cout << "Unsupported temperature profile type for the disk!" << "\n"; return ERROR; }

    // -------------------- The velocity profile

    temp_param_var = Common_paramaters_element->FirstChildElement("Velocity_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk velocity profile type!" << "\n"; return ERROR; }
    std::string Velocity_profile_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Velocity_profile_string.c_str()), "Theta Dependant")) {

        Disk_params->Velocity_profile_type = e_Theta_dependant;

    }
    else if (0 == strcmp(static_cast<const char*>(Velocity_profile_string.c_str()), "Keplarian")) {

        Disk_params->Velocity_profile_type = e_Keplarian;

    }
    else { std::cout << "Unsupported velocity profile type for the disk!" << "\n"; return ERROR; }


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

        if (Power_law_profile_element == nullptr) { std::cout << "Failed to find the power law profile element!" << "\n"; return ERROR; }

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

        if (Exponential_law_profile_element == nullptr) { std::cout << "Failed to find the exponential law profile element!" << "\n"; return ERROR; }

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

    // -------------------- PID controller I gain
    temp_param_var = Integrator_element->FirstChildElement("PID_controller_I_gain");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the PID controller I gain!" << "\n"; return ERROR; }
    Integrator_params->PID_gain_I = std::stod(temp_param_var->GetText());

    // -------------------- PID controller P gain
    temp_param_var = Integrator_element->FirstChildElement("PID_controller_P_gain");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the PID controller P gain!" << "\n"; return ERROR; }
    Integrator_params->PID_gain_P = std::stod(temp_param_var->GetText());

    // -------------------- PID controller D gain
    temp_param_var = Integrator_element->FirstChildElement("PID_controller_D_gain");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the PID controller D gain!" << "\n"; return ERROR; }
    Integrator_params->PID_gain_D = std::stod(temp_param_var->GetText());

    // -------------------- Gustafsson controller k_1 gain
    temp_param_var = Integrator_element->FirstChildElement("Gustafsson_controller_k_1");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the Gustafsson controller k_1 gain!" << "\n"; return ERROR; }
    Integrator_params->Gustafsson_k1 = std::stod(temp_param_var->GetText());

    // -------------------- Gustafsson controller k_2 gain
    temp_param_var = Integrator_element->FirstChildElement("Gustafsson_controller_k_2");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the Gustafsson controller k_2 gain!" << "\n"; return ERROR; }
    Integrator_params->Gustafsson_k2 = std::stod(temp_param_var->GetText());

    // -------------------- Max relative step increase
    temp_param_var = Integrator_element->FirstChildElement("Max_rel_step_increase");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the max relative step increase!" << "\n"; return ERROR; }
    Integrator_params->Max_rel_step_increase = std::stod(temp_param_var->GetText());

    // -------------------- Min relative step increase
    temp_param_var = Integrator_element->FirstChildElement("Min_rel_step_increase");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the min relative step increase!" << "\n"; return ERROR; }
    Integrator_params->Min_rel_step_increase = std::stod(temp_param_var->GetText());

    // -------------------- Max integration count 
    temp_param_var = Integrator_element->FirstChildElement("max_integration_count");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the maximum integration count!" << "\n"; return ERROR; }
    Integrator_params->Max_integration_count = std::stoi(temp_param_var->GetText());

    // -------------------- Max integration count 
    temp_param_var = Integrator_element->FirstChildElement("simpson_method_accuracy");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the simpson method accuracy parameter!" << "\n"; return ERROR; }
    Integrator_params->Simpson_accuracy = std::stod(temp_param_var->GetText());

    // -------------------- The step controller type
    temp_param_var = Integrator_element->FirstChildElement("Step_controller_type");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the step controller type!" << "\n"; return ERROR; }
    std::string Step_controller_type = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Step_controller_type.c_str()), "Gustafsson")) {

        Integrator_params->Simpson_accuracy = Gustafsson;

    }
    else if (0 == strcmp(static_cast<const char*>(Step_controller_type.c_str()), "PID")) {

        Integrator_params->Simpson_accuracy = PID;

    }
    else {

        std::cout << "Unsupported step controller type!" << "\n";

        return ERROR;

    }

    return OK;

}

Return_Values static parse_emission_model_params(tinyxml2::XMLElement* Emission_model_element, Initial_conditions_type* p_Init_conditions) {

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
    Observer_params->resolution_y = std::stoi(temp_param_var->GetText());

    // -------------------- Image resolution X
    temp_param_var = Observer_element->FirstChildElement("Resolution_x");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse X resolution!" << "\n"; return ERROR; }
    Observer_params->resolution_x = std::stoi(temp_param_var->GetText());

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
        NT_params->r_in = std::stod(temp_param_var->GetText());

        // -------------------- Outer disk radius
        temp_param_var = NT_element->FirstChildElement("r_out");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the outer NT disk radius!" << "\n"; return ERROR; }
        NT_params->r_out = std::stod(temp_param_var->GetText());

    }

    return OK;

}

Return_Values static parse_metric_parameters(tinyxml2::XMLElement* Metric_element, Metric_parameters_type* Metric_params) {

    tinyxml2::XMLElement* temp_param_var;

    std::string Metric_type = Metric_element->FirstChildElement("Metric_type")->GetText();

    if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Kerr")) {

        Metric_params->e_Spacetime = Kerr;

        temp_param_var = Metric_element->FirstChildElement("Spin_parameter");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the metric spin parameter!" << "\n"; return ERROR; }
        Metric_params->Spin = std::stod(temp_param_var->GetText());

    }
    else if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Wormhole")) {

        Metric_params->e_Spacetime = Wormhole;

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

        Metric_params->e_Spacetime = Janis_Newman_Winicour;

        temp_param_var = Metric_element->FirstChildElement("JNW_gamma");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Janis-Newman-Winicour metric parameter!" << "\n"; return ERROR; }
        Metric_params->JNW_Gamma_Parameter = std::stod(temp_param_var->GetText());

    }
    else if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Einstein-Gauss-Bonnet")) {

        Metric_params->e_Spacetime = Einstein_Gauss_Bonnet;

        temp_param_var = Metric_element->FirstChildElement("EGB_gamma");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Einstein-Gauss-Bonnet metric parameter!" << "\n"; return ERROR; }
        Metric_params->GB_Gamma_Parameter = std::stod(temp_param_var->GetText());

    }
    else if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Regular-Black-Hole")) {

        Metric_params->e_Spacetime = Reg_Black_Hole;

        temp_param_var = Metric_element->FirstChildElement("RBH_param");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the regular black hole metric parameter!" << "\n"; return ERROR; }
        Metric_params->RBH_Parameter = std::stod(temp_param_var->GetText());

    }
    else if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Black-Hole-w-Dark-Matter")) {

        Metric_params->e_Spacetime = BH_w_Dark_Matter;

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

Return_Values static parse_file_manager_params(tinyxml2::XMLElement* File_manager_element, File_manager_parameters* File_manager_params) {

    tinyxml2::XMLElement* temp_param_var;

    // -------------------- Sim mode 2 input file path
    temp_param_var = File_manager_element->FirstChildElement("Sim_mode_2_input_file_path");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the simulation mode 2 input file path!" << "\n"; return ERROR; }
    if (temp_param_var->GetText() != nullptr) { File_manager_params->Sim_mode_2_imput_path = temp_param_var->GetText(); }

    // -------------------- Output file path
    temp_param_var = File_manager_element->FirstChildElement("Output_file_directory");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the output file directory!" << "\n"; return ERROR; }
    File_manager_params->Output_file_directory = temp_param_var->GetText();

    // -------------------- Common file names
    temp_param_var = File_manager_element->FirstChildElement("Common_file_names");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the common output file name!" << "\n"; return ERROR; }
    if (temp_param_var->GetText() != nullptr) { File_manager_params->Common_file_names = temp_param_var->GetText(); }

    // -------------------- Vertex shader path
    temp_param_var = File_manager_element->FirstChildElement("Vert_shader_path");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the vertex shader path!" << "\n"; return ERROR; }
    File_manager_params->Vert_shader_path = temp_param_var->GetText();

    // -------------------- Fragment shader path
    temp_param_var = File_manager_element->FirstChildElement("Frag_shader_path");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the fragment shader path!" << "\n"; return ERROR; }
    File_manager_params->Frag_shader_path = temp_param_var->GetText();

    // -------------------- Truncate files flag
    temp_param_var = File_manager_element->FirstChildElement("Truncate_files");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the Truncate Files flag!" << "\n"; return ERROR; }
    File_manager_params->Truncate_files = std::stoi(temp_param_var->GetText());

    return OK;

}

Return_Values parse_simulation_input_XML(const std::string input_file_path, Initial_conditions_type* const p_Initial_conditions) {

    tinyxml2::XMLDocument xml_doc;

    tinyxml2::XMLError e_parse_result = xml_doc.LoadFile(static_cast<const char*>(input_file_path.c_str()));
    if (e_parse_result != tinyxml2::XML_SUCCESS) { return ERROR; }

    tinyxml2::XMLElement* Root_node = xml_doc.FirstChildElement("Simulation_Input");
    if (Root_node == nullptr) { return ERROR; }

    p_Initial_conditions->File_manager_params.Simulation_name = Root_node->Attribute("Simulation_Name");

    /* ====================================== Parse the average pitch angle flag and sample number ====================================== */

    tinyxml2::XMLElement* Average_pitch_angle_flag_element = Root_node->FirstChildElement("Average_emission_pitch_angle");
    if (Average_pitch_angle_flag_element == nullptr) { std::cout << "Failed to find the pitch angle averaging flag!" << "\n"; return ERROR; }
    p_Initial_conditions->Average_electron_pitch_angle = std::stoi(Average_pitch_angle_flag_element->GetText());

    tinyxml2::XMLElement* Average_pitch_angle_number_element = Root_node->FirstChildElement("Emission_pitch_angle_samples_to_average");
    if (Average_pitch_angle_number_element == nullptr) { std::cout << "Failed to find the number of pitch angle samples to average!" << "\n"; return ERROR; }
    p_Initial_conditions->Emission_pitch_angle_samples_to_average = std::stoi(Average_pitch_angle_number_element->GetText());

    /* ====================================== Parse the simulation mode specific settings ====================================== */

    tinyxml2::XMLElement* Simulation_mode_element = Root_node->FirstChildElement("Simulation_mode");
    if (Simulation_mode_element == nullptr) { std::cout << "Failed to find the simulation mode!" << "\n"; return ERROR; }
    p_Initial_conditions->Simulation_mode = std::stoi(Simulation_mode_element->GetText());

    tinyxml2::XMLElement* Sim_mode_2_param_number_element = Root_node->FirstChildElement("Sim_mode_2_param_value_number");
    if (Sim_mode_2_param_number_element == nullptr) { std::cout << "Failed to find the simulation mode 2 number of metric parameter values!" << "\n"; return ERROR; }
    p_Initial_conditions->Sim_mode_2_param_value_number = std::stoi(Sim_mode_2_param_number_element->GetText());

    tinyxml2::XMLElement* Sim_mode_3_X_init = Root_node->FirstChildElement("Sim_mode_3_X_init");
    if (Sim_mode_3_X_init == nullptr) { std::cout << "Failed to find sim mode 3 X init!" << "\n"; return ERROR; }
    p_Initial_conditions->Sim_mode_3_X_init = std::stod(Sim_mode_3_X_init->GetText());

    tinyxml2::XMLElement* Sim_mode_3_Y_init = Root_node->FirstChildElement("Sim_mode_3_Y_init");
    if (Sim_mode_3_Y_init == nullptr) { std::cout << "Failed to find sim mode 3 Y init!" << "\n"; return ERROR; }
    p_Initial_conditions->Sim_mode_3_Y_init = std::stod(Sim_mode_3_Y_init->GetText());

    /* ====================================== Parse the central object mass ====================================== */

    tinyxml2::XMLElement* Central_object_mass_element = Root_node->FirstChildElement("Central_object_mass");
    if (Central_object_mass_element == nullptr) { std::cout << "Failed to find the central object mass!" << "\n"; return ERROR; }
    p_Initial_conditions->central_object_mass = std::stod(Central_object_mass_element->GetText());

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

    /* ====================================== Parse the file paths ====================================== */

    tinyxml2::XMLElement* File_manager_element = Root_node->FirstChildElement("File_Manager");
    if (File_manager_element == nullptr) { std::cout << "Failed to find the File paths node!" << "\n"; return ERROR; }
    if (OK != parse_file_manager_params(File_manager_element, &p_Initial_conditions->File_manager_params)) { return ERROR; };

    return OK;

}