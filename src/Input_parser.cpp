#include "Input_parser.h"
#include "tinyxml2.h"

Return_Values static parse_hotspot_params(tinyxml2::XMLElement* Hotspot_element, Hotspot_model_parameters_type *Hotspot_params) {

    tinyxml2::XMLElement* temp_param_var;

    // -------------------- The threshold relative density
    temp_param_var = Hotspot_element->FirstChildElement("Threshold_relative_density");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot threshold relative density!" << "\n"; return ERROR; }
    Hotspot_params->Threshold_relative_density = std::stod(temp_param_var->GetText());

    // -------------------- The ensamble type

    temp_param_var = Hotspot_element->FirstChildElement("Ensamble_type");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot ensamble type!" << "\n"; return ERROR; }
    std::string Ensamble_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Kappa")) { Hotspot_params->Ensamble_type = e_Kappa_ensamble; }

    else if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Thermal")) { Hotspot_params->Ensamble_type = e_Thermal_ensamble; }

    else if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Phenomenological")) { Hotspot_params->Ensamble_type = e_Phenomenological_ensamble; }

    else { std::cout << "Unsupported ensamble type for the hotspot!" << "\n"; return ERROR; }
    
    // -------------------- The density profile

    temp_param_var = Hotspot_element->FirstChildElement("Density_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot density profile type!" << "\n"; return ERROR; }
    std::string Density_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Density_type_string.c_str()), "Gaussian")) { Hotspot_params->Density_profile_type = e_Gaussian; }

    else if (0 == strcmp(static_cast<const char*>(Density_type_string.c_str()), "Sphere")) { Hotspot_params->Density_profile_type = e_Spherical; }

    else if (0 == strcmp(static_cast<const char*>(Density_type_string.c_str()), "Hybrid_power_law_gaussian")) { Hotspot_params->Density_profile_type = e_Hybrid_power_gaussian; }

    else { std::cout << "Unsupported density profile type for the hotspot!" << "\n"; return ERROR; }

    // -------------------- The temperature profile

    temp_param_var = Hotspot_element->FirstChildElement("Temperature_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot temperature profile type!" << "\n"; return ERROR; }
    std::string Temperature_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Temperature_type_string.c_str()), "Gaussian")) { Hotspot_params->Temperature_profile_type = e_Gaussian; }

    else if (0 == strcmp(static_cast<const char*>(Temperature_type_string.c_str()), "Sphere")) { Hotspot_params->Temperature_profile_type = e_Spherical; }

    else if (0 == strcmp(static_cast<const char*>(Temperature_type_string.c_str()), "Hybrid_power_law_gaussian")) { Hotspot_params->Temperature_profile_type = e_Hybrid_power_gaussian; }

    else { std::cout << "Unsupported temperature profile type for the hotspot!" << "\n"; return ERROR; }

    // -------------------- The velocity profile

    temp_param_var = Hotspot_element->FirstChildElement("Velocity_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot velocity profile node!" << "\n"; return ERROR; }

    temp_param_var = Hotspot_element->FirstChildElement("Velocity_profile")->FirstChildElement("Type");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot velocity profile type!" << "\n"; return ERROR; }
    std::string Velocity_profile_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Velocity_profile_string.c_str()), "Theta Dependant")) { Hotspot_params->Velocity_profile_type = e_Theta_dependant; }

    else if (0 == strcmp(static_cast<const char*>(Velocity_profile_string.c_str()), "Keplarian")) { Hotspot_params->Velocity_profile_type = e_Keplarian; }

    else if (0 == strcmp(static_cast<const char*>(Velocity_profile_string.c_str()), "Circular Fixed Rate")) { Hotspot_params->Velocity_profile_type = e_Circular_fixed_rate; }

    else{std::cout << "Unsupported velocity profile type for the hotspot!" << "\n"; return ERROR; }

    temp_param_var = Hotspot_element->FirstChildElement("Velocity_profile")->FirstChildElement("Radial_velocity_fraction");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot radial velocity fraction!" << "\n"; return ERROR; }
    Hotspot_params->Radial_velocity_fraction = std::stod(temp_param_var->GetText());

    // -------------------- The density sclae factor
    temp_param_var = Hotspot_element->FirstChildElement("Density_scale_factor");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot density scale factor!" << "\n"; return ERROR; }
    Hotspot_params->Electron_density_scale = std::stod(temp_param_var->GetText());

    // -------------------- The density sclae factor
    temp_param_var = Hotspot_element->FirstChildElement("Temperature_scale_factor");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot temperature scale factor!" << "\n"; return ERROR; }
    Hotspot_params->Electron_temperature_scale = std::stod(temp_param_var->GetText());

    // -------------------- The magnetic field geometry r component
    temp_param_var = Hotspot_element->FirstChildElement("Mag_field_geometry_r");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot magnetic field geometry r component!" << "\n"; return ERROR; }
    Hotspot_params->Mag_field_geometry[e_r - 1] = std::stod(temp_param_var->GetText());

    // -------------------- The magnetic field geometry theta component
    temp_param_var = Hotspot_element->FirstChildElement("Mag_field_geometry_theta");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot magnetic field geometry theta component!" << "\n"; return ERROR; }
    Hotspot_params->Mag_field_geometry[e_theta - 1] = std::stod(temp_param_var->GetText());

    // -------------------- The magnetic field geometry phi component
    temp_param_var = Hotspot_element->FirstChildElement("Mag_field_geometry_phi");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot magnetic field geometry phi component!" << "\n"; return ERROR; }
    Hotspot_params->Mag_field_geometry[e_phi - 1] = std::stod(temp_param_var->GetText());

    Hotspot_params->Position[e_t] = 0.0;

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
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot magnetization!" << "\n"; return ERROR; }
    Hotspot_params->Magnetization = std::stod(temp_param_var->GetText());

    /* ================================== The gaussian profile parameters ================================== */

    tinyxml2::XMLElement* Gaussian_node = Hotspot_element->FirstChildElement("Gaussian_profile");
    tinyxml2::XMLElement* Spherical_node = Hotspot_element->FirstChildElement("Spherical_profile");
    tinyxml2::XMLElement* Hybrid_node = Hotspot_element->FirstChildElement("Hybrid_power_law_gaussian_profile");

    if (e_Gaussian == Hotspot_params->Density_profile_type) {

        if (Gaussian_node == nullptr) { std::cout << "Failed to parse the hotspot gaussian profile node!" << "\n"; return ERROR; }

        // -------------------- The Gaussian density standard deviation
        temp_param_var = Gaussian_node->FirstChildElement("Density_spread");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot density Gaussian spread!" << "\n"; return ERROR; }
        Hotspot_params->Profile_params.Density_gaussian_spread = std::stod(temp_param_var->GetText());

    }
    else if (e_Hybrid_power_gaussian == Hotspot_params->Density_profile_type) {

        if (Hybrid_node == nullptr) { std::cout << "Failed to parse the hotspot hybrid profile node!" << "\n"; return ERROR; }

        // -------------------- The Gaussian density standard deviation
        temp_param_var = Hybrid_node->FirstChildElement("Density_spread");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot density Gaussian spread!" << "\n"; return ERROR; }
        Hotspot_params->Profile_params.Density_gaussian_spread = std::stod(temp_param_var->GetText());

        // -------------------- The radial power law exponent
        temp_param_var = Hybrid_node->FirstChildElement("Density_power_law_power");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot density power law exponent!" << "\n"; return ERROR; }
        Hotspot_params->Profile_params.Density_power_law_power = std::stod(temp_param_var->GetText());

        // -------------------- The radial power law scale
        temp_param_var = Hybrid_node->FirstChildElement("Density_power_law_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot density power law scale spread!" << "\n"; return ERROR; }
        Hotspot_params->Profile_params.Density_power_law_scale = std::stod(temp_param_var->GetText());

    }

    if (e_Gaussian == Hotspot_params->Temperature_profile_type) {

        if (Gaussian_node == nullptr) { std::cout << "Failed to parse the temperature Gaussian profile node!" << "\n"; return ERROR; }

        // -------------------- The Gaussian temperature standard devivation
        temp_param_var = Gaussian_node->FirstChildElement("Temperature_spread");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the temperature Gaussian spread!" << "\n"; return ERROR; }
        Hotspot_params->Profile_params.Temperature_gaussian_spread = std::stod(temp_param_var->GetText());


    }
    else if (e_Hybrid_power_gaussian == Hotspot_params->Temperature_profile_type) {

        if (Hybrid_node == nullptr) { std::cout << "Failed to parse the hotspot hybrid profile node!" << "\n"; return ERROR; }

        // -------------------- The Gaussian density standard deviation
        temp_param_var = Hybrid_node->FirstChildElement("Temperature_spread");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot temperature Gaussian spread!" << "\n"; return ERROR; }
        Hotspot_params->Profile_params.Temperature_gaussian_spread = std::stod(temp_param_var->GetText());

        // -------------------- The radial power law exponent
        temp_param_var = Hybrid_node->FirstChildElement("Temperature_power_law_power");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot temperature power law exponent!" << "\n"; return ERROR; }
        Hotspot_params->Profile_params.Temperature_power_law_power = std::stod(temp_param_var->GetText());

        // -------------------- The radial power law scale
        temp_param_var = Hybrid_node->FirstChildElement("Temperature_power_law_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot temperature power law scale spread!" << "\n"; return ERROR; }
        Hotspot_params->Profile_params.Temperature_power_law_scale = std::stod(temp_param_var->GetText());

    }

    if (e_Spherical == Hotspot_params->Temperature_profile_type || e_Spherical == Hotspot_params->Density_profile_type) {

        if (Spherical_node == nullptr) { std::cout << "Failed to parse the temperature spherical profile node!" << "\n"; return ERROR; }

        // -------------------- The Spherical radius
        temp_param_var = Spherical_node->FirstChildElement("Radius");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Spherical hotspot radius!" << "\n"; return ERROR; }
        Hotspot_params->Profile_params.Radius = std::stod(temp_param_var->GetText());

    }

    // -------------------- The Temporal spread
    temp_param_var = Hotspot_element->FirstChildElement("Temporal_spread");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the temporal spread!" << "\n"; return ERROR; }
    Hotspot_params->Profile_params.Temporal_gaussian_spread = std::stod(temp_param_var->GetText());

    // -------------------- The Coordiante time at max emission
    temp_param_var = Hotspot_element->FirstChildElement("Coord_time_at_max");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the coordiante time at max emission!" << "\n"; return ERROR; }
    Hotspot_params->Profile_params.Coord_time_offset = std::stod(temp_param_var->GetText());

    // -------------------- The hotspot Magnetic field magnitude profile enum
    temp_param_var = Hotspot_element->FirstChildElement("Mag_field_magnitude_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot magnetic field geometry enum!" << "\n"; return ERROR; }
    std::string Mag_field_magnitude_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Mag_field_magnitude_string.c_str()), "Power_law_based")) { Hotspot_params->e_Mag_field_magnitude_profile = Power_law_based; }

    else if (0 == strcmp(static_cast<const char*>(Mag_field_magnitude_string.c_str()), "Magnetization_based")) { Hotspot_params->e_Mag_field_magnitude_profile = Magnetization_based; }

    else { std::cout << "Unsupported velocity profile type for the disk!" << "\n"; return ERROR; }

    // -------------------- The disk Magnetic field radial scale
    temp_param_var = Hotspot_element->FirstChildElement("Mag_field_radial_scale");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot magnetic field radial scale!" << "\n"; return ERROR; }
    Hotspot_params->Mag_field_radial_scale = std::stod(temp_param_var->GetText());

    // -------------------- The disk Magnetic field geometry enum
    temp_param_var = Hotspot_element->FirstChildElement("Mag_field_geometry");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot magnetic field geometry enum!" << "\n"; return ERROR; }
    std::string Mag_field_geometry_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Mag_field_geometry_string.c_str()), "Toroidal")) { Hotspot_params->e_Mag_field_geometry = Toroidal; }

    else if (0 == strcmp(static_cast<const char*>(Mag_field_geometry_string.c_str()), "Vertical")) { Hotspot_params->e_Mag_field_geometry = Vertical; }

    else if (0 == strcmp(static_cast<const char*>(Mag_field_geometry_string.c_str()), "Constant")) { Hotspot_params->e_Mag_field_geometry = Constant; }

    else { std::cout << "Unsupported velocity profile type for the disk!" << "\n"; return ERROR; }

    // -------------------- The hotspot Magnetic field magnitude scale
    temp_param_var = Hotspot_element->FirstChildElement("Mag_field_magnitude_scale");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot magnetic field scale!" << "\n"; return ERROR; }
    Hotspot_params->Mag_field_magnitude_scale = std::stod(temp_param_var->GetText());

    // -------------------- The hotspot Magnetic field power law power
    temp_param_var = Hotspot_element->FirstChildElement("Mag_field_power");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the hotspot magnetic field power law power!" << "\n"; return ERROR; }
    Hotspot_params->Mag_field_power = std::stod(temp_param_var->GetText());

    return OK;

}

Return_Values static parse_disk_params(tinyxml2::XMLElement* Accretion_disk_element, Disk_model_parameters_type* Disk_params) {

    tinyxml2::XMLElement* temp_param_var;

    // -------------------- The disk model

    temp_param_var = Accretion_disk_element->FirstChildElement("Disk_Model");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk model!" << "\n"; return ERROR; }
    std::string Profile_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Profile_type_string.c_str()), "Phenom_RIAF_1")) { Disk_params->e_Disk_model = e_Phenom_RIAF_1; }

    else if (0 == strcmp(static_cast<const char*>(Profile_type_string.c_str()), "Phenom_RIAF_2")) { Disk_params->e_Disk_model = e_Phenom_RIAF_2; }

    else if (0 == strcmp(static_cast<const char*>(Profile_type_string.c_str()), "Colab_test_1")) { Disk_params->e_Disk_model = e_Colab_test_1; }

    else if (0 == strcmp(static_cast<const char*>(Profile_type_string.c_str()), "Page-Thorne")) { Disk_params->e_Disk_model = e_Page_Thorne; }

    else if (0 == strcmp(static_cast<const char*>(Profile_type_string.c_str()), "Debug_constant_density")) { Disk_params->e_Disk_model = e_Debug_constant_density; }

    else { std::cout << "Unsupported disk model!" << "\n"; return ERROR; }

    // ----------------- Parse the Page-Thorne parameters

    if (e_Page_Thorne == Disk_params->e_Disk_model) {

        temp_param_var = Accretion_disk_element->FirstChildElement("r_in");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the inner Page-Thorne radius!" << "\n"; return ERROR; }
        Disk_params->Page_Thorne_params.r_in = std::stod(temp_param_var->GetText());

        temp_param_var = Accretion_disk_element->FirstChildElement("r_out");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the inner Page-Thorne radius!" << "\n"; return ERROR; }
        Disk_params->Page_Thorne_params.r_out = std::stod(temp_param_var->GetText());

        temp_param_var = Accretion_disk_element->FirstChildElement("Mag_field_geometry");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field geometry enum!" << "\n"; return ERROR; }
        std::string Mag_field_geometry_string = temp_param_var->GetText();

        if (0 == strcmp(static_cast<const char*>(Mag_field_geometry_string.c_str()), "Toroidal")) { Disk_params->e_Mag_field_geometry = Toroidal; }

        else if (0 == strcmp(static_cast<const char*>(Mag_field_geometry_string.c_str()), "Vertical")) { Disk_params->e_Mag_field_geometry = Vertical; }

        else if (0 == strcmp(static_cast<const char*>(Mag_field_geometry_string.c_str()), "Constant")) { Disk_params->e_Mag_field_geometry = Constant; }

        else { std::cout << "Unsupported magnetic field profile type for the disk!" << "\n"; return ERROR; }

        // -------------------- The disk magnetic field geometry r component
        temp_param_var = Accretion_disk_element->FirstChildElement("Mag_field_geometry_r");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field geometry r component!" << "\n"; return ERROR; }
        Disk_params->Mag_field_geometry[e_r - 1] = std::stod(temp_param_var->GetText());

        // -------------------- The disk magnetic field geometry theta component
        temp_param_var = Accretion_disk_element->FirstChildElement("Mag_field_geometry_theta");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field geometry theta component!" << "\n"; return ERROR; }
        Disk_params->Mag_field_geometry[e_theta - 1] = std::stod(temp_param_var->GetText());

        // -------------------- The disk magnetic field geometry phi component 
        temp_param_var = Accretion_disk_element->FirstChildElement("Mag_field_geometry_phi");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field geometry phi component!" << "\n"; return ERROR; }
        Disk_params->Mag_field_geometry[e_phi - 1] = std::stod(temp_param_var->GetText());

        return OK;

    }

    /* ======================================== Common paramaters ======================================== */

    tinyxml2::XMLElement* Common_paramaters_element = Accretion_disk_element->FirstChildElement("Common_parameters");
    if (Common_paramaters_element == nullptr) { std::cout << "Failed to find the Common parameters element!" << "\n"; return ERROR; }

    // -------------------- The threshold relative density
    temp_param_var = Common_paramaters_element->FirstChildElement("Threshold_relative_density");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk threshold relative density!" << "\n"; return ERROR; }
    Disk_params->Threshold_relative_density = std::stod(temp_param_var->GetText());

    // -------------------- The ensamble type

    temp_param_var = Common_paramaters_element->FirstChildElement("Ensamble_type");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk ensamble type!" << "\n"; return ERROR; }
    std::string Ensamble_type_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Kappa")) { Disk_params->Ensamble_type = e_Kappa_ensamble; }

    else if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Thermal")) { Disk_params->Ensamble_type = e_Thermal_ensamble; }

    else if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Phenomenological")) { Disk_params->Ensamble_type = e_Phenomenological_ensamble; }

    else if (0 == strcmp(static_cast<const char*>(Ensamble_type_string.c_str()), "Debug_constant_functions")) { Disk_params->Ensamble_type = e_Debug_constant_functions; }

    else { std::cout << "Unsupported ensamble type for the disk!" << "\n"; return ERROR; }

    // -------------------- The velocity profile

    temp_param_var = Common_paramaters_element->FirstChildElement("Velocity_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk velocity profile node!" << "\n"; return ERROR; }

    temp_param_var = Common_paramaters_element->FirstChildElement("Velocity_profile")->FirstChildElement("Type");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk velocity profile type!" << "\n"; return ERROR; }
    std::string Velocity_profile_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Velocity_profile_string.c_str()), "Theta Dependant")) { Disk_params->Velocity_profile_type = e_Theta_dependant; }

    else if (0 == strcmp(static_cast<const char*>(Velocity_profile_string.c_str()), "Keplarian")) { Disk_params->Velocity_profile_type = e_Keplarian; }

    else { std::cout << "Unsupported velocity profile type for the disk!" << "\n"; return ERROR; }

    temp_param_var = Common_paramaters_element->FirstChildElement("Velocity_profile")->FirstChildElement("Radial_velocity_fraction");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the radial velocity fraction!" << "\n"; return ERROR; }
    Disk_params->Radial_velocity_fraction = std::stod(temp_param_var->GetText());

    // -------------------- The disk electron density scale factor
    temp_param_var = Common_paramaters_element->FirstChildElement("Density_scale_factor");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk density scale factor!" << "\n"; return ERROR; }
    Disk_params->Electron_density_scale = std::stod(temp_param_var->GetText());

    // -------------------- The disk electron temperature scale factor
    temp_param_var = Common_paramaters_element->FirstChildElement("Temperature_scale_factor");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk temperature scale factor!" << "\n"; return ERROR; }
    Disk_params->Electron_temperature_scale = std::stod(temp_param_var->GetText());

    // -------------------- The disk magnetic field geometry r component
    temp_param_var = Common_paramaters_element->FirstChildElement("Mag_field_geometry_r");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field geometry r component!" << "\n"; return ERROR; }
    Disk_params->Mag_field_geometry[e_r - 1] = std::stod(temp_param_var->GetText());

    // -------------------- The disk magnetic field geometry theta component
    temp_param_var = Common_paramaters_element->FirstChildElement("Mag_field_geometry_theta");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field geometry theta component!" << "\n"; return ERROR; }
    Disk_params->Mag_field_geometry[e_theta - 1] = std::stod(temp_param_var->GetText());

    // -------------------- The disk magnetic field geometry phi component 
    temp_param_var = Common_paramaters_element->FirstChildElement("Mag_field_geometry_phi");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field geometry phi component!" << "\n"; return ERROR; }
    Disk_params->Mag_field_geometry[e_phi - 1] = std::stod(temp_param_var->GetText());

    // -------------------- The disk magnetization
    temp_param_var = Common_paramaters_element->FirstChildElement("Magnetization"); 
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk Magnetization!" << "\n"; return ERROR; }
    Disk_params->Magnetization = std::stod(temp_param_var->GetText());

    // -------------------- The disk Magnetic field magnitude scale
    temp_param_var = Common_paramaters_element->FirstChildElement("Mag_field_magnitude_scale");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field scale!" << "\n"; return ERROR; }
    Disk_params->Mag_field_magnitude_scale = std::stod(temp_param_var->GetText());

    // -------------------- The disk Magnetic field power law power
    temp_param_var = Common_paramaters_element->FirstChildElement("Mag_field_power");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field power law power!" << "\n"; return ERROR; }
    Disk_params->Mag_field_power = std::stod(temp_param_var->GetText());

    // -------------------- The disk Magnetic field radial scale
    temp_param_var = Common_paramaters_element->FirstChildElement("Mag_field_radial_scale");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field radial scale!" << "\n"; return ERROR; }
    Disk_params->Mag_field_radial_scale = std::stod(temp_param_var->GetText());

    // -------------------- The disk Magnetic field geometry enum
    temp_param_var = Common_paramaters_element->FirstChildElement("Mag_field_geometry");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field geometry enum!" << "\n"; return ERROR; }
    std::string Mag_field_geometry_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Mag_field_geometry_string.c_str()), "Toroidal")) { Disk_params->e_Mag_field_geometry = Toroidal;}

    else if (0 == strcmp(static_cast<const char*>(Mag_field_geometry_string.c_str()), "Vertical")) { Disk_params->e_Mag_field_geometry = Vertical;}

    else if (0 == strcmp(static_cast<const char*>(Mag_field_geometry_string.c_str()), "Constant")) {Disk_params->e_Mag_field_geometry = Constant;}

    else { std::cout << "Unsupported velocity profile type for the disk!" << "\n"; return ERROR; }

    // -------------------- The disk Magnetic field magnitude profile enum
    temp_param_var = Common_paramaters_element->FirstChildElement("Mag_field_magnitude_profile");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the disk magnetic field geometry enum!" << "\n"; return ERROR; }
    std::string Mag_field_magnitude_string = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Mag_field_magnitude_string.c_str()), "Power_law_based")) { Disk_params->e_Mag_field_magnitude_profile = Power_law_based; }

    else if (0 == strcmp(static_cast<const char*>(Mag_field_magnitude_string.c_str()), "Magnetization_based")) { Disk_params->e_Mag_field_magnitude_profile = Magnetization_based; }

    else { std::cout << "Unsupported velocity profile type for the disk!" << "\n"; return ERROR; }

    /* ======================================== Disk model parameters ======================================== */

    tinyxml2::XMLElement* Common_RIAF_element = Accretion_disk_element->FirstChildElement("Common_RIAF_profile");
    tinyxml2::XMLElement* Colab_test_1_element = Accretion_disk_element->FirstChildElement("Colab_test_1_profile");

    if (e_Colab_test_1 != Disk_params->e_Disk_model) {

        if (Common_RIAF_element == nullptr) { std::cout << "Failed to find the common RIAF profile element!" << "\n"; return ERROR; }

        // -------------------- The density radial power law
        temp_param_var = Common_RIAF_element->FirstChildElement("Density_power_law_power");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the common RIAF disk density radial power law!" << "\n"; return ERROR; }
        Disk_params->Common_RIAF_params.Density_power_law_power = std::stod(temp_param_var->GetText());

        // -------------------- The temperature radial power law
        temp_param_var = Common_RIAF_element->FirstChildElement("Temperature_power_law_power");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the common RIAF disk temperature radial power law!" << "\n"; return ERROR; }
        Disk_params->Common_RIAF_params.Temperature_power_law_power = std::stod(temp_param_var->GetText());

        // -------------------- The density radial cutoff scale
        temp_param_var = Common_RIAF_element->FirstChildElement("Density_cutoff_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the common RIAF disk density radial cutoff scale!" << "\n"; return ERROR; }
        Disk_params->Common_RIAF_params.Density_cutoff_scale = std::stod(temp_param_var->GetText());

        // -------------------- The temperature radial cutoff scale
        temp_param_var = Common_RIAF_element->FirstChildElement("Temperature_cutoff_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the common RIAF disk temperature radial cutoff scale!" << "\n"; return ERROR; }
        Disk_params->Common_RIAF_params.Temperature_cutoff_scale = std::stod(temp_param_var->GetText());

        // -------------------- The density r_cutoff
        temp_param_var = Common_RIAF_element->FirstChildElement("Density_cutoff_radius");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the common RIAF disk density r_cutoff!" << "\n"; return ERROR; }
        Disk_params->Common_RIAF_params.Density_cutoff_radius = std::stod(temp_param_var->GetText());

        // -------------------- The temperature r_cutoff
        temp_param_var = Common_RIAF_element->FirstChildElement("Temperature_cutoff_radius");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the common RIAF disk temperature r_cutoff!" << "\n"; return ERROR; }
        Disk_params->Common_RIAF_params.Temperature_cutoff_radius = std::stod(temp_param_var->GetText());

        // -------------------- The density r_0
        temp_param_var = Common_RIAF_element->FirstChildElement("Density_power_law_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the common RIAF disk density power law scale!" << "\n"; return ERROR; }
        Disk_params->Common_RIAF_params.Density_power_law_scale = std::stod(temp_param_var->GetText());

        // -------------------- The temperature r_0
        temp_param_var = Common_RIAF_element->FirstChildElement("Density_power_law_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the common RIAF disk temperature power law scale!" << "\n"; return ERROR; }
        Disk_params->Common_RIAF_params.Temperature_power_law_scale = std::stod(temp_param_var->GetText());

        // -------------------- The disk opening angle
        temp_param_var = Common_RIAF_element->FirstChildElement("Opening_angle");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the common RIAF disk opening angle parameter!" << "\n"; return ERROR; }
        Disk_params->Common_RIAF_params.Disk_opening_angle = std::stod(temp_param_var->GetText());

    }
    else {

        if (Colab_test_1_element == nullptr) { std::cout << "Failed to find the Colab Test 1 profile element !" << "\n"; return ERROR; }

        // -------------------- The density radial power law
        temp_param_var = Colab_test_1_element->FirstChildElement("Radial_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Colab Test 1 disk radial scale!" << "\n"; return ERROR; }
        Disk_params->Colab_test_1_params.Radial_scale = std::stod(temp_param_var->GetText());

        // -------------------- The density radial power law
        temp_param_var = Colab_test_1_element->FirstChildElement("Vertical_scale");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Colab Test 1 disk vertical scale!" << "\n"; return ERROR; }
        Disk_params->Colab_test_1_params.Vertical_scale = std::stod(temp_param_var->GetText());

    }

    return OK;

}

Return_Values static parse_integrator_params(tinyxml2::XMLElement* Integrator_element, Integrator_parameters_type* Integrator_params) {

    tinyxml2::XMLElement* temp_param_var;

    // -------------------- Init stepsize
    temp_param_var = Integrator_element->FirstChildElement("init_stepsize");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the initial step size!" << "\n"; return ERROR; }
    Integrator_params->Init_stepzie = std::stod(temp_param_var->GetText());

    // -------------------- RK78 abs accuracy
    temp_param_var = Integrator_element->FirstChildElement("RK78_abs_accuracy");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the RK78 abs accuracy parameter!" << "\n"; return ERROR; }
    Integrator_params->RK_78_abs_accuracy = std::stod(temp_param_var->GetText());

    // -------------------- RK78 rel accuracy
    temp_param_var = Integrator_element->FirstChildElement("RK78_rel_accuracy");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the RK78 rel accuracy parameter!" << "\n"; return ERROR; }
    Integrator_params->RK_78_rel_accuracy = std::stod(temp_param_var->GetText());

    // -------------------- ESDIRK54 abs accuracy
    temp_param_var = Integrator_element->FirstChildElement("ESDIRK54_abs_accuracy");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the ESDIRK54 abs accuracy parameter!" << "\n"; return ERROR; }
    Integrator_params->ESDIRK54_abs_accuracy = std::stod(temp_param_var->GetText());

    // -------------------- ESDIRK54 rel accuracy
    temp_param_var = Integrator_element->FirstChildElement("ESDIRK54_rel_accuracy");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the ESDIRK54 rel accuracy parameter!" << "\n"; return ERROR; }
    Integrator_params->ESDIRK54_rel_accuracy = std::stod(temp_param_var->GetText());

    // -------------------- Step controller safety 1
    temp_param_var = Integrator_element->FirstChildElement("step_controller_safety_factor_1");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the step controller safety parameter 1!" << "\n"; return ERROR; }
    Integrator_params->Safety_1 = std::stod(temp_param_var->GetText());

    // -------------------- Step controller safety 2
    temp_param_var = Integrator_element->FirstChildElement("step_controller_safety_factor_2");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the step controller safety parameter 2!" << "\n"; return ERROR; }
    Integrator_params->Safety_2 = std::stod(temp_param_var->GetText());

    // -------------------- RK78 PID controller I gain
    temp_param_var = Integrator_element->FirstChildElement("RK78_PID_controller_I_gain");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the RK78 PID controller I gain!" << "\n"; return ERROR; }
    Integrator_params->RK78_PID_gain_I = std::stod(temp_param_var->GetText());

    // -------------------- RK78 PID controller P gain
    temp_param_var = Integrator_element->FirstChildElement("RK78_PID_controller_P_gain");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the RK78 PID controller P gain!" << "\n"; return ERROR; }
    Integrator_params->RK78_PID_gain_P = std::stod(temp_param_var->GetText());

    // -------------------- RK78 PID controller D gain
    temp_param_var = Integrator_element->FirstChildElement("RK78_PID_controller_D_gain");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the RK78 PID controller D gain!" << "\n"; return ERROR; }
    Integrator_params->RK78_PID_gain_D = std::stod(temp_param_var->GetText());

    // -------------------- RK78 Gustafsson controller k_1 gain
    temp_param_var = Integrator_element->FirstChildElement("RK78_Gustafsson_controller_k_1");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the RK78 Gustafsson controller k_1 gain!" << "\n"; return ERROR; }
    Integrator_params->RK78_Gustafsson_k1 = std::stod(temp_param_var->GetText());

    // -------------------- RK78 Gustafsson controller k_2 gain
    temp_param_var = Integrator_element->FirstChildElement("RK78_Gustafsson_controller_k_2");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the RK78 Gustafsson controller k_2 gain!" << "\n"; return ERROR; }
    Integrator_params->RK78_Gustafsson_k2 = std::stod(temp_param_var->GetText());

    // -------------------- ESDIRK54 PID controller I gain
    temp_param_var = Integrator_element->FirstChildElement("ESDIRK54_PID_controller_I_gain");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the ESDIRK54 PID controller I gain!" << "\n"; return ERROR; }
    Integrator_params->ESDIRK54_PID_gain_I = std::stod(temp_param_var->GetText());

    // -------------------- ESDIRK54 PID controller P gain
    temp_param_var = Integrator_element->FirstChildElement("ESDIRK54_PID_controller_P_gain");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the ESDIRK54 PID controller P gain!" << "\n"; return ERROR; }
    Integrator_params->ESDIRK54_PID_gain_P = std::stod(temp_param_var->GetText());

    // -------------------- ESDIRK54 PID controller D gain
    temp_param_var = Integrator_element->FirstChildElement("ESDIRK54_PID_controller_D_gain");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the ESDIRK54 PID controller D gain!" << "\n"; return ERROR; }
    Integrator_params->ESDIRK54_PID_gain_D = std::stod(temp_param_var->GetText());

    // -------------------- ESDIRK54 Gustafsson controller k_1 gain
    temp_param_var = Integrator_element->FirstChildElement("ESDIRK54_Gustafsson_controller_k_1");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the ESDIRK54 Gustafsson controller k_1 gain!" << "\n"; return ERROR; }
    Integrator_params->ESDIRK54_Gustafsson_k1 = std::stod(temp_param_var->GetText());

    // -------------------- ESDIRK54 Gustafsson controller k_2 gain
    temp_param_var = Integrator_element->FirstChildElement("ESDIRK54_Gustafsson_controller_k_2");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the ESDIRK54 Gustafsson controller k_2 gain!" << "\n"; return ERROR; }
    Integrator_params->ESDIRK54_Gustafsson_k2 = std::stod(temp_param_var->GetText());


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

    // -------------------- The adaptive Simpson integral solver accuracy parameter
    temp_param_var = Integrator_element->FirstChildElement("simpson_method_accuracy");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the simpson method accuracy parameter!" << "\n"; return ERROR; }
    Integrator_params->Simpson_accuracy = std::stod(temp_param_var->GetText());

    // -------------------- Max affine parameter value
    temp_param_var = Integrator_element->FirstChildElement("max_affine_parameter");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the max affine parameter value!" << "\n"; return ERROR; }
    Integrator_params->Max_affine_param = std::stod(temp_param_var->GetText());

    // -------------------- Use adaptive step flag
    temp_param_var = Integrator_element->FirstChildElement("use_adaptive_step");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the use adaptive step flag!" << "\n"; return ERROR; }
    Integrator_params->Use_adaptive_step = std::stoi(temp_param_var->GetText());

    // -------------------- The step controller type
    temp_param_var = Integrator_element->FirstChildElement("Step_controller_type");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the step controller type!" << "\n"; return ERROR; }
    std::string Step_controller_type = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Step_controller_type.c_str()), "Gustafsson")) {

        Integrator_params->Controller_type = Gustafsson;

    }
    else if (0 == strcmp(static_cast<const char*>(Step_controller_type.c_str()), "PID")) {

        Integrator_params->Controller_type = PID;

    }
    else {

        std::cout << "Unsupported step controller type!" << "\n";

        return ERROR;

    }

    // -------------------- Max stepsize
    temp_param_var = Integrator_element->FirstChildElement("max_stepsize");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the max step size!" << "\n"; return ERROR; }
    Integrator_params->Max_stepsize = std::stod(temp_param_var->GetText());

    // -------------------- Radiative transfer integrator
    temp_param_var = Integrator_element->FirstChildElement("radiative_transfer_integrator_type");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the radiative transfer integrator type!" << "\n"; return ERROR; }
    std::string Radiative_transfer_integrator_type = temp_param_var->GetText();

    if (0 == strcmp(static_cast<const char*>(Radiative_transfer_integrator_type.c_str()), "Analytic")) {

        Integrator_params->e_Radiative_transfer_integrator = Analytic;

    }
    else if (0 == strcmp(static_cast<const char*>(Radiative_transfer_integrator_type.c_str()), "Implicit Trapezoid")) {

        Integrator_params->e_Radiative_transfer_integrator = Implicit_Trapezoid;

    }
    else if (0 == strcmp(static_cast<const char*>(Radiative_transfer_integrator_type.c_str()), "RK5")) {

        Integrator_params->e_Radiative_transfer_integrator = RK5;
    
    }
    else {

        std::cout << "Unsupported radiative transfer integrator type!" << "\n";

        return ERROR;

    }

    // -------------------- Default gegodesic integrator
    temp_param_var = Integrator_element->FirstChildElement("default_geodesic_integrator_type");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the geodesic integrator type!" << "\n"; return ERROR; }
    std::string Default_geodesic_integrator_type = temp_param_var->GetText();


    if (0 == strcmp(static_cast<const char*>(Default_geodesic_integrator_type.c_str()), "RK78_Fehlberg")) {

        Integrator_params->e_Default_geodesic_integrator = RK78_Fehlberg;

    }
    else if (0 == strcmp(static_cast<const char*>(Default_geodesic_integrator_type.c_str()), "RK78_DP")) {

        Integrator_params->e_Default_geodesic_integrator = RK78_DP;

    }
    else if (0 == strcmp(static_cast<const char*>(Default_geodesic_integrator_type.c_str()), "ESDIRK54")) {

        Integrator_params->e_Default_geodesic_integrator = ESDIRK54;

    }
    else {

        std::cout << "Unsupported geodesic integrator type!" << "\n";

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

    if ((e_Debug_constant_functions == p_Init_conditions->Disk_params.Ensamble_type && 0 != p_Init_conditions->Disk_params.Electron_density_scale) ||
         e_Debug_constant_functions == p_Init_conditions->Hotspot_params.Ensamble_type && 0 != p_Init_conditions->Hotspot_params.Electron_density_scale) {

        // -------------------- Debug_j_I_value
        temp_param_var = Emission_model_element->FirstChildElement("Debug_j_I_value");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Debug_j_I_value value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Debug_j_I_value = std::stod(temp_param_var->GetText());

        // -------------------- Debug_j_Q_value
        temp_param_var = Emission_model_element->FirstChildElement("Debug_j_Q_value");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Debug_j_Q_value value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Debug_j_Q_value = std::stod(temp_param_var->GetText());

        // -------------------- Debug_j_U_value
        temp_param_var = Emission_model_element->FirstChildElement("Debug_j_U_value");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Debug_j_U_value value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Debug_j_U_value = std::stod(temp_param_var->GetText());

        // -------------------- Debug_j_V_value
        temp_param_var = Emission_model_element->FirstChildElement("Debug_j_V_value");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Debug_j_V_value value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Debug_j_V_value = std::stod(temp_param_var->GetText());

        // -------------------- Debug_alpha_I_value
        temp_param_var = Emission_model_element->FirstChildElement("Debug_alpha_I_value");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Debug_alpha_I_value value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Debug_alpha_I_value = std::stod(temp_param_var->GetText());

        // -------------------- Debug_alpha_Q_value
        temp_param_var = Emission_model_element->FirstChildElement("Debug_alpha_Q_value");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Debug_alpha_Q_value value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Debug_alpha_Q_value = std::stod(temp_param_var->GetText());

        // -------------------- Debug_alpha_U_value
        temp_param_var = Emission_model_element->FirstChildElement("Debug_alpha_U_value");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Debug_alpha_U_value value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Debug_alpha_U_value = std::stod(temp_param_var->GetText());

        // -------------------- Debug_alpha_V_value
        temp_param_var = Emission_model_element->FirstChildElement("Debug_alpha_V_value");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Debug_alpha_V_value value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Debug_alpha_V_value = std::stod(temp_param_var->GetText());

        // -------------------- Debug_rho_I_value
        temp_param_var = Emission_model_element->FirstChildElement("Debug_rho_I_value");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Debug_rho_I_value value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Debug_rho_I_value = std::stod(temp_param_var->GetText());

        // -------------------- Debug_rho_Q_value
        temp_param_var = Emission_model_element->FirstChildElement("Debug_rho_Q_value");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Debug_rho_Q_value value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Debug_rho_Q_value = std::stod(temp_param_var->GetText());

        // -------------------- Debug_rho_U_value
        temp_param_var = Emission_model_element->FirstChildElement("Debug_rho_U_value");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Debug_rho_U_value value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Debug_rho_U_value = std::stod(temp_param_var->GetText());

        // -------------------- Debug_rho_V_value
        temp_param_var = Emission_model_element->FirstChildElement("Debug_rho_V_value");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the Debug_rho_V_value value!" << "\n"; return ERROR; }
        p_Init_conditions->Emission_params.Debug_rho_V_value = std::stod(temp_param_var->GetText());

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
    temp_param_var = Observer_element->FirstChildElement("Init_time");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observer init_time!" << "\n"; return ERROR; }
    Observer_params->init_time = std::stod(temp_param_var->GetText());

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
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observer azimuth!" << "\n"; return ERROR; }
    Observer_params->azimuth = std::stod(temp_param_var->GetText());

    // -------------------- Camera Rotation Angle
    temp_param_var = Observer_element->FirstChildElement("Cam_rotation_angle");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the camera rotation angle!" << "\n"; return ERROR; }
    Observer_params->cam_rotation_angle = std::stod(temp_param_var->GetText());

    // -------------------- Image linear Y min
    temp_param_var = Observer_element->FirstChildElement("Image_y_min");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observation window linear Y min!" << "\n"; return ERROR; }
    Observer_params->y_min = std::stod(temp_param_var->GetText());

    // -------------------- Image linear Y max
    temp_param_var = Observer_element->FirstChildElement("Image_y_max");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observation window linear Y max!" << "\n"; return ERROR; }
    Observer_params->y_max = std::stod(temp_param_var->GetText());

    // -------------------- Image linear X min
    temp_param_var = Observer_element->FirstChildElement("Image_x_min");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observation window linear X min!" << "\n"; return ERROR; }
    Observer_params->x_min = std::stod(temp_param_var->GetText());

    // -------------------- Image linear X max
    temp_param_var = Observer_element->FirstChildElement("Image_x_max");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observation window linear X max!" << "\n"; return ERROR; }
    Observer_params->x_max = std::stod(temp_param_var->GetText());

    // -------------------- Image angular Y min
    temp_param_var = Observer_element->FirstChildElement("Image_y_angle_min");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observation window angular Y min!" << "\n"; return ERROR; }
    Observer_params->y_angle_min = std::stod(temp_param_var->GetText());

    // -------------------- Image angular Y max
    temp_param_var = Observer_element->FirstChildElement("Image_y_angle_max");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observation window angular Y max!" << "\n"; return ERROR; }
    Observer_params->y_angle_max = std::stod(temp_param_var->GetText());

    // -------------------- Image angular X min
    temp_param_var = Observer_element->FirstChildElement("Image_x_angle_min");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observation window angular X min!" << "\n"; return ERROR; }
    Observer_params->x_angle_min = std::stod(temp_param_var->GetText());

    // -------------------- Image angular X max
    temp_param_var = Observer_element->FirstChildElement("Image_x_angle_max");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse observation window angular X max!" << "\n"; return ERROR; }
    Observer_params->x_angle_max = std::stod(temp_param_var->GetText());

    // -------------------- Use angular coords flag
    temp_param_var = Observer_element->FirstChildElement("Use_angular_coords");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the \"use angular coords\" flag!" << "\n"; return ERROR; }
    Observer_params->Use_angular_coords = bool(std::stoi(temp_param_var->GetText()));

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

Return_Values static parse_numerical_metric_XML(tinyxml2::XMLElement* Spline_XML, Metric_parameters_type* Metric_params) {

    tinyxml2::XMLElement* F0_control_vector_element;
    tinyxml2::XMLElement* F1_control_vector_element;
    tinyxml2::XMLElement* F2_control_vector_element;
    tinyxml2::XMLElement* W_control_vector_element;

    tinyxml2::XMLElement* Compactified_radial_grid_element;
    tinyxml2::XMLElement* Compactified_radial_grid_control_vector_element;
    tinyxml2::XMLElement* Theta_grid_element;
    tinyxml2::XMLElement* Theta_grid_control_vector_element;

    tinyxml2::XMLElement* temp_param_var;

    // -------------------- The F_0 control vector

    F0_control_vector_element = Spline_XML->FirstChildElement("F_0")->FirstChildElement("Control_vector");
    if (F0_control_vector_element == nullptr) { std::cout << "Failed to parse the numerical F_0 potential control vector node!" << "\n"; return ERROR; }

    /* Parse the length of the controll vector and allocate an array to hold it. */
    int Control_vector_size = std::stoi(F0_control_vector_element->Attribute("Component_number"));
    Metric_params->Numerical_metric_params.F_0_control_vector = new double[Control_vector_size];
    Metric_params->Numerical_metric_params.Control_vector_size = Control_vector_size;

    for (int idx = 0; idx <= Control_vector_size - 1; idx++) {

        temp_param_var = F0_control_vector_element->FirstChildElement(static_cast<const char*>(("Component_idx_" + std::to_string(idx)).c_str()));
        if (temp_param_var == nullptr) { std::cout << std::format("Failed to parse the numerical F_0 potential control vector component at idx {}! \n", idx) ; return ERROR; }

        Metric_params->Numerical_metric_params.F_0_control_vector[idx] = std::stod(temp_param_var->GetText());

    }

    // -------------------- The g_rr control vector

    F1_control_vector_element = Spline_XML->FirstChildElement("F_1")->FirstChildElement("Control_vector");
    if (F1_control_vector_element == nullptr) { std::cout << "Failed to parse the numerical F_1 potential control vector node!" << "\n"; return ERROR; }

    /* Parse the length of the controll vector and allocate an array to hold it. */
    Control_vector_size = std::stoi(F1_control_vector_element->Attribute("Component_number"));
    Metric_params->Numerical_metric_params.F_1_control_vector = new double[Control_vector_size];

    for (int idx = 0; idx <= Control_vector_size - 1; idx++) {

        temp_param_var = F1_control_vector_element->FirstChildElement(static_cast<const char*>(("Component_idx_" + std::to_string(idx)).c_str()));
        if (temp_param_var == nullptr) { std::cout << std::format("Failed to parse the numerical F_1 potential control vector component at idx {}! \n", idx); return ERROR; }

        Metric_params->Numerical_metric_params.F_1_control_vector[idx] = std::stod(temp_param_var->GetText());

    }

    // -------------------- The g_thth control vector

    F2_control_vector_element = Spline_XML->FirstChildElement("F_2")->FirstChildElement("Control_vector");
    if (F2_control_vector_element == nullptr) { std::cout << "Failed to parse the numerical F_2 potential control vector node!" << "\n"; return ERROR; }

    /* Parse the length of the controll vector and allocate an array to hold it. */
    Control_vector_size = std::stoi(F2_control_vector_element->Attribute("Component_number"));
    Metric_params->Numerical_metric_params.F_2_control_vector = new double[Control_vector_size];

    for (int idx = 0; idx <= Control_vector_size - 1; idx++) {

        temp_param_var = F2_control_vector_element->FirstChildElement(static_cast<const char*>(("Component_idx_" + std::to_string(idx)).c_str()));
        if (temp_param_var == nullptr) { std::cout << std::format("Failed to parse the numerical F_2 potential control vector component at idx {}! \n", idx); return ERROR; }

        Metric_params->Numerical_metric_params.F_2_control_vector[idx] = std::stod(temp_param_var->GetText());

    }

    // -------------------- The g_phiphi control vector

    W_control_vector_element = Spline_XML->FirstChildElement("W")->FirstChildElement("Control_vector");
    if (W_control_vector_element == nullptr) { std::cout << "Failed to parse the numerical W potential control vector node!" << "\n"; return ERROR; }

    /* Parse the length of the controll vector and allocate an array to hold it. */
    Control_vector_size = std::stoi(W_control_vector_element->Attribute("Component_number"));
    Metric_params->Numerical_metric_params.W_control_vector = new double[Control_vector_size];

    for (int idx = 0; idx <= Control_vector_size - 1; idx++) {

        temp_param_var = W_control_vector_element->FirstChildElement(static_cast<const char*>(("Component_idx_" + std::to_string(idx)).c_str()));
        if (temp_param_var == nullptr) { std::cout << std::format("Failed to parse the numerical W potential control vector component at idx {}! \n", idx); return ERROR; }

        Metric_params->Numerical_metric_params.W_control_vector[idx] = std::stod(temp_param_var->GetText());

    }

    // -------------------- The compactified radial coordinate grid 

    Compactified_radial_grid_element = Spline_XML->FirstChildElement("Coordinate_grid")->FirstChildElement("Compactified_radial_coordinate_grid")->FirstChildElement("Grid_knots");
    if (Compactified_radial_grid_element == nullptr) { std::cout << "Failed to parse the compactified radial coordinate grid knots node!" << "\n"; return ERROR; }

    /* Parse the length of the controll vector and allocate an array to hold it. */
    int Grid_size = std::stoi(Compactified_radial_grid_element->Attribute("Grid_size"));
    Metric_params->Numerical_metric_params.Compactified_radial_grid = new double[Grid_size];
    Metric_params->Numerical_metric_params.Radial_grid_size = Grid_size;


    for (int idx = 0; idx <= Grid_size - 1; idx++) {

        temp_param_var = Compactified_radial_grid_element->FirstChildElement(static_cast<const char*>(("Grid_point_idx_" + std::to_string(idx)).c_str()));
        if (temp_param_var == nullptr) { std::cout << std::format("Failed to parse the compactified radial coordinate knot component at idx {}! \n", idx); return ERROR; }

        Metric_params->Numerical_metric_params.Compactified_radial_grid[idx] = std::stod(temp_param_var->GetText());

    }

    // -------------------- The compactified radial coordinate control vector - not strictly needed, 
    //                      but I use it as a sanity check to make sure the spline does not break somewhere 

    Compactified_radial_grid_control_vector_element = Spline_XML->FirstChildElement("Coordinate_grid")->FirstChildElement("Compactified_radial_coordinate_grid")->FirstChildElement("Control_vector");
    if (Compactified_radial_grid_control_vector_element == nullptr) { std::cout << "Failed to parse the compactified radial coordinate grid knots node!" << "\n"; return ERROR; }

    /* Parse the length of the controll vector and allocate an array to hold it. */
    Control_vector_size = std::stoi(Compactified_radial_grid_control_vector_element->Attribute("Component_number"));
    Metric_params->Numerical_metric_params.Compactified_radial_grid_control_vector = new double[Control_vector_size];

    for (int idx = 0; idx <= Control_vector_size - 1; idx++) {

        temp_param_var = Compactified_radial_grid_control_vector_element->FirstChildElement(static_cast<const char*>(("Component_idx_" + std::to_string(idx)).c_str()));
        if (temp_param_var == nullptr) { std::cout << std::format("Failed to parse the compactified radial coordinate control vector component at idx {}! \n", idx); return ERROR; }

        Metric_params->Numerical_metric_params.Compactified_radial_grid_control_vector[idx] = std::stod(temp_param_var->GetText());

    }

    // -------------------- The theta grid 

    Theta_grid_element = Spline_XML->FirstChildElement("Coordinate_grid")->FirstChildElement("Theta_coordinate_grid")->FirstChildElement("Grid_knots");
    if (Theta_grid_element == nullptr) { std::cout << "Failed to parse the theta grid knots node!" << "\n"; return ERROR; }

    /* Parse the length of the controll vector and allocate an array to hold it. */
    Grid_size = std::stoi(Theta_grid_element->Attribute("Grid_size"));
    Metric_params->Numerical_metric_params.Theta_grid = new double[Grid_size];
    Metric_params->Numerical_metric_params.Theta_grid_size = Grid_size;

    for (int idx = 0; idx <= Grid_size - 1; idx++) {

        temp_param_var = Theta_grid_element->FirstChildElement(static_cast<const char*>(("Grid_point_idx_" + std::to_string(idx)).c_str()));
        if (temp_param_var == nullptr) { std::cout << std::format("Failed to parse the theta knot component at idx {}! \n", idx); return ERROR; }

        Metric_params->Numerical_metric_params.Theta_grid[idx] = std::stod(temp_param_var->GetText());

    }

    // -------------------- The theta control vector - not strictly needed, 
    //                      but I use it as a sanity check to make sure the spline does not break somewhere 

    Theta_grid_control_vector_element = Spline_XML->FirstChildElement("Coordinate_grid")->FirstChildElement("Theta_coordinate_grid")->FirstChildElement("Control_vector");
    if (Theta_grid_control_vector_element == nullptr) { std::cout << "Failed to parse the theta coordinate grid knots node!" << "\n"; return ERROR; }

    /* Parse the length of the controll vector and allocate an array to hold it. */
    Control_vector_size = std::stoi(Theta_grid_control_vector_element->Attribute("Component_number"));
    Metric_params->Numerical_metric_params.Theta_grid_control_vector = new double[Control_vector_size];

    for (int idx = 0; idx <= Control_vector_size - 1; idx++) {

        temp_param_var = Theta_grid_control_vector_element->FirstChildElement(static_cast<const char*>(("Component_idx_" + std::to_string(idx)).c_str()));
        if (temp_param_var == nullptr) { std::cout << std::format("Failed to parse the theta coordinate control vector component at idx {}! \n", idx); return ERROR; }

        Metric_params->Numerical_metric_params.Theta_grid_control_vector[idx] = std::stod(temp_param_var->GetText());

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

    }else if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Numerical")) {

        Metric_params->e_Spacetime = Numerical;

        temp_param_var = Metric_element->FirstChildElement("Numerical_metric_spline_path");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the numerical metric XML file path!" << "\n"; return ERROR; }
        Metric_params->Numerical_metric_params.Metric_file_path = temp_param_var->GetText();

        tinyxml2::XMLDocument Spline_XML;
        tinyxml2::XMLError e_parse_result = Spline_XML.LoadFile(temp_param_var->GetText());
        if (e_parse_result != tinyxml2::XML_SUCCESS) { return ERROR; }

        temp_param_var = Spline_XML.FirstChildElement("Metric_spline_coefficients");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the metric spline coefficients node!" << "\n"; return ERROR; }

        if (OK != parse_numerical_metric_XML(temp_param_var, Metric_params)) { std::cout << "Failed to parse the metric XML!" << "\n"; return ERROR; }

        temp_param_var = Metric_element->FirstChildElement("Horizon_radius");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the horizon radius!" << "\n"; return ERROR; }
        Metric_params->Numerical_metric_params.Horizon_radius = std::stod(temp_param_var->GetText());

        temp_param_var = Metric_element->FirstChildElement("ADM_Mass");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the ADM mass!" << "\n"; return ERROR; }
        Metric_params->Numerical_metric_params.M_ADM = std::stod(temp_param_var->GetText());

        temp_param_var = Metric_element->FirstChildElement("ADM_ang_momentum");
        if (temp_param_var == nullptr) { std::cout << "Failed to parse the ADM angular momentum!" << "\n"; return ERROR; }
        Metric_params->Numerical_metric_params.a_ADM = std::stod(temp_param_var->GetText());

        std::string Anzatz_type = Metric_element->FirstChildElement("Numerical_metric_anzatz_type")->GetText();

        if (0 == strcmp(static_cast<const char*>(Anzatz_type.c_str()), "Anzatz_1")) {

            Metric_params->Numerical_metric_params.e_Anzatz = e_Anzatz_1;

        }
        else if (0 == strcmp(static_cast<const char*>(Anzatz_type.c_str()), "Anzatz_2")) {

            Metric_params->Numerical_metric_params.e_Anzatz = e_Anzatz_2;

        }
        else { std::cout << "Unsuppored metric anzatz type! \n"; return ERROR; }

    }else if (0 == strcmp(static_cast<const char*>(Metric_type.c_str()), "Minkowski")) {

        Metric_params->e_Spacetime = Minkowski;

    }
    else { std::cout << "Unsupported metric type! \n"; return ERROR; }

    temp_param_var = Metric_element->FirstChildElement("Scattering_radius");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the scattering radius!" << "\n"; return ERROR; }
    Metric_params->Scattering_radius = std::stod(temp_param_var->GetText());

    temp_param_var = Metric_element->FirstChildElement("Distance_to_singular_point");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the distance to the singular point!" << "\n"; return ERROR; }
    Metric_params->Min_distance_to_singular_point = std::stod(temp_param_var->GetText());

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

    tinyxml2::XMLElement* temp_param_var;

    /* ====================================== Parse the average pitch angle flag and sample number ====================================== */

    temp_param_var = Root_node->FirstChildElement("Average_emission_pitch_angle");
    if (temp_param_var == nullptr) { std::cout << "Failed to find the pitch angle averaging flag!" << "\n"; return ERROR; }
    p_Initial_conditions->Average_electron_pitch_angle = std::stoi(temp_param_var->GetText());

    temp_param_var = Root_node->FirstChildElement("Emission_pitch_angle_samples_to_average");
    if (temp_param_var == nullptr) { std::cout << "Failed to find the number of pitch angle samples to average!" << "\n"; return ERROR; }
    p_Initial_conditions->Emission_pitch_angle_samples_to_average = std::stoi(temp_param_var->GetText());

    temp_param_var = Root_node->FirstChildElement("Thermalize_emission_medium");
    if (temp_param_var == nullptr) { std::cout << "Failed to parse the Thermalize_emission_medium flag!" << "\n"; return ERROR; }
    p_Initial_conditions->Thermalize_emission_medium = std::stoi(temp_param_var->GetText());

    /* ====================================== Parse the simulation mode specific settings ====================================== */

    temp_param_var = Root_node->FirstChildElement("Simulation_mode");
    if (temp_param_var == nullptr) { std::cout << "Failed to find the simulation mode!" << "\n"; return ERROR; }
    p_Initial_conditions->Simulation_mode = static_cast<Simulation_mode_enums>(std::stoi(temp_param_var->GetText()));

    temp_param_var = Root_node->FirstChildElement("Sim_mode_2_param_value_number");
    if (temp_param_var == nullptr) { std::cout << "Failed to find the simulation mode 2 number of metric parameter values!" << "\n"; return ERROR; }
    p_Initial_conditions->Sim_mode_2_param_value_number = std::stoi(temp_param_var->GetText());

    temp_param_var = Root_node->FirstChildElement("Sim_mode_3_X_init");
    if (temp_param_var == nullptr) { std::cout << "Failed to find sim mode 3 X init!" << "\n"; return ERROR; }
    p_Initial_conditions->Sim_mode_3_X_init = std::stod(temp_param_var->GetText());

    temp_param_var = Root_node->FirstChildElement("Sim_mode_3_Y_init");
    if (temp_param_var == nullptr) { std::cout << "Failed to find sim mode 3 Y init!" << "\n"; return ERROR; }
    p_Initial_conditions->Sim_mode_3_Y_init = std::stod(temp_param_var->GetText());

    temp_param_var = Root_node->FirstChildElement("Max_image_order");
    if (temp_param_var == nullptr) { std::cout << "Failed to find max image order!" << "\n"; return ERROR; }
    p_Initial_conditions->Max_order = std::stoi(temp_param_var->GetText());

    /* ====================================== Parse the central object mass ====================================== */

    temp_param_var = Root_node->FirstChildElement("Central_object_mass");
    if (temp_param_var == nullptr) { std::cout << "Failed to find the central object mass!" << "\n"; return ERROR; }
    p_Initial_conditions->central_object_mass = std::stod(temp_param_var->GetText());

    /* ====================================== Parse the observer parameters ====================================== */

    temp_param_var = Root_node->FirstChildElement("Observer");
    if (temp_param_var == nullptr) { std::cout << "Failed to find the Observer node!" << "\n"; return ERROR; }
    if (OK != parse_observer_parameters(temp_param_var, &p_Initial_conditions->Observer_params)) { return ERROR; }

    /* ====================================== Parse the metric parameters ====================================== */

    temp_param_var = Root_node->FirstChildElement("Metric");
    if (temp_param_var == nullptr) { std::cout << "Failed to find the Metric node!" << "\n"; return ERROR; }
    if (OK != parse_metric_parameters(temp_param_var, &p_Initial_conditions->Metric_parameters)) { return ERROR; };

    /* ====================================== Parse the integrator parameters ====================================== */

    temp_param_var = Root_node->FirstChildElement("Integrator");
    if (temp_param_var == nullptr) { std::cout << "Failed to find the integrator node!" << "\n"; return ERROR; }
    if (OK != parse_integrator_params(temp_param_var, &p_Initial_conditions->Integrator_params)) { return ERROR; };

    /* ====================================== Parse the accretion disk parameters ====================================== */

    temp_param_var = Root_node->FirstChildElement("Accretion_Disk");
    if (temp_param_var == nullptr) { std::cout << "Failed to find the Accretion Disk node!" << "\n"; return ERROR; }
    if (OK != parse_disk_params(temp_param_var, &p_Initial_conditions->Disk_params)) { return ERROR; };

    /* ====================================== Parse the hotspot parameters ====================================== */

    temp_param_var = Root_node->FirstChildElement("Hotspot");
    if (temp_param_var == nullptr) { std::cout << "Failed to find the Hotspot node!" << "\n"; return ERROR; }
    if (OK != parse_hotspot_params(temp_param_var, &p_Initial_conditions->Hotspot_params)) { return ERROR; };

    /* ====================================== Parse the emission model parameters ====================================== */

    temp_param_var = Root_node->FirstChildElement("Emission_models");
    if (temp_param_var == nullptr) { std::cout << "Failed to find the Hotspot node!" << "\n"; return ERROR; }
    if (OK != parse_emission_model_params(temp_param_var, p_Initial_conditions)) { return ERROR; };

    /* ====================================== Parse the file paths ====================================== */

    temp_param_var = Root_node->FirstChildElement("File_Manager");
    if (temp_param_var == nullptr) { std::cout << "Failed to find the File paths node!" << "\n"; return ERROR; }
    if (OK != parse_file_manager_params(temp_param_var, &p_Initial_conditions->File_manager_params)) { return ERROR; };

    return OK;

}