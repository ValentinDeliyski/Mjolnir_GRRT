#include "Console_printing.h"
#include "Spacetimes.h"

void Console_Printer_class::print_ASCII_art() {

        std::cout <<

            " ######   ########     ###    ##     ## #### ########    ###    ######## ####  #######  ##    ##    ###    ##          ########     ###    ##    ##    ######## ########     ###     ######  ######## ########  \n"
            "##    ##  ##     ##   ## ##   ##     ##  ##     ##      ## ##      ##     ##  ##     ## ###   ##   ## ##   ##          ##     ##   ## ##    ##  ##        ##    ##     ##   ## ##   ##    ## ##       ##     ## \n"
            "##        ##     ##  ##   ##  ##     ##  ##     ##     ##   ##     ##     ##  ##     ## ####  ##  ##   ##  ##          ##     ##  ##   ##    ####         ##    ##     ##  ##   ##  ##       ##       ##     ## \n"
            "##   #### ########  ##     ## ##     ##  ##     ##    ##     ##    ##     ##  ##     ## ## ## ## ##     ## ##          ########  ##     ##    ##          ##    ########  ##     ## ##       ######   ########  \n"
            "##    ##  ##   ##   #########  ##   ##   ##     ##    #########    ##     ##  ##     ## ##  #### ######### ##          ##   ##   #########    ##          ##    ##   ##   ######### ##       ##       ##   ##   \n"
            "##    ##  ##    ##  ##     ##   ## ##    ##     ##    ##     ##    ##     ##  ##     ## ##   ### ##     ## ##          ##    ##  ##     ##    ##          ##    ##    ##  ##     ## ##    ## ##       ##    ##  \n"
            " ######   ##     ## ##     ##    ###    ####    ##    ##     ##    ##    ####  #######  ##    ## ##     ## ########    ##     ## ##     ##    ##          ##    ##     ## ##     ##  ######  ######## ##     ## \n";

        std::cout << '\n' << '\n';

}

void Console_Printer_class::print_sim_parameters() {


    std::cout << "============================================================ SIMULATION METADATA ============================================================"
              << "\n"
              << "Spacetime: "
              << Metric_strings[e_metric]
              << "\n";

    switch (e_metric) {

    case Kerr:

        std::cout << "Spin Parameter = " << SPIN << '\n';

        break;

    case Wormhole:

        std::cout << "Redshit Parameter = " << WH_REDSHIFT << '\n'
            << "Spin Parameter = " << SPIN << '\n';
        break;

    case Reg_Black_Hole:

        std::cout << "Metric Parameter = " << RBH_PARAM << '\n';

        break;

    case Naked_Singularity:

        std::cout << "Metric Parameter = " << JNW_GAMMA << '\n';

        break;

    case Gauss_Bonnet:

        std::cout << "Metric Parameter = " << GAUSS_BONNET_GAMMA << '\n';

        break;

    case BH_w_Dark_Matter:

        std::cout << "Compactness = " << COMPACTNESS << '\n'
            << "Halo Mass = " << M_HALO << '\n';
        break;

    default:

        std::cout << "Unknown Metric Enum - Something Broke!" << '\n';

    };

    std::cout << "Observer Radial Position [M]: " << r_obs << '\n'
              << "Observer Inclination [Deg]: " << theta_obs * 180 / M_PI << '\n'
              << "Observation Frequency [Hz]: " << OBS_FREQUENCY_CGS
              << '\n' << '\n'
              << "------------------------------------------------------- Optically Thin Disk Metadata -------------------------------------------------------"
              << "\n"
              << "Active disk profile: "
              << this->Disk_model_strings[e_disk_model] << '\n';

    switch (e_disk_model) {

    case Power_law:

        std::cout << "Disk Opening Angle [tan(angle)]: "
            << DISK_OPENING_ANGLE
            << "\n"
            << "R_0 [M]: "
            << DISK_R_0
            << "\n"
            << "R_Cutoff [M]: ";

        std::cout << DISK_R_Cutoff << "\n";

        std::cout << "Cutoff Scale [M]: " << DISK_CUTOFF_SCALE << "\n";

        break;

    case Exponential_law:

        std::cout << "Exponential law (rho ~ exp(-r^2))"
            << "\n"
            << "Disk Height Scale [M]: "
            << DISK_HEIGHT_SCALE
            << "\n"
            << "Disk Radial Scale [M]: "
            << DISK_RADIAL_SCALE
            << "\n";

        break;

    }

    std::cout << '\n';
    std::cout << "Active Emission Function: "
              << this->Emission_model_strings[e_disk_emission] << '\n';

    switch (e_disk_emission) {

    case Phenomenological_synchotron:

        std::cout << "Emission Power Law Exponent [-]: "
                  << EMISSION_POWER_LAW
                  << "\n"
                  << "Absorbtion Coefficient [-]: "
                  << DISK_ABSORBTION_COEFF
                  << "\n"
                  << "Source Function Power Law Exponent [-]: "
                  << SOURCE_F_POWER_LAW
                  << "\n"
                  << "Emission Scale [Jy / sRad]: "
                  << EMISSION_SCALE_PHENOMENOLOGICAL
                  << "\n";

        break;

    case Thermal_dist_synchotron:

        std::cout << "Disk Magnetization [-]: "
                  << DISK_MAGNETIZATION
                  << "\n"
                  << "Max Electron Temperature [K]: "
                  << DISK_T_ELECTRON_SCALE_CGS
                  << "\n"
                  << "Max Disk Density [g / cm^3]: "
                  << DISK_N_ELECTRON_SCALE_CGS
                  << "\n";

        break;

    }

    std::cout << '\n';

    std::cout << "------------------------------------------------------------- Hotpost Metadata -------------------------------------------------------------"
             << "\n"
             << "Max Hotspot Density[g / cm ^ 3] "
             << HOTSPOT_N_ELECTRON_SLACE_CGS
             << "\n"
             << "Hotspot Radial Scale [M]: "
             << HOTSPOT_SPREAD
             << "\n"
             << "Hotspot r_center [M]: "
             << HOTSPOT_R_COORD
             << "\n";

    std::cout << '\n';

    std::cout << "Simulation Resolution: " << NUM_RAYS_X << "X" << NUM_RAYS_Y << '\n';

}



