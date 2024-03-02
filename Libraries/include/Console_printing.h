#ifndef CONSOLE_PRINTER
    
    #include "Constants.h"
    #include <string>
    #include <iostream>

	#define CONSOLE_PRINTER


    class Console_Printer_class {

    private:

        std::string Metric_strings[SPACETIME_NUMBER] = {

            "Kerr Black Hole",
            "Teo Wormhole",
            "Regular Black Hole",
            "Janis-Newman-Winicour Naked Singularity",
            "Gauss-Bonnet Naked Singularity",
            "Black Hole With A Dark Matter Halo"
        };

        std::string Emission_model_strings[2] = {

            "Thermally Averaged Synchotron Emission",
            "Phenomenological Model Of The M87* Emission"

        };

        std::string Disk_model_strings[2] = {

            "Power Law Of The Form: Density ~ 1 / r^2",
            "Exponential Law Of The Form: Density ~ exp(-r^2 / radial_scale^2)"

        };

    public:

        void print_ASCII_art();

        void print_sim_parameters();

	};


#endif