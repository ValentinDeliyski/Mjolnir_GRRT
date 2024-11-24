#pragma once
#define _USE_MATH_DEFINES
#include "Enumerations.h"
#include "Spacetimes.h"
#include "Constants.h"
#include "Structs.h"
#include <iostream>
#include <string>
#include <cmath>

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

public:

    void print_ASCII_art();

    void print_sim_parameters(Initial_conditions_type* p_Initial_Conditions);;

};