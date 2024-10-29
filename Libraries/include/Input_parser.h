#pragma once
#include <iostream>
#include "Enumerations.h"
#include "Structs.h"
#include <filesystem>
#include <string>
#include <fstream>

Return_Values parse_simulation_input_XML(std::string input_file_path, Initial_conditions_type* Initial_conditions);
