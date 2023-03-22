#pragma once

#ifndef SIM_MODES

	#define SIM_MODES

	#include "Constants.h"
	#include <iostream>

	void run_simulation_mode_1(Initial_conditions_type* s_Initial_Conditions);

	void run_simulation_mode_2(Initial_conditions_type* s_Initial_Conditions);

	void run_simulation_mode_3(Initial_conditions_type* s_Initial_Conditions);

	void print_progress(int current, int max, bool lens_from_file, bool Normalizing_colormap);

#endif

