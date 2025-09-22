#pragma once

struct Simulation_Context_type;
struct Results_type;

void run_image_generation(const Simulation_Context_type* const p_Sim_Context, Results_type* const p_Ray_results);

void run_geodesic_sweep(const Simulation_Context_type* const p_Sim_Context, Results_type* const p_Ray_results);

void make_geodesic_log(const Simulation_Context_type* const p_Sim_Context, Results_type* const p_Ray_results);


