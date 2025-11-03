#include "General_math_functions.h"

double vector_norm(const double* const Vector, const int Vector_size) {

	double norm{};

	for (int index = 0; index < Vector_size ; index++) {

		norm += Vector[index] * Vector[index];

	}
	
	norm = std::sqrt(norm);

	return norm;

};

void matrix_matrix_multiply(const double Mat_A[4][4], const double Mat_B[4][4], double Result[4][4]) {

	/* Note to self: This memset is here for a reason, so do not remove it! */
	memset(Result, 0, 16 * sizeof(double));

	for (int row = 0; row < 4; row++) {

		for (int column = 0; column < 4; column++){

			Result[row][column] = 0.0;

			for (int k = 0; k < 4; k++) {

				Result[row][column] += Mat_A[row][k] * Mat_B[k][column];

			}
		}
	}
}

void mat_vec_multiply_4D(double const Matrix[4][4], const double* const Vector, double* const Result) {

	/* Note to self: This memset is here for a reason, so do not remove it! */
	memset(Result, 0, 4 * sizeof(double));

	for (int row = 0; row < 4; row++) {

		for (int column = 0; column < 4; column++) {

			Result[row] += Matrix[row][column] * Vector[column];

		}
	}

}

double get_max_element(const double* const Vector, int const Element_number) {

	double max = fabs(Vector[0]);
	double max_candidate = 0;

	for (int index = 1; index < Element_number; index++) {

		max_candidate = fabs(Vector[index]);

		if (max_candidate > max) {
			
			max = max_candidate;
		
		}

	}

	return max;
}

double get_max_relative_error(const double* const Error_state, const double* const Current_state) {

	double max_rel_error = fabs(Error_state[e_t] / Current_state[e_t]);

	if (isinf(max_rel_error) || isnan(max_rel_error)) { max_rel_error = 0.0; }

	for (int index = e_r; index < e_Dynamic_state_size; index += 1) {

		double temp_error = fabs(Error_state[index] / (Current_state[index]));

		if (temp_error > max_rel_error && !isinf(temp_error)) {

			max_rel_error = temp_error;

		}
	}

	return max_rel_error;
}

bool interpolate_equatorial_crossing(const double* const State_Vector, 
									 const double* const Old_State_Vector, 
									 double* const Crossing_State) {

	// Check weather the equator has been crossed
	if ((State_Vector[e_theta] - M_PI_2) * (Old_State_Vector[e_theta] - M_PI_2) > 0) { return false; }

	/* ---------- Interpolate the equatorial crossing coorinates between the two photon state vectors (including the coordinate time) ----------  */

	double Current_position_cartesian[3]{};
	double Previous_position_cartesian[3]{};
	double Crossing_state_cartesian[3]{};

	convert_spherical_to_cartesian(State_Vector, Current_position_cartesian);
	convert_spherical_to_cartesian(Old_State_Vector, Previous_position_cartesian);

	double Direction_vector[3]{};

	for (int idx = e_x; idx <= e_z; idx++) {
		
		Direction_vector[idx] = Current_position_cartesian[idx] - Previous_position_cartesian[idx]; 
	
	}

	const double crossing_param_position = -Previous_position_cartesian[e_z] / Direction_vector[e_z];

	for (int idx = e_x; idx <= e_z; idx++) {
		
		Crossing_state_cartesian[idx] = Previous_position_cartesian[idx] + crossing_param_position * Direction_vector[idx];
	
	}

	Crossing_State[e_t] = Old_State_Vector[e_t] + crossing_param_position * (State_Vector[e_t] - Old_State_Vector[e_t]);

	/* ---------- The crossing coordinates are currently in cartesian form -> convert them to spherical ---------- */

	convert_cartesian_to_spherical(Crossing_state_cartesian, Crossing_State);

	/* ---------- Interpolate the covariant momenta at those coorinates ---------- */

	const double crossing_param_momentum = (M_PI_2 - Old_State_Vector[e_theta]) / (State_Vector[e_theta] - Old_State_Vector[e_theta]);

	for (int index = e_p_t; index <= e_p_phi; index++) {

		Crossing_State[index] = crossing_param_momentum * State_Vector[index] + (1 - crossing_param_momentum) * Old_State_Vector[index];

	}

	return true;
}

void interpolate_celestial_sphere_crossing(const double* const Current_State_Vector_Spherical,
										   const double* const Old_State_Vector_Spherical,
										   const double Celestial_Sphere_Raius,
										   double* const Crossing_State_Spherical) {

	/* This should only happen to rays that either hit the central object or run out of integration iterations / affine parameter. */
	if (Current_State_Vector_Spherical[e_r] < Celestial_Sphere_Raius) {

		Crossing_State_Spherical[e_r] = Current_State_Vector_Spherical[e_r];
		Crossing_State_Spherical[e_theta] = 1e100;
		Crossing_State_Spherical[e_phi] = 1e100;

		return;

	}

	/* ----------- Get the current and old state in cartesian components ----------- */

	double Current_State_Cartesian[3]{}, Old_State_Cartesian[3]{};

	convert_spherical_to_cartesian(Current_State_Vector_Spherical, Current_State_Cartesian);
	convert_spherical_to_cartesian(Old_State_Vector_Spherical, Old_State_Cartesian);

	/* -- Construct the difference vector between the two and intersect it with the celestial sphere -- */

	const double Difference_Vector[3] = { Current_State_Cartesian[e_x] - Old_State_Cartesian[e_x],
										  Current_State_Cartesian[e_y] - Old_State_Cartesian[e_y],
										  Current_State_Cartesian[e_z] - Old_State_Cartesian[e_z] };

	/* The expression for the intersection of a line with a sphere is given here https ://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
	   NOTE: In our case the sphere has its center point at (0, 0, 0). */

	const double U_dot_O   = dot_product(Difference_Vector, Old_State_Cartesian, 3);
	const double U_squared = dot_product(Difference_Vector, Difference_Vector, 3);
	const double O_squared = dot_product(Old_State_Cartesian, Old_State_Cartesian, 3);

	const double Crossing_param = (-U_dot_O + sqrt(U_dot_O * U_dot_O - U_squared * (O_squared - Celestial_Sphere_Raius * Celestial_Sphere_Raius))) / U_squared;

	const double Crossing_State_Cartesian[3] = { Old_State_Cartesian[e_x] + Crossing_param * Difference_Vector[e_x],
												 Old_State_Cartesian[e_y] + Crossing_param * Difference_Vector[e_y],
												 Old_State_Cartesian[e_z] + Crossing_param * Difference_Vector[e_z] };

	convert_cartesian_to_spherical(Crossing_State_Cartesian, Crossing_State_Spherical);

}

double dot_product(const double* const Vector_1, const double* const Vector_2, int Vector_size) {

	double result{};

	for (int index = 0; index < Vector_size; index++) {

		result += Vector_1[index] * Vector_2[index];

	}

	return  result;

}

void cross_product(const double* const Vec_1, const double* const Vec_2, double* const Result) {

	Result[e_x] = Vec_1[e_y] * Vec_2[e_z] - Vec_1[e_z] * Vec_2[e_y];
	Result[e_y] = Vec_1[e_z] * Vec_2[e_x] - Vec_1[e_x] * Vec_2[e_z];
	Result[e_z] = Vec_1[e_x] * Vec_2[e_y] - Vec_1[e_y] * Vec_2[e_x];

}

void convert_spherical_to_cartesian(const double* const Spherical_Coords, double* const Cartesian_Coords) {

	double sin_theta = sin(Spherical_Coords[e_theta]);
	double cos_theta = cos(Spherical_Coords[e_theta]);

	double sin_phi = sin(Spherical_Coords[e_phi]);
	double cos_phi = cos(Spherical_Coords[e_phi]);

	Cartesian_Coords[e_x] = Spherical_Coords[e_r] * sin_theta * cos_phi;
	Cartesian_Coords[e_y] = Spherical_Coords[e_r] * sin_theta * sin_phi;
	Cartesian_Coords[e_z] = Spherical_Coords[e_r] * cos_theta;

}

void convert_cartesian_to_spherical(const double* const Cartesian_Coords, double* const Spherical_Coords) {

	Spherical_Coords[e_r]  = Cartesian_Coords[e_x] * Cartesian_Coords[e_x];
	Spherical_Coords[e_r] += Cartesian_Coords[e_y] * Cartesian_Coords[e_y];
	Spherical_Coords[e_r] += Cartesian_Coords[e_z] * Cartesian_Coords[e_z];
	Spherical_Coords[e_r]  = sqrt(Spherical_Coords[e_r]);

	Spherical_Coords[e_theta] = acos(Cartesian_Coords[e_z] / Spherical_Coords[e_r]);

	Spherical_Coords[e_phi] = atan2(Cartesian_Coords[e_y], Cartesian_Coords[e_x]);

}

void add_vectors(const double* const Vec_1, const double* const Vec_2, const int size, double* const Result) {

	for (int idx = 0; idx < size; idx++) {

		Result[idx] = Vec_1[idx] + Vec_2[idx];

	}

}

double int_power(const double base, const int exponent) {

	double result = 1.0;

	for (int counter = exponent; counter > 0; counter--) { result *= base; }

	return result;

}