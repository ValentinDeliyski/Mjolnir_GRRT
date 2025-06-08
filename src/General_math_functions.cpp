#include "General_math_functions.h"

double vector_norm(const double* const Vector, const int Vector_size) {

	double norm{};

	for (int index = 0; index <= Vector_size - 1; index++) {

		norm += Vector[index] * Vector[index];

	}
	
	norm = sqrt(norm);

	return norm;

};

void matrix_matrix_multiply(const double Mat_A[4][4], const double Mat_B[4][4], double Result[4][4]) {

	memset(Result, 0.0, 16 * sizeof(double));

	for (int row = 0; row < 4; row++) {

		for (int column = 0; column < 4; column++){

			for (int k = 0; k < 4; k++) {

				Result[row][column] += Mat_A[row][k] * Mat_B[k][column];

			}
		}
	}
}

void mat_vec_multiply_4D(double const Matrix[4][4], const double* const Vector, double* const Result) {

	memset(Result, 0.0, 4 * sizeof(double));

	for (int row = 0; row <= 3; row += 1) {

		Result[row] = 0.0;

		for (int column = 0; column <= 3; column += 1) {

			Result[row] += Matrix[row][column] * Vector[column];

		}
	}

}

double get_max_element(const double* const Vector, int const Element_number) {

	double max = fabs(Vector[e_r]);
	double max_candidate = 0;

	for (int index = e_theta; index <= Element_number; index += 1) {

		max_candidate = fabs(Vector[index]);

		if (max_candidate > max) {
			
			max = max_candidate;
		
		}

	}

	return max;
}

double get_max_relative_error(const double* const Error_state, const double* const Current_state) {

	double max_rel_error = fabs(Error_state[e_t] / Current_state[e_t]);

	for (int index = e_r; index <= e_State_Number - 2; index += 1) {

		double temp_error = fabs(Error_state[index] / (Current_state[index]));

		if (temp_error > max_rel_error && !isinf(temp_error)) {

			max_rel_error = temp_error;

		}
	}

	return max_rel_error;
}

bool interpolate_crossing(const double* const State_Vector, 
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

double dot_product(const double* const Vector_1, const double* const Vector_2, int Vector_size) {

	double result{};

	for (int index = 0; index <= Vector_size - 1; index++) {

		result += Vector_1[index] * Vector_2[index];

	}

	return  result;

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