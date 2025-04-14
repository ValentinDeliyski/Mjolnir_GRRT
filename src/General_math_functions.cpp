#include "General_math_functions.h"

double vector_norm(const double* const Vector, const int Vector_size) {

	double norm{};

	for (int index = 0; index <= Vector_size - 1; index++) {

		norm += Vector[index] * Vector[index];

	}

	norm = sqrt(norm);

	return norm;

};

void mat_vec_multiply_4D(double const Matrix[4][4], const double* const Vector, double* const Result) {

	for (int row = 0; row <= 3; row += 1) {

		Result[row] = 0.0;

		for (int column = 0; column <= 3; column += 1) {

			Result[row] += Matrix[row][column] * Vector[column];

		}
	}

};

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
						  double* const Crossing_coords, 
						  double* const Crossing_momenta) {

	// Check weather the equator has been crossed
	if ((State_Vector[e_theta] - M_PI_2) * (Old_State_Vector[e_theta] - M_PI_2) > 0)
	{

		memset(Crossing_coords, 0, 3 * sizeof(double));
		memset(Crossing_momenta, 0, 3 * sizeof(double));

		return false;

	}

	/*
	
	Interpolate the equatorial crossing coorinates
	
	*/

	double x = State_Vector[e_r] * sin(State_Vector[e_theta]) * cos(State_Vector[e_phi]);
	double y = State_Vector[e_r] * sin(State_Vector[e_theta]) * sin(State_Vector[e_phi]);
	double z = State_Vector[e_r] * cos(State_Vector[e_theta]);

	double x_old = Old_State_Vector[e_r] * sin(Old_State_Vector[e_theta]) * cos(Old_State_Vector[e_phi]);
	double y_old = Old_State_Vector[e_r] * sin(Old_State_Vector[e_theta]) * sin(Old_State_Vector[e_phi]);
	double z_old = Old_State_Vector[e_r] * cos(Old_State_Vector[e_theta]);

	double gradient[3]   = { x - x_old, y - y_old, z - z_old };
	double const_term[3] = { x_old, y_old, z_old };

	double crossing_param = -const_term[2] / gradient[2];

	for (int index = 0; index <= 2; index++) {

		Crossing_coords[index] = gradient[index] * crossing_param + const_term[index];

	}

	/*
	
	Interpolate the covariant momenta at those coorinates
	
	*/

	double momentum_param = (M_PI_2 - Old_State_Vector[e_theta]) / (State_Vector[e_theta] - Old_State_Vector[e_theta]);

	for (int index = e_t; index <= e_phi; index++) {

		Crossing_momenta[index] = momentum_param * State_Vector[e_p_t + index] + (1 - momentum_param) * Old_State_Vector[e_p_t + index];
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

void convert_spherical_to_cartesian(double* Spherical_Coords, double* Cartesian_Coords) {

	double sin_theta = sin(Spherical_Coords[e_theta]);
	double cos_theta = cos(Spherical_Coords[e_theta]);

	double sin_phi = sin(Spherical_Coords[e_phi]);
	double cos_phi = cos(Spherical_Coords[e_phi]);

	Cartesian_Coords[x] = Spherical_Coords[e_r] * sin_theta * cos_phi;
	Cartesian_Coords[y] = Spherical_Coords[e_r] * sin_theta * sin_phi;
	Cartesian_Coords[z] = Spherical_Coords[e_r] * cos_theta;

}

void convert_cartesian_to_spherical(double* Cartesian_Coords, double* Spherical_Coords) {

	Spherical_Coords[e_r]  = Cartesian_Coords[x] * Cartesian_Coords[x];
	Spherical_Coords[e_r] += Cartesian_Coords[y] * Cartesian_Coords[y];
	Spherical_Coords[e_r] += Cartesian_Coords[z] * Cartesian_Coords[z];
	Spherical_Coords[e_r]  = sqrt(Spherical_Coords[e_r]);

	Spherical_Coords[e_theta] = acos(Cartesian_Coords[z] / Spherical_Coords[e_r]);

	Spherical_Coords[e_phi] = atan2(Cartesian_Coords[y], Cartesian_Coords[x]);

}

void add_4D_vectors(const double* const Vec_1, const double* const Vec_2, double* const Result) {

	for (int idx = 0; idx <= 3; idx++) {

		Result[idx] = Vec_1[idx] + Vec_2[idx];

	}

}

double int_power(const double base, const int exponent) {

	double result = 1.0;

	for (int counter = exponent; counter > 0; counter--) { result *= base; }

	return result;

}