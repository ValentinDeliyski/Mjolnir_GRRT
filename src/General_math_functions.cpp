#define _USE_MATH_DEFINES

#include "Disk_Models.h"
#include "Enumerations.h"

#include <cmath>

extern Novikov_Thorne_Model NT_Model;

double vector_norm(double Vector[], int Vector_size) {

	/**********************************************************************
	|                                                                     |
	|   @ Description: Computes the Euclidian norm of the vector Vector   |
	|                                                                     |
	|   @ Inputs:                                                         |
		  * Vector: Pointer to the vector, whose norm we will find		  |
	|	  * Vector_size: The dimention of the vector					  |
	|																	  |
	|   @ Ouput: Euclidian norm of the vector							  |
	|                                                                     |
	**********************************************************************/

	double norm{};

	for (int index = 0; index <= Vector_size - 1; index++) {

		norm += Vector[index] * Vector[index];

	}

	norm = sqrt(norm);

	return norm;

};

int mat_vec_multiply_4D(double Matrix[4][4], double Vector[4], double result[4]) {

	/******************************************************************************
	|                                                                             |
	|   @ Description: Multiplies the 4D vector Vector by the 4x4 matrix Matrix   |
	|                                                                             |
	|   @ Inputs:                                                                 |
	|     * Matrix: Pointer to the input matrix we are multiplying by			  |
	|	  * Vector: Pointer to the vector to be multiplied						  |
	|	  * Result: Pointer to the vector, which stores the result of the		  |
	|	  multiplication														  |
	|																			  |
	|   @ Ouput: Pointer to the result from the multiplication					  |
	|                                                                             |
	******************************************************************************/


	for (int row = 0; row <= 3; row += 1) {

		for (int column = 0; column <= 3; column += 1) {

			result[row] += Matrix[row][column] * Vector[column];

		}
	}

	return OK;

};

double my_max(double vector[]) {

	/*****************************************************************************
	|                                                                            |
	|   @ Description: Returns the largest by absolute value element of vector   |
	|                                                                            |
	|   @ Inputs:                                                                |
	|     * vector: Pointer to the input array									 |
	|																			 |
	|   @ Ouput: The largest (unsigned) element of vector by absolute value 	 |
	|                                                                            |
	*****************************************************************************/

	int index_max = e_State_Number - 1;

	for (int index = 0; index <= index_max; index += 1) {

		if (vector[index] < 0) {

			vector[index] = -1.0 * vector[index];
		}

	}

	double max = vector[0];

	for (int index = 1; index <= index_max; index += 1) {

		if (vector[index] > max) {

			max = vector[index];

		}
	}

	return max;
}

bool interpolate_crossing(double State_Vector[], double Old_State_Vector[], double Crossing_coords[], double crossing_momenta[]) {

	/***********************************************************************************************
	|                                                                                              |
	|   @ Description: Interpolates the equatorial crossing point from two state vecrors -	       |
	|	One at z < 0 and the other at z > 0														   |
	|                                                                                              |
	|   @ Inputs:                                                                                  |
	|     * State_Vector: Pointer to the current state vector									   |
	|     * Old_State_Vector: Pointer to the old state vector									   |
	|     * Crossing_coords: Pointer to the interpolated x, y coordinates of equatorial crossing   |
	|                                                                                              |
	|   @ Ouput: Boolean: Weather the equatorial crossing point is within the NT disk              |
	|                                                                                              |
	***********************************************************************************************/

	if (cos(State_Vector[e_theta]) * cos(Old_State_Vector[e_theta]) > 0)
	{

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

	for (int index = 0; index < 3; index++) {

		Crossing_coords[index] = gradient[index] * crossing_param + const_term[index];

	}

	double r_crossing_squared = Crossing_coords[0] * Crossing_coords[0] + Crossing_coords[1] * Crossing_coords[1];

	double crossing_coords_spherical[2]{};

		crossing_coords_spherical[e_r]     = sqrt(r_crossing_squared);
		crossing_coords_spherical[e_theta] = M_PI_2;

	/*
	
	Interpolate the covariant momenta at those coorinates
	
	*/

	double momentum_param[3]{};

	for (int index = e_r; index <= e_theta; index++) {

		momentum_param[index]   = (crossing_coords_spherical[index] - Old_State_Vector[index]) / (State_Vector[index] - Old_State_Vector[index]);
		crossing_momenta[index] = momentum_param[index] * State_Vector[5 - index] + (1 - momentum_param[index]) * Old_State_Vector[5 - index];
	}

	double r_out = NT_Model.get_r_out();
	double r_in  = NT_Model.get_r_in();

	return r_crossing_squared < r_out * r_out && r_crossing_squared > r_in * r_in;

}

double dot_product(double vector_1[3], double vector_2[3]) {

	/************************************************************************************
	|                                                                                   |
	|   @ Description: Computes a simple dot product between two Euclidian 3D vectors   |
	|                                                                                   |
	|   @ Inputs:                                                                       |
	|     * vector_1: Pointer to the first vector										|
	|     * vector_2: Pointer to the second												|
	|																					|
	|   @ Ouput: The dot product between vector_1 and vector_2		                    |
	|                                                                                   |
	************************************************************************************/

	return vector_1[0] * vector_2[0] + vector_1[1] * vector_2[1] + vector_1[2] * vector_2[2];

}
