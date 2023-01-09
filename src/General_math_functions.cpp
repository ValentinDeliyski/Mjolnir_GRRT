#define _USE_MATH_DEFINES

#include "Disk_Models.h"

#include <cmath>

#include "Enumerations.h"


extern Novikov_Thorne_Model NT_Model;

double vector_norm(double Vector[], int Vector_size) {

    double norm{};

    for (int index = 0; index <= Vector_size - 1; index++) {

        norm += Vector[index] * Vector[index];

    }

    norm = sqrt(norm);

    return norm;

};

int mat_vec_multiply_4D(double Matrix[4][4], double Vector[4], double result[4]) {

    for (int row = 0; row <= 3; row += 1) {

        for (int column = 0; column <= 3; column += 1) {

            result[row] += Matrix[row][column] * Vector[column];

        }
    }

    return OK;

};

int Rorate_to_obs_plane(double theta_obs, double phi_obs, double Image_point[3], double rotated_Image_point[3]) {

    const double R_theta[3][3] =
    {
      {1,             0,                       0           },
      {0,  cos(M_PI_2 - theta_obs), sin(M_PI_2 - theta_obs)},
      {0, -sin(M_PI_2 - theta_obs), cos(M_PI_2 - theta_obs)}
    };

    const double R_phi[3][3] =
    {
      {cos(M_PI_2 - phi_obs), -sin(M_PI_2 - phi_obs), 0},
      {sin(M_PI_2 - phi_obs),  cos(M_PI_2 - phi_obs), 0},
      {           0,                      0,          1}
    };

    double Rotation_matrix[3][3]{};

    for (int row = 0; row <= 3 - 1; row += 1) {

        for (int column = 0; column <= 3 - 1; column += 1) {

            for (int k = 0; k <= 3 - 1; k += 1) {

                Rotation_matrix[row][column] += R_theta[row][k] * R_phi[k][column];

            }

        }
    }

    for (int k = x; k <= z; k += 1) {

        rotated_Image_point[k] = 0;

    }

    for (int vector_index = 0; vector_index <= 3 - 1; vector_index += 1) {

        for (int k = 0; k <= 3 - 1; k += 1) {

            rotated_Image_point[vector_index] += Rotation_matrix[vector_index][k] * Image_point[k];

        }
    }

    return OK;
}

double my_max(double vector[]) {

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

bool crossed_equatior(double State_vector[], double Old_State_Vector[]) {

    return cos(State_vector[e_theta]) * cos(Old_State_Vector[e_theta]) < 0;

}

bool interpolate_crossing(double State_Vector[], double Old_State_Vector[], double Crossing_coords[]) {

    if (cos(State_Vector[e_theta]) * cos(Old_State_Vector[e_theta]) > 0)
    {

        return 0;

    }

    double x = State_Vector[e_r] * sin(State_Vector[e_theta]) * cos(State_Vector[e_phi]);
    double y = State_Vector[e_r] * sin(State_Vector[e_theta]) * sin(State_Vector[e_phi]);
    double z = State_Vector[e_r] * cos(State_Vector[e_theta]);

    double x_old = Old_State_Vector[e_r] * sin(Old_State_Vector[e_theta]) * cos(Old_State_Vector[e_phi]);
    double y_old = Old_State_Vector[e_r] * sin(Old_State_Vector[e_theta]) * sin(Old_State_Vector[e_phi]);
    double z_old = Old_State_Vector[e_r] * cos(Old_State_Vector[e_theta]);

    double gradient[3] = { x - x_old, y - y_old, z - z_old };
    double const_term[3] = { x_old, y_old, z_old };

    double crossing_param = -const_term[2] / gradient[2];

    for (int index = 0; index < 3; index++) {

        Crossing_coords[index] = gradient[index] * crossing_param + const_term[index];

    }

    double r_crossing_squared = Crossing_coords[0] * Crossing_coords[0] + Crossing_coords[1] * Crossing_coords[1];
    double r_out = NT_Model.get_r_out();
    double r_in = NT_Model.get_r_in();

    return r_crossing_squared < r_out* r_out && r_crossing_squared > r_in * r_in;

}


double dot_product(double vector_1[3], double vector_2[3]) {

    return vector_1[0] * vector_2[0] + vector_1[1] * vector_2[1] + vector_1[2] * vector_2[2];

}
