from numpy import array, concatenate, tile, zeros, linspace, meshgrid, einsum, reshape, concatenate, outer, digitize, dot, sqrt, float64, roll, append, ones
from numpy.typing import NDArray
from numpy.linalg import inv
import matplotlib.pyplot as plt

class Surface_Cubic_B_spline():
    
    def __init__(self, x_grid: NDArray[float64], y_grid: NDArray[float64], z_grid: NDArray[float64], X_patch_number: int = 5, Y_patch_number: int = 5) -> None:

        """ === The source for this script is https://hal.science/hal-03017566/document === """

        """ ==== Specify the number of patches in each coordinate direction ==== """
        self.X_patch_number = X_patch_number
        self.Y_patch_number = Y_patch_number

        self.Control_point_number = (Y_patch_number + 2) * (X_patch_number + 2)

        """ ==== Setup the fit knot positions - these are the (x, y, z) coordinate pairs of the suraface ==== """
        self.x_grid_steps = abs(x_grid.T[0] - roll(x_grid.T[0], 1))[1:]
        self.x_grid_points = x_grid.flatten()
        
        self.y_grid_steps = abs(y_grid[0] - roll(y_grid[0], 1))[1:]
        self.y_grid_points = y_grid.flatten()
        
        self.z_grid_points = z_grid.flatten()

        """ Construct the knot vectors - they contain the known surface points + appended 0's for the end conditions. 
            In order to fully specify the splined surface, a total of (X_patch_number + 2) * (Y_patch_number + 2) - X_patch_number * Y_patch_number 
            additional conditions must be specified. They take on the form "Sum_ij a_ij Q_ij = 0", where Q_ij are the control points of the splined
            surface. The "= 0" part is where the appended 0's in the knot vector come from """

        Boundary_condition_number = self.Control_point_number - X_patch_number * Y_patch_number

        Knot_vector_X = concatenate((self.x_grid_points, tile([0], Boundary_condition_number)))
        Knot_vector_Y = concatenate((self.y_grid_points, tile([0], Boundary_condition_number)))
        Knot_vector_Z = concatenate((self.z_grid_points, tile([0], Boundary_condition_number)))

        """ Construct the mapping matrix between the knots, and control points -> Knot_vector = Mapping_matrix * Control_point_vector.
            It is a sparse matrix, so we initialize it to zero and fill in the few non -zero components. """
        Mapping_matrix = zeros((self.Control_point_number, self.Control_point_number))

        self.x_grid_steps = append(self.x_grid_steps, [self.x_grid_steps[-1], self.x_grid_steps[-1], self.x_grid_steps[-1]])
        self.x_grid_steps = append([self.x_grid_steps[0], self.x_grid_steps[0], self.x_grid_steps[0]], self.x_grid_steps)
        
        self.y_grid_steps = append(self.y_grid_steps, [self.y_grid_steps[-1], self.y_grid_steps[-1], self.y_grid_steps[-1]])
        self.y_grid_steps = append([self.y_grid_steps[0], self.y_grid_steps[0], self.y_grid_steps[0]], self.y_grid_steps)
        
        """ These two loops specify the Knot coordinate <-> Control point mappting part of the matrix. The loop runs row by row. """
        for x_idx in range(X_patch_number):
            
            a_coeff_X, _, _, _, _, f_coeff_X = self.get_delta_coeffs(self.x_grid_steps, x_idx + 3)
            
            for y_idx in range(Y_patch_number):
                           
                a_coeff_Y, _, _, _, _, f_coeff_Y = self.get_delta_coeffs(self.y_grid_steps, y_idx + 3)
                
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 0) * (Y_patch_number + 2) + 0] = a_coeff_Y * a_coeff_X 
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 0) * (Y_patch_number + 2) + 1] = a_coeff_Y * (1 - f_coeff_X - a_coeff_X)
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 0) * (Y_patch_number + 2) + 2] = a_coeff_Y * f_coeff_X
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 1) * (Y_patch_number + 2) + 0] = (1 - a_coeff_Y - f_coeff_Y) * a_coeff_X
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 1) * (Y_patch_number + 2) + 1] = (1 - a_coeff_Y - f_coeff_Y) * (1 - f_coeff_X - a_coeff_X)
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 1) * (Y_patch_number + 2) + 2] = (1 - a_coeff_Y - f_coeff_Y) * f_coeff_X
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 2) * (Y_patch_number + 2) + 0] = f_coeff_Y * a_coeff_X
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 2) * (Y_patch_number + 2) + 1] = f_coeff_Y * (1 - a_coeff_X - f_coeff_X)
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 2) * (Y_patch_number + 2) + 2] = f_coeff_Y * f_coeff_X

        self.specify_free_end_conditions(Mapping_matrix = Mapping_matrix)
       
        """ === Compute the control points by inverting the mapping matrix === """
        Inverse_mapping_matrix = inv(Mapping_matrix)

        self.Control_vector_X = Inverse_mapping_matrix.dot(Knot_vector_X)
        self.Control_vector_Y = Inverse_mapping_matrix.dot(Knot_vector_Y)
        self.Control_vector_Z = Inverse_mapping_matrix.dot(Knot_vector_Z)
    
    def get_delta_coeffs(self, Grid_steps: NDArray[float64], idx: int) -> tuple[float, float, float, float, float, float]:
        
        a_coeff = Grid_steps[idx]**2 / ((Grid_steps[idx - 2] + Grid_steps[idx - 1] + Grid_steps[idx]) * (Grid_steps[idx - 1] + Grid_steps[idx]))

        b_coeff = Grid_steps[idx]**2 / ((Grid_steps[idx - 1] + Grid_steps[idx] + Grid_steps[idx + 1]) * (Grid_steps[idx - 1] + Grid_steps[idx]))
   
        c_coeff = Grid_steps[idx]**2 / ((Grid_steps[idx - 1] + Grid_steps[idx] + Grid_steps[idx + 1]) * (Grid_steps[idx] + Grid_steps[idx + 1]))
        
        d_coeff = Grid_steps[idx]**2 / ((Grid_steps[idx] + Grid_steps[idx + 1] + Grid_steps[idx + 2]) * (Grid_steps[idx] + Grid_steps[idx + 1]))
        
        e_coeff = Grid_steps[idx] * Grid_steps[idx - 1] / ((Grid_steps[idx - 1] + Grid_steps[idx] + Grid_steps[idx + 1]) * (Grid_steps[idx - 1] + Grid_steps[idx]))
        
        f_coeff =  Grid_steps[idx - 1]**2 / ((Grid_steps[idx - 1] + Grid_steps[idx] + Grid_steps[idx + 1]) * (Grid_steps[idx - 1] + Grid_steps[idx]))
   
        return a_coeff, b_coeff, c_coeff, d_coeff, e_coeff, f_coeff
    
    def specify_not_a_knot_conditions(self, Mapping_matrix: NDArray[float64]) -> None:
        
        """ === Specify the boundary conditions on the [-X, -Y] corner === """    
        Mapping_matrix[self.X_patch_number * self.Y_patch_number, 0] = 1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number, 1] = -1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number, self.Y_patch_number + 2] = -1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number, self.Y_patch_number + 2 + 1] = 1
     
        """ === Specify the boundary conditions on the [-X, +Y] corner === """ 
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 1, self.Y_patch_number] = -1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 1, self.Y_patch_number + 1] = 1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 1, self.Y_patch_number + self.Y_patch_number + 2] = 1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 1, self.Y_patch_number + self.Y_patch_number + 2 + 1] = -1
        
        """ === Specify the boundary conditions on the [+X, -Y] corner === """ 
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 2, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 1 - self.Y_patch_number - 3 - self.Y_patch_number] = -1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 2, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 1 - self.Y_patch_number - 2 - self.Y_patch_number] = 1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 2, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 1 - 1 - self.Y_patch_number] = 1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 2, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 1 - self.Y_patch_number] = -1
        
        """ === Specify the boundary conditions on the [+X, +Y] corner === """ 
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 3, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 1 - self.Y_patch_number - 2 - 1] = 1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 3, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 1 - self.Y_patch_number - 2] = -1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 3, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 1 - 1] = -1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 3, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 1] = 1
        
        """ === Specify the boundary condition ay -Y === """     
        for idx in range(self.Y_patch_number):
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 0] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 1] = -4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 2] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 0 + 1 * (self.Y_patch_number+2)] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 1 + 1 * (self.Y_patch_number+2)] = 16
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 2 + 1 * (self.Y_patch_number+2)] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 0 + 2 * (self.Y_patch_number+2)] = -6
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 1 + 2 * (self.Y_patch_number+2)] = -24
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 2 + 2 * (self.Y_patch_number+2)] = -6
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 0 + 3 * (self.Y_patch_number+2)] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 1 + 3 * (self.Y_patch_number+2)] = 16
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 2 + 3 * (self.Y_patch_number+2)] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 0 + 4 * (self.Y_patch_number+2)] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 1 + 4 * (self.Y_patch_number+2)] = -4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + idx, idx + 2 + 4 * (self.Y_patch_number+2)] = -1
            
        """ === Specify the boundary condition ay +Y === """        
        for idx in range(self.Y_patch_number):
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + 0] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + 1] = -4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + 2] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + (self.Y_patch_number + 2) + 0] = 4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + (self.Y_patch_number + 2) + 1] = 16
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + (self.Y_patch_number + 2) + 2] = 4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 0] = -6
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 1] = -24
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 2] = -6
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + 3 * (self.Y_patch_number + 2) + 0] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + 3 * (self.Y_patch_number + 2) + 1] = 16
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + 3 * (self.Y_patch_number + 2) + 2] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + 4 * (self.Y_patch_number + 2) + 0] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + 4 * (self.Y_patch_number + 2) + 1] = -4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.Y_patch_number + idx, idx + (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 5 * (self.Y_patch_number + 2) + 4 * (self.Y_patch_number + 2) + 2] = -1
            
        """ === Specify the boundary condition ay -X === """        
        for idx in range(self.X_patch_number):
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 0] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2] = -6
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 3] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 4] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 0] = -4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 1] =  16
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 2] = -24
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 3] =  16
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 4] = -4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 0] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 1] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 2] = -6
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 3] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 4] = -1
            
        """ === Specify the boundary condition ay +X === """      
        for idx in range(self.X_patch_number):
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 0 + self.Y_patch_number + 2 - 5] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 + self.Y_patch_number + 2 - 5] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 + self.Y_patch_number + 2 - 5] = -6
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 3 + self.Y_patch_number + 2 - 5] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 4 + self.Y_patch_number + 2 - 5] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 0 + self.Y_patch_number + 2 - 5] = -4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 1 + self.Y_patch_number + 2 - 5] =  16
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 2 + self.Y_patch_number + 2 - 5] = -24
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 3 + self.Y_patch_number + 2 - 5] =  16
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 4 + self.Y_patch_number + 2 - 5] = -4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 0 + self.Y_patch_number + 2 - 5] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 1 + self.Y_patch_number + 2 - 5] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 2 + self.Y_patch_number + 2 - 5] = -6
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 3 + self.Y_patch_number + 2 - 5] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 4 + self.Y_patch_number + 2 - 5] = -1 
        
    def specify_free_end_conditions(self, Mapping_matrix: NDArray[float64]) -> None:
        
        """ === Specify the boundary condition ay -Y === """
        for idx in range(self.X_patch_number):
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + idx, (idx + 1) * (self.Y_patch_number + 2) + 0] =  1
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + idx, (idx + 1) * (self.Y_patch_number + 2) + 1] = -2
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + idx, (idx + 1) * (self.Y_patch_number + 2) + 2] =  1
            
        """ === Specify the boundary condition ay +Y === """
        for idx in range(self.X_patch_number):
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + self.X_patch_number + idx, (idx + 1) * (self.Y_patch_number + 2) + (self.Y_patch_number + 2 - 3)] =  1
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + self.X_patch_number + idx, (idx + 1) * (self.Y_patch_number + 2) + (self.Y_patch_number + 2 - 2)] = -2
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + self.X_patch_number + idx, (idx + 1) * (self.Y_patch_number + 2) + (self.Y_patch_number + 2 - 1)] =  1
            
        """ === Specify the boundary condition ay +X === """
        for idx in range(self.Y_patch_number):
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + idx, 1 + idx + 0 * (self.Y_patch_number + 2)] =  1
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + idx, 1 + idx + 1 * (self.Y_patch_number + 2)] = -2
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + idx, 1 + idx + 2 * (self.Y_patch_number + 2)] =  1

        """ === Specify the boundary condition ay -X === """
        for idx in range(self.Y_patch_number):
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + self.Y_patch_number + idx,-1 - (1 + idx + 0 * (self.Y_patch_number + 2))] =  1
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + self.Y_patch_number + idx,-1 - (1 + idx + 1 * (self.Y_patch_number + 2))] = -2
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + self.Y_patch_number + idx,-1 - (1 + idx + 2 * (self.Y_patch_number + 2))] =  1  
            
        """ === Specify the boundary conditions on the [-X, -Y] corner === """   
        Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + 2 * self.Y_patch_number, 0 * (self.Y_patch_number + 2) + 0] =  1
        Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + 2 * self.Y_patch_number, 1 * (self.Y_patch_number + 2) + 1] = -2
        Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + 2 * self.Y_patch_number, 2 * (self.Y_patch_number + 2) + 2] =  1

        """ === Specify the boundary conditions on the [-X, +Y] corner === """     
        Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + 2 * self.Y_patch_number + 1, 1 * (self.Y_patch_number + 2) - 1] = 1
        Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + 2 * self.Y_patch_number + 1, 1 * (self.Y_patch_number + 2) + (self.Y_patch_number + 2 - 1) - 1] = -2
        Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + 2 * self.Y_patch_number + 1, 2 * (self.Y_patch_number + 2) + (self.Y_patch_number + 2 - 2) - 1] =  1

        """ === Specify the boundary conditions on the [+X, -Y] corner === """   
        Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + 2 * self.Y_patch_number + 2, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 1 * (self.Y_patch_number + 2) + 0] =  1
        Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + 2 * self.Y_patch_number + 2, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 2 * (self.Y_patch_number + 2) + 1] = -2
        Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + 2 * self.Y_patch_number + 2, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 3 * (self.Y_patch_number + 2) + 2] =  1

        """ === Specify the boundary conditions on the [+X, +Y] corner === """     
        Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + 2 * self.Y_patch_number + 3, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 2 * (self.Y_patch_number + 2) - 3] =  1
        Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + 2 * self.Y_patch_number + 3, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 1 * (self.Y_patch_number + 2) - 2] = -2
        Mapping_matrix[self.Y_patch_number * self.X_patch_number + 2 * self.X_patch_number + 2 * self.Y_patch_number + 3, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 0 * (self.Y_patch_number + 2) - 1] =  1

    def get_control_point_matrix(self, Control_vector: NDArray[float64], V_idx: int, U_idx: int) -> NDArray[float64]:
    
        Q_i0_j0 = Control_vector[V_idx + 0 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        Q_i1_j0 = Control_vector[V_idx + 1 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        Q_i2_j0 = Control_vector[V_idx + 2 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        Q_i3_j0 = Control_vector[V_idx + 3 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        
        Q_i0_j1 = Control_vector[V_idx + 1 + 0 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        Q_i1_j1 = Control_vector[V_idx + 1 + 1 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        Q_i2_j1 = Control_vector[V_idx + 1 + 2 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        Q_i3_j1 = Control_vector[V_idx + 1 + 3 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        
        Q_i0_j2 = Control_vector[V_idx + 2 + 0 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        Q_i1_j2 = Control_vector[V_idx + 2 + 1 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        Q_i2_j2 = Control_vector[V_idx + 2 + 2 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        Q_i3_j2 = Control_vector[V_idx + 2 + 3 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        
        Q_i0_j3 = Control_vector[V_idx + 3 + 0 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        Q_i1_j3 = Control_vector[V_idx + 3 + 1 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        Q_i2_j3 = Control_vector[V_idx + 3 + 2 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        Q_i3_j3 = Control_vector[V_idx + 3 + 3 * (self.Y_patch_number + 2) + U_idx * (self.Y_patch_number + 2)]
        
        Control_point_matrix = array([[Q_i0_j0, Q_i1_j0, Q_i2_j0, Q_i3_j0],
                                      [Q_i0_j1, Q_i1_j1, Q_i2_j1, Q_i3_j1],
                                      [Q_i0_j2, Q_i1_j2, Q_i2_j2, Q_i3_j2],
                                      [Q_i0_j3, Q_i1_j3, Q_i2_j3, Q_i3_j3]])
        
        return Control_point_matrix

    def evaluate_spline(self, Patch_discretization: int = 15) -> tuple[NDArray[float64], NDArray[float64], NDArray[float64]]:
        
        X_surface: NDArray[float64] = array([[]])
        Y_surface: NDArray[float64] = array([[]])
        Z_surface: NDArray[float64] = array([[]])
        
        U, V = meshgrid(linspace(0, 1, num = Patch_discretization), linspace(0, 1, num = Patch_discretization))

        for V_idx in range(0, self.Y_patch_number - 1):
            
            Partial_X_grid: NDArray[float64] = array([])
            Partial_Y_grid: NDArray[float64] = array([])
            Partial_Z_grid: NDArray[float64] = array([])
            
            Patch_X_coords: list[float] = []
            Patch_Y_coords: list[float] = []
            Patch_Z_coords: list[float] = []
            
            a_coeff_V, b_coeff_V, c_coeff_V, d_coeff_V, e_coeff_V, f_coeff_V = self.get_delta_coeffs(self.y_grid_steps, V_idx + 3)
            
            for U_idx in range(0, self.X_patch_number - 1):
                
                
                a_coeff_U, b_coeff_U, c_coeff_U, d_coeff_U, e_coeff_U, f_coeff_U = self.get_delta_coeffs(self.x_grid_steps, U_idx + 3)
                
                """ === Compute the basis polynomials vectors === """
                
                V1 = -a_coeff_V * V**3 + 3 * a_coeff_V * V**2 - 3 * a_coeff_V * V + a_coeff_V 
                V2 = (a_coeff_V + b_coeff_V + c_coeff_V) * V**3 + (-3 * a_coeff_V - 3 * b_coeff_V) * V**2 + (3 * a_coeff_V - 3 * e_coeff_V) * V + 1 - a_coeff_V - f_coeff_V
                V3 = (-b_coeff_V - c_coeff_V - d_coeff_V) * V**3 + 3 * b_coeff_V * V**2 + 3 * e_coeff_V * V + f_coeff_V
                V4 = V**3 * d_coeff_V
                
                Basis_V_vector = array([V1, V2, V3, V4])
                
                U1 = -a_coeff_U * U**3 + 3 * a_coeff_U * U**2 - 3 * a_coeff_U * U + a_coeff_U 
                U2 = (a_coeff_U + b_coeff_U + c_coeff_U) * U**3 + (-3 * a_coeff_U - 3 * b_coeff_U) * U**2 + (3 * a_coeff_U - 3 * e_coeff_U) * U + 1 - a_coeff_U - f_coeff_U
                U3 = (-b_coeff_U - c_coeff_U - d_coeff_U) * U**3 + 3 * b_coeff_U * U**2 + 3 * e_coeff_U * U + f_coeff_U
                U4 = U**3 * d_coeff_U
                
                Basis_U_vector = array([U1, U2, U3, U4])
                
                """ Evaluate the actual spline -> this uses the Knot vector and basais polynomials to compte the (x, y, z) points of the parametric surface """

                Control_point_matrix = self.get_control_point_matrix(Control_vector = self.Control_vector_X, U_idx = U_idx, V_idx = V_idx)
                Patch_X_coords.append(sum([x * y for x, y in zip(Basis_V_vector, einsum("ij,jlk->ilk", Control_point_matrix, Basis_U_vector))]))
                    
                Control_point_matrix = self.get_control_point_matrix(Control_vector = self.Control_vector_Y, U_idx = U_idx, V_idx = V_idx)
                Patch_Y_coords.append(sum([x * y for x, y in zip(Basis_V_vector, einsum("ij,jlk->ilk", Control_point_matrix, Basis_U_vector))]))

                Control_point_matrix = self.get_control_point_matrix(Control_vector = self.Control_vector_Z, U_idx = U_idx, V_idx = V_idx)
                Patch_Z_coords.append(sum([x * y for x, y in zip(Basis_V_vector, einsum("ij,jlk->ilk", Control_point_matrix, Basis_U_vector))]))
 
            for Patch_X_coords_element, Patch_Y_coords_element, Patch_Z_coords_element in zip(Patch_X_coords, Patch_Y_coords, Patch_Z_coords):
                
                try:
                    Partial_X_grid = concatenate((Partial_X_grid, Patch_X_coords_element), axis = 1)
                    Partial_Y_grid = concatenate((Partial_Y_grid, Patch_Y_coords_element), axis = 1)
                    Partial_Z_grid = concatenate((Partial_Z_grid, Patch_Z_coords_element), axis = 1)
                    
                except:
                    Partial_X_grid = array(Patch_X_coords_element)
                    Partial_Y_grid = array(Patch_Y_coords_element)
                    Partial_Z_grid = array(Patch_Z_coords_element)
                    
            try:
                X_surface = concatenate((X_surface, Partial_X_grid), axis = 0)
                Y_surface = concatenate((Y_surface, Partial_Y_grid), axis = 0)
                Z_surface = concatenate((Z_surface, Partial_Z_grid), axis = 0)
                
            except:
                X_surface = array(Partial_X_grid)
                Y_surface = array(Partial_Y_grid)
                Z_surface = array(Partial_Z_grid)

        return X_surface, Y_surface, Z_surface