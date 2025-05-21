from numpy import array, concatenate, tile, zeros, linspace, meshgrid, einsum, reshape, random, concatenate, tanh, outer
from numpy.linalg import inv

class Surface_Cubic_B_spline():
    
    def __init__(self, x_grid: array, y_grid: array, z_grid: array, X_patch_number: int = 5, Y_patch_number: int = 5) -> None:

        """ === The source for this script is https://hal.science/hal-03017566/document === """

        """ ==== Specify the number of patches in each coordinate direction ==== """
        self.X_patch_number = X_patch_number
        self.Y_patch_number = Y_patch_number

        self.Control_point_number = (Y_patch_number + 2) * (X_patch_number + 2)

        """ ==== Setup the fit knot positions - these are the (x, y, z) coordinate pairs of the suraface ==== """
        self.x_grid_points = x_grid.flatten()
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

        """ These two loops specify the Knot coordinate <-> Control point mappting part of the matrix. The loop runs row by row. """
        for x_idx in range(X_patch_number):
            
            for y_idx in range(Y_patch_number):
                
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 0) * (Y_patch_number + 2) + 0] = 1
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 0) * (Y_patch_number + 2) + 1] = 4
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 0) * (Y_patch_number + 2) + 2] = 1
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 1) * (Y_patch_number + 2) + 0] = 4
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 1) * (Y_patch_number + 2) + 1] = 16
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 1) * (Y_patch_number + 2) + 2] = 4
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 2) * (Y_patch_number + 2) + 0] = 1
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 2) * (Y_patch_number + 2) + 1] = 4
                Mapping_matrix[y_idx + x_idx * Y_patch_number, y_idx + (x_idx + 2) * (Y_patch_number + 2) + 2] = 1

        self.specify_free_end_conditions(Mapping_matrix = Mapping_matrix)
       
        """ === Compute the control points by inverting the mapping matrix === """
        Inverse_mapping_matrix = inv(Mapping_matrix)

        self.Control_vector_X = 36 * Inverse_mapping_matrix.dot(Knot_vector_X)
        self.Control_vector_Y = 36 * Inverse_mapping_matrix.dot(Knot_vector_Y)
        self.Control_vector_Z = 36 * Inverse_mapping_matrix.dot(Knot_vector_Z)
        
    def specify_not_a_knot_conditions(self, Mapping_matrix: array) -> None:
        
        """ === Specify the boundary conditions on the [-X, -Y] corner === """    
        Mapping_matrix[self.X_patch_number * self.Y_patch_number, 0] = 1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number, 1] = -1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number, self.X_patch_number + 2] = -1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number, self.X_patch_number + 2 + 1] = 1
     
        """ === Specify the boundary conditions on the [-X, +Y] corner === """ 
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 1, self.Y_patch_number] = -1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 1, self.Y_patch_number + 1] = 1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 1, self.Y_patch_number + self.Y_patch_number + 2] = 1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 1, self.Y_patch_number + self.Y_patch_number + 2 + 1] = -1
        
        """ === Specify the boundary conditions on the [+X, -Y] corner === """ 
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 2, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 1 - self.Y_patch_number - 3 - self.X_patch_number] = -1
        Mapping_matrix[self.X_patch_number * self.Y_patch_number + 2, (self.X_patch_number + 2) * (self.Y_patch_number + 2) - 1 - self.Y_patch_number - 2 - self.X_patch_number] = 1
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
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 1 + self.Y_patch_number + 2 - 5] = 16
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 2 + self.Y_patch_number + 2 - 5] = 24
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 3 + self.Y_patch_number + 2 - 5] = 16
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 1 * (self.Y_patch_number + 2) + 4 + self.Y_patch_number + 2 - 5] = -4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 0 + self.Y_patch_number + 2 - 5] = -1
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 1 + self.Y_patch_number + 2 - 5] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 2 + self.Y_patch_number + 2 - 5] = -6
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 3 + self.Y_patch_number + 2 - 5] =  4
            Mapping_matrix[self.X_patch_number * self.Y_patch_number + 4 + self.X_patch_number + 2 * self.Y_patch_number + idx, idx * (self.Y_patch_number + 2) + 2 * (self.Y_patch_number + 2) + 4 + self.Y_patch_number + 2 - 5] = -1 
        
    def specify_free_end_conditions(self, Mapping_matrix: array) -> None:
        
        """ === Specify the boundary condition ay -Y === """
        for idx in range(self.X_patch_number):
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + idx,(idx + 1) * (self.Y_patch_number + 2) + 0] =  1
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + idx,(idx + 1) * (self.Y_patch_number + 2) + 1] = -2
            Mapping_matrix[self.Y_patch_number * self.X_patch_number + idx,(idx + 1) * (self.Y_patch_number + 2) + 2] =  1
            
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

    def get_control_point_matrix(self, Control_vector: array, V_idx: int, U_idx: int) -> array:
    
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

    def evaluate_spline(self, Patch_discretization: int = 15) -> tuple[array, array, array]:
        
        X_surface = array([[]])
        Y_surface = array([[]])
        Z_surface = array([[]])
        
        U, V = meshgrid(linspace(0, 1, num = Patch_discretization), linspace(0, 1, num = Patch_discretization))

        for V_idx in range(0, self.Y_patch_number - 1):
            
            Partial_X_grid = []
            Partial_Y_grid = []
            Partial_Z_grid = []
            
            Patch_X_coords = []
            Patch_Y_coords = []
            Patch_Z_coords = []
            
            for U_idx in range (0, self.X_patch_number - 1):
                
                """ === Compute the basis polynomials vectors === """
                
                V1 = (1 - V)**3
                V2 =  3 * V**3 - 6 * V**2 + 4
                V3 = -3 * V**3 + 3 * V**2 + 3 * V + 1
                V4 =  V**3
                
                Basis_V_vector = array([V1, V2, V3, V4])
                
                U1 = (1 - U)**3    
                U2 = 3 * U**3 - 6 * U**2 + 4
                U3 = -3 * U**3 + 3 * U**2 + 3 * U + 1
                U4 = U**3
                
                Basis_U_vector = array([U1, U2, U3, U4])
                
                """ Evaluate the actual spline -> this uses the Knot vector and basais polynomials to compte the (x, y, z) points of the parametric surface """

                Control_point_matrix = self.get_control_point_matrix(Control_vector = self.Control_vector_X, U_idx = U_idx, V_idx = V_idx)
                Patch_X_coords.append(sum([x * y for x, y in zip(Basis_V_vector, einsum("ij,jlk->ilk", Control_point_matrix, Basis_U_vector))]) / 36)
                    
                Control_point_matrix = self.get_control_point_matrix(Control_vector = self.Control_vector_Y, U_idx = U_idx, V_idx = V_idx)
                Patch_Y_coords.append(sum([x * y for x, y in zip(Basis_V_vector, einsum("ij,jlk->ilk", Control_point_matrix, Basis_U_vector))]) / 36)

                Control_point_matrix = self.get_control_point_matrix(Control_vector = self.Control_vector_Z, U_idx = U_idx, V_idx = V_idx)
                Patch_Z_coords.append(sum([x * y for x, y in zip(Basis_V_vector, einsum("ij,jlk->ilk", Control_point_matrix, Basis_U_vector))]) / 36)
 
            for Patch_X_coords, Patch_Y_coords, Patch_Z_coords in zip(Patch_X_coords, Patch_Y_coords, Patch_Z_coords):
                
                try:
                    Partial_X_grid = concatenate((Partial_X_grid, Patch_X_coords), axis = 1)
                    Partial_Y_grid = concatenate((Partial_Y_grid, Patch_Y_coords), axis = 1)
                    Partial_Z_grid = concatenate((Partial_Z_grid, Patch_Z_coords), axis = 1)
                    
                except:
                    Partial_X_grid = Patch_X_coords
                    Partial_Y_grid = Patch_Y_coords
                    Partial_Z_grid = Patch_Z_coords
                    
            try:
                X_surface = concatenate((X_surface, Partial_X_grid), axis = 0)
                Y_surface = concatenate((Y_surface, Partial_Y_grid), axis = 0)
                Z_surface = concatenate((Z_surface, Partial_Z_grid), axis = 0)
                
            except:
                X_surface = Partial_X_grid
                Y_surface = Partial_Y_grid
                Z_surface = Partial_Z_grid

        return X_surface, Y_surface, Z_surface

if __name__ == "__main__":  
    
    Y_patch_number = 25
    X_patch_number = 25
    """ ==== Setup the fit knot positions - these are the (x, y, z) coordinate pairs of the suraface ==== """
    x_span = linspace(0, 10, X_patch_number)
    y_span = linspace(0, 15, Y_patch_number)

    # ==== This makes a grid with the x values for each point in the surface domain, then reshapes it into a 1D array
    x_grid = reshape(tile(x_span, (Y_patch_number,1)).T, (Y_patch_number * X_patch_number))

    # ==== This makes a grid with the y values for each point in the surface domain, then reshapes it into a 1D array
    y_grid = tile(y_span, X_patch_number)
    
    import numpy as np
    # ==== Placeholder for actual z values -> this will be the actual metric functions
    z = np.sqrt((outer(x_span, y_span) / 1.254))
    
    Spline_class_instance = Surface_Cubic_B_spline(x_grid = x_grid, y_grid = y_grid, z_grid = z, X_patch_number = X_patch_number, Y_patch_number= Y_patch_number)

    import matplotlib.pyplot as plt

    """ === Plot the resulting parametric surface === """
    Fig = plt.figure(figsize = (8, 8))
    Surface_subplot = Fig.add_subplot(111, projection = '3d')
    # Surface_subplot.scatter(x_grid, y_grid, z, color = 'black')
    
    x_test = 1
    y_test = 2.5

    x_bin_idx = np.digitize(x_test,x_span) - 1
    y_bin_idx = np.digitize(y_test,y_span) - 1
     
    V = (y_test - y_span[y_bin_idx]) / (y_span[y_bin_idx + 1] - y_span[y_bin_idx])
    
    V1 = (1 - V)**3
    V2 =  3 * V**3 - 6 * V**2 + 4
    V3 = -3 * V**3 + 3 * V**2 + 3 * V + 1
    V4 =  V**3
                
    Basis_V_vector = array([V1, V2, V3, V4])
    
    U = (x_test - x_span[x_bin_idx]) / (x_span[y_bin_idx + 1] - x_span[y_bin_idx])
    
    U1 = (1 - U)**3    
    U2 = 3 * U**3 - 6 * U**2 + 4
    U3 = -3 * U**3 + 3 * U**2 + 3 * U + 1
    U4 = U**3
                
    Basis_U_vector = array([U1, U2, U3, U4])
    
    Control_point_matrix = Spline_class_instance.get_control_point_matrix(Spline_class_instance.Control_vector_X, y_bin_idx, x_bin_idx)
    x_interp = np.dot(Basis_V_vector, np.dot(Control_point_matrix, Basis_U_vector)) / 36
    
    Control_point_matrix = Spline_class_instance.get_control_point_matrix(Spline_class_instance.Control_vector_Y, y_bin_idx, x_bin_idx)
    y_interp = np.dot(Basis_V_vector, np.dot(Control_point_matrix, Basis_U_vector)) / 36
    
    Control_point_matrix = Spline_class_instance.get_control_point_matrix(Spline_class_instance.Control_vector_Z, y_bin_idx, x_bin_idx)
    z_interp = np.dot(Basis_V_vector, np.dot(Control_point_matrix, Basis_U_vector)) / 36
    
    print(x_interp, y_interp, z_interp, np.sqrt(x_interp * y_interp / 1.254))
    
    X_surface, Y_surface, Z_surface = Spline_class_instance.evaluate_spline(Patch_discretization = 5)

    # Surface_subplot.plot_surface(X_surface, Y_surface, Z_surface, color = 'orange', alpha = 0.5) 
    Surface_subplot.plot_surface(reshape(x_grid, (X_patch_number , Y_patch_number)), 
                                 reshape(y_grid, (X_patch_number , Y_patch_number)), 
                                 z, color = 'orange', alpha = 0.5) 
    
    Surface_subplot.scatter(x_interp, y_interp, z_interp, color = 'black')

    Surface_subplot.set_xlabel('x')
    Surface_subplot.set_ylabel('y')
    Surface_subplot.set_zlabel('z')
    Surface_subplot.set_zlim(5,-5)
    plt.show()