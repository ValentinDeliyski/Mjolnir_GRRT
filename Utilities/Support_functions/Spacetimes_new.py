from numpy import sqrt, sin, diag, pi, cos, arccos, infty, sign, exp, cosh, arccosh, roots, array, append, flip, digitize, dot
from numpy.typing import NDArray

from Support_functions.Surface_Cubic_B_spline import Surface_Cubic_B_spline

from enum import Enum

class Coords(Enum):

    e_t = 0
    e_r = 1
    e_theta = 2
    e_phi = 3

    e_x = 0
    e_y = 1
    e_z = 2

class Spacetime:
    
    class Anzatz_enums(Enum):
        
        Anzatz_0 = 0
        Anzatz_1 = 1
        
    class Derivative_selector(Enum):
        
        No_deriv = 0
        First_radial = 1
    
    class Spacetime_enums(Enum):
        
        Kerr = 0
        Wormhole = 1
        Gauss_Bonnet = 2
        Janis_Newman_Winicour = 3
        Regular_Black_Hole = 4
        Numerical = 5

    def __init__(self) -> None:
        
        self.M: float
        self.a: float
        self.Param: float
        self.R_throat: float
        self.Use_global_coords: float

        print("Initializing a base spacetime class - something broke!")
        exit(1)

    def get_metric(self, r: float, theta: float) -> NDArray:

        print("Using base spacetime class - something broke in the get_metric function!")
        exit(1)

    def get_dr_metric(self, r: float, theta: float) -> NDArray:

        print("Using base spacetime class - something broke in the get_dr_metric function!")
        exit(1)

    def get_photon_orbits(self) -> tuple[float, float]:

        print("Using base spacetime class - something broke in the get_photon_orbits function!")
        exit(1)

    def get_ISCO(self) -> tuple[float, float]:
        
        print("Using base spacetime class - something broke in the get_ISCO function!")
        exit(1)
    
    def Identify(self) -> int:
        
        print("Using base spacetime class - something broke in the Identify function!")
        exit(1)
        

class Numerical(Spacetime):
    
    def __init__(self, Anzatz: Spacetime.Anzatz_enums, Horizon_radius: float, Radial_grid_size: int, Theta_grid_size: int, a_ADM: float, M_ADM: float, XML_spline_path: str | None = None, Raw_data_path: str | None = None) -> None:
        
        self.r_H = Horizon_radius
        self.a = a_ADM
        self.M = M_ADM
        self.Anzatz = Anzatz.value
        self.Raw_data_path = Raw_data_path
        self.XML_spline_path = XML_spline_path
        
        self.Radial_grid_size = Radial_grid_size
        self.Theta_grid_size = Theta_grid_size
        
        if self.XML_spline_path == None:
            
            self.__parse_raw_metric__()
            
            self.F_0_spline_instance = Surface_Cubic_B_spline(x_grid = self.theta_coord, y_grid = self.x_coord, z_grid = self.F_0, X_patch_number = 2 * self.Theta_grid_size - 1, Y_patch_number = self.Radial_grid_size)
            self.F_1_spline_instance = Surface_Cubic_B_spline(x_grid = self.theta_coord, y_grid = self.x_coord, z_grid = self.F_1, X_patch_number = 2 * self.Theta_grid_size - 1, Y_patch_number = self.Radial_grid_size)
            self.F_2_spline_instance = Surface_Cubic_B_spline(x_grid = self.theta_coord, y_grid = self.x_coord, z_grid = self.F_2, X_patch_number = 2 * self.Theta_grid_size - 1, Y_patch_number = self.Radial_grid_size)
            self.W_spline_instance = Surface_Cubic_B_spline(x_grid = self.theta_coord, y_grid = self.x_coord, z_grid = self.W, X_patch_number = 2 * self.Theta_grid_size - 1, Y_patch_number = self.Radial_grid_size)

    def __parse_raw_metric__(self) -> None:
        
        """ ====================== Initialize the arrays that hold the metric functions ====================== """
        x_coord = []
        theta_coord = []
        F_0 = []
        F_1 = []
        F_2 = []
        W = []

        if self.Raw_data_path == None:
            
            print("ERROR: No file path to the raw metric data was provided!")
            exit(1)

        with open(self.Raw_data_path, "r") as file:
            
            for line in file:
                
                Line_contents = line.strip().split(" ")
                Line_contents = [x for x in Line_contents if "" != x]
                
                """ Different theta values on the grid are seperated by a "\n" character - parsing it results in an empty list. """
                if len(Line_contents) == 0:
                    continue

                x_coord.append(float(Line_contents[0]))
                theta_coord.append(float(Line_contents[1]))
                
                match self.Anzatz:
                    
                    case self.Anzatz_enums.Anzatz_0.value:
                        
                        F_0.append(float(Line_contents[2]))
                        F_2.append(float(Line_contents[3]))
                        F_2.append(float(Line_contents[4]))
                        W.append(-float(Line_contents[5]))
                    
                    case self.Anzatz_enums.Anzatz_1.value:
                
                        F_1.append(float(Line_contents[2]))
                        F_2.append(float(Line_contents[3]))
                        F_0.append(float(Line_contents[4]))
                        W.append(float(Line_contents[5]))
                    
                    case _:
                        
                        print("ERROR: Invalid metric anzatz!")
                        exit(1)
                        
                        
        """ The metric is calculated only for theta values in the range [0, pi / 2]. We use the reflection symmetry of the problem to get the rest of the grid. """
        x_coord = array(x_coord).reshape(self.Theta_grid_size, self.Radial_grid_size)
        self.x_coord = append(x_coord, x_coord[1:], axis = 0)
        
        theta_coord = array(theta_coord).reshape(self.Theta_grid_size, self.Radial_grid_size)
        self.theta_coord = append(theta_coord, flip(pi - theta_coord)[1:], axis = 0)
         
        F_0 = array(F_0).reshape(self.Theta_grid_size, self.Radial_grid_size)  
        self.F_0 = append(F_0, flip(F_0, axis = 0)[1:], axis = 0)  
        
        F_1 = array(F_1).reshape(self.Theta_grid_size, self.Radial_grid_size)  
        self.F_1 = append(F_1, flip(F_1, axis = 0)[1:], axis = 0)  
        
        F_2 = array(F_2).reshape(self.Theta_grid_size, self.Radial_grid_size) 
        self.F_2 = append(F_2, flip(F_2, axis = 0)[1:], axis = 0)  
        
        W = array(W).reshape(self.Theta_grid_size, self.Radial_grid_size)
        self.W = append(W, flip(W, axis = 0)[1:], axis = 0)
        
        """ Convert the compactified coordinate x to the unbounded radial coordinate """
        self.r_coord = sqrt((self.x_coord / (1 - self.x_coord))**2 + self.r_H**2)
        
    def __compactify_radial_coordinate__(self, r: float) -> float:
        
        """ Comute the shifted coordinate that the paper (http://gravitation.web.ua.pt/node/416) uses in the numerical implementation (they label this little "x"). """
        x_uncompactified = sqrt(r**2 - self.r_H**2)
        
        return x_uncompactified / (1 + x_uncompactified)
        
    def __evaluate_spline__(self, r: float, theta: float, Derivative_selector: Spacetime.Derivative_selector) -> tuple[float, float, float, float]:
        
        X_coord_span = self.x_coord[0]
        X_coord = self.__compactify_radial_coordinate__(r)   
        
        Theta_coord_span = self.theta_coord.T[0]
        
        theta_bin_idx: int = int(digitize(theta, Theta_coord_span, right = True))
        X_bin_idx: int = int(digitize(X_coord, X_coord_span, right = True))
             
        V = (X_coord - X_coord_span[X_bin_idx - 1]) / (X_coord_span[X_bin_idx] - X_coord_span[X_bin_idx - 1])
        
        match Derivative_selector.value:
            
            case Spacetime.Derivative_selector.No_deriv.value:
        
                V1 = (1 - V)**3
                V2 =  3 * V**3 - 6 * V**2 + 4
                V3 = -3 * V**3 + 3 * V**2 + 3 * V + 1
                V4 =  V**3
           
            case Spacetime.Derivative_selector.First_radial.value:
                
                V1 = -3 * (1 - V)**2
                V2 =  9 * V**2 - 12 * V
                V3 = -9 * V**2 + 6 * V + 3
                V4 =  3 * V**2
                
            case _:
                
                print("ERROR: Invalid derivaive enum!")
                exit(1)
                    
        Basis_V_vector = array([V1, V2, V3, V4])
        
        U = (theta - Theta_coord_span[theta_bin_idx - 1]) / (Theta_coord_span[theta_bin_idx] - Theta_coord_span[theta_bin_idx - 1])
        
        U1 = (1 - U)**3    
        U2 = 3 * U**3 - 6 * U**2 + 4
        U3 = -3 * U**3 + 3 * U**2 + 3 * U + 1
        U4 = U**3
                    
        Basis_U_vector = array([U1, U2, U3, U4])
        
        F_0_Control_point_matrix = self.F_0_spline_instance.get_control_point_matrix(self.F_0_spline_instance.Control_vector_Z, X_bin_idx - 1, theta_bin_idx - 1)
        F_0_value = dot(Basis_V_vector, dot(F_0_Control_point_matrix, Basis_U_vector)) / 36
        
        F_1_Control_point_matrix = self.F_1_spline_instance.get_control_point_matrix(self.F_1_spline_instance.Control_vector_Z, X_bin_idx - 1, theta_bin_idx - 1)
        F_1_value = dot(Basis_V_vector, dot(F_1_Control_point_matrix, Basis_U_vector)) / 36
        
        F_2_Control_point_matrix = self.F_2_spline_instance.get_control_point_matrix(self.F_2_spline_instance.Control_vector_Z, X_bin_idx - 1, theta_bin_idx - 1)
        F_2_value = dot(Basis_V_vector, dot(F_2_Control_point_matrix, Basis_U_vector)) / 36
        
        W_Control_point_matrix = self.W_spline_instance.get_control_point_matrix(self.W_spline_instance.Control_vector_Z, X_bin_idx - 1, theta_bin_idx - 1)
        W_value = dot(Basis_V_vector, dot(W_Control_point_matrix, Basis_U_vector)) / 36
        
        match Derivative_selector.value:
            
            case Spacetime.Derivative_selector.No_deriv.value:
                return F_0_value, F_1_value, F_2_value, W_value
           
            case Spacetime.Derivative_selector.First_radial.value: 
                
                Deriv_correction_factor = (1 - X_coord)**3 / (X_coord_span[X_bin_idx] - X_coord_span[X_bin_idx - 1]) * r / X_coord
                 
                return F_0_value * Deriv_correction_factor, F_1_value * Deriv_correction_factor, F_2_value * Deriv_correction_factor, W_value * Deriv_correction_factor
                
    def get_metric(self, r: float, theta: float) -> NDArray:

        F_0, F_1, F_2, W = self.__evaluate_spline__(r = r, theta = theta, Derivative_selector = self.Derivative_selector.No_deriv)

        N = (1 - self.r_H / r)

        match self.Anzatz:
            
            case self.Anzatz_enums.Anzatz_0.value:
                
                g_tt = -exp(2 * F_0) * N + exp(2 * F_2) * W**2 * sin(theta)**2
                g_tphi = -exp(2 * F_2) * r * W * sin(theta)**2
                
            case self.Anzatz_enums.Anzatz_1.value:
                
                g_tt = -exp(2 * F_0) * N + exp(2 * F_2) * r**2 * W**2 * sin(theta)**2
                g_tphi = -exp(2 * F_2) * r**2 * W * sin(theta)**2
                
            case _:
                
                print("ERROR: Invalid metric anzatz")
                exit(1)
                
        g_rr = exp(2 * F_1) / N
        g_thth = exp(2 * F_1) * r**2
        g_phiphi = exp(2 * F_2) * r**2 * sin(theta)**2
        
        Metric = diag([g_tt, g_rr, g_thth, g_phiphi])
        Metric[Coords.e_t.value][Coords.e_phi.value] = g_tphi
        Metric[Coords.e_phi.value][Coords.e_t.value] = g_tphi
        
        return Metric
    
    def Identify(self) -> int:
        
        return self.Spacetime_enums.Numerical.value
    
    def get_dr_metric(self, r: float, theta: float) -> NDArray:

        F_0, F_1, F_2, W = self.__evaluate_spline__(r = r, theta = theta, Derivative_selector = self.Derivative_selector.No_deriv)
        dr_F_0, dr_F_1, dr_F_2, dr_W = self.__evaluate_spline__(r = r, theta = theta, Derivative_selector = self.Derivative_selector.First_radial)

        N = (1 - self.r_H / r)
        dr_N = self.r_H / r**2

        match self.Anzatz:
            
            case self.Anzatz_enums.Anzatz_0.value:
                
                dr_g_tt =  -exp(2 * F_0) * (2 * N * dr_F_0 + dr_N) + 2 * exp(2 * F_2) * sin(theta)**2 * W * (W * dr_F_2 + dr_W)
                dr_g_tphi = -exp(2 * F_2) * r * sin(theta)**2 * (2 * r * W * dr_F_2 + 2 * W + r * dr_W)
                
            case self.Anzatz_enums.Anzatz_1.value:
                
                dr_g_tt = -exp(2 * F_0) * (2 * N * dr_F_0 + dr_N) + 2 * exp(2 * F_2) * r * sin(theta)**2 * W * (r * W * dr_F_2 + W + r * dr_W)
                dr_g_tphi = -exp(2 * F_2) * r * sin(theta)**2 * (2 * r * W * dr_F_2 + 2 * W + r * dr_W)
                
            case _:
                
                print("ERROR: Invalid metric anzatz")
                exit(1)
                
        dr_g_rr = exp(2 * F_1) / N * (2 * dr_F_1 - dr_N / N)
        dr_g_thth = 2 * r * exp(2 * F_1) * (r * dr_F_1 + 1)
        dr_g_phiphi = 2 * r * exp(2 * F_2) * (r * dr_F_2 + 1) * sin(theta)**2

        dr_Metric = diag([dr_g_tt, dr_g_rr, dr_g_thth, dr_g_phiphi])
        dr_Metric[Coords.e_t.value][Coords.e_phi.value] = dr_g_tphi
        dr_Metric[Coords.e_phi.value][Coords.e_t.value] = dr_g_tphi
        
        return dr_Metric
    
class Kerr(Spacetime):

    def __init__(self, mass: float = 1, spin_param: float = 0) -> None:
        self.M = mass
        self.a = spin_param

    def get_metric(self, r: float, theta: float) -> NDArray:

        Sigma = r**2 + self.a**2 * cos(theta)**2
        Delta = r**2 - 2 * self.M * r + self.a**2

        """ Handles division by zero that happens at the singularity. """
        g_tt = -(1 - 2 * self.M * r / Sigma) if Sigma else infty

        """ Handles division by zero that happens at the horizons. """
        g_rr = Sigma / Delta if Delta else infty

        g_thth = Sigma

        """ Handles division by zero that happens at the singularity. """
        g_phiphi = (r**2 + self.a**2 + 2 * self.M * r * self.a**2 / Sigma * sin(theta)**2) * sin(theta)**2 if Sigma else infty

        """ Handles division by zero that happens at the singularity. """
        g_tphi = -2 * self.M * r * self.a * sin(theta)**2 / Sigma if Sigma else -sign(self.a) * infty

        Metric = diag([g_tt, g_rr, g_thth, g_phiphi])
        Metric[Coords.e_t.value][Coords.e_phi.value] = g_tphi
        Metric[Coords.e_phi.value][Coords.e_t.value] = g_tphi
        
        return Metric
    
    def Identify(self) -> int:
        
        return self.Spacetime_enums.Kerr.value
    
    def get_dr_metric(self, r: float, theta: float) -> NDArray:

        Sigma = r**2 + self.a**2 * cos(theta)**2
        Delta = r**2 - 2 * self.M * r + self.a**2

        """ Handles division by zero that happens at the singularity. """
        dr_g_tt = -2 * self.M / Sigma * (2 * r**2 / Sigma - 1) if Sigma else -infty

        """ Handles division by zero that happens at the horizons. """
        dr_g_rr = 2 * r / Delta * (1 - Sigma / Delta * (1 - self.M / r)) if Delta else infty

        dr_g_thth = 2 * r

        """ Handles division by zero that happens at the singularity. """
        dr_g_phiphi = 2 * (r - self.M * self.a**2 / Sigma * (2 * r**2 / Sigma - 1) * sin(theta)**2) * sin(theta)**2 if Sigma else infty

        """ Handles division by zero that happens at the singularity. """
        dr_g_tphi = 2 * self.M * self.a * sin(theta)**2 / Sigma * (2 * r**2 / Sigma - 1) if Sigma else sign(self.a) * infty

        dr_Metric = diag([dr_g_tt, dr_g_rr, dr_g_thth, dr_g_phiphi])
        dr_Metric[Coords.e_t.value][Coords.e_phi.value] = dr_g_tphi
        dr_Metric[Coords.e_phi.value][Coords.e_t.value] = dr_g_tphi
        
        return dr_Metric
    
    def get_photon_orbits(self) -> tuple[float, float]:

        r_prograde = 2 * self.M * (1 + cos(2 / 3 * arccos(self.a / self.M)))
        r_retrograde = 2 * self.M * (1 + cos(2 / 3 * arccos(-self.a / self.M)))

        return r_prograde, r_retrograde
    
    def get_ISCO(self) -> tuple[float, float]:

        Z_1 = 1 + pow(1 - (self.a / self.M)**2, 1. / 3) * (pow(1 + self.a / self.M, 1. / 3) + pow(1 - self.a / self.M, 1 / 3))
        Z_2 = sqrt(3 * (self.a / self.M)**2 + Z_1**2)

        r_prograde   = self.M * (3 + Z_2 - sqrt((3 - Z_1) * (3 + Z_1 + 2 * Z_2)))
        r_retrograde = self.M * (3 + Z_2 + sqrt((3 - Z_1) * (3 + Z_1 + 2 * Z_2)))

        return r_prograde, r_retrograde

class Regular_Black_Hole(Spacetime):

    def __init__(self, Param: float) -> None:
        self.M = 1
        self.Param = Param

    def get_metric(self, r: float, theta: float) -> NDArray:

        r_eff = sqrt(r**2 + self.Param**2)

        g_tt = -(1 - 2 * self.M / r_eff)

        """ Handles division by zero that happens at the horizons. """
        g_rr = -1 / g_tt if g_tt else -infty

        g_thth   = r_eff**2
        g_phiphi = g_thth * sin(theta)**2

        return diag([g_tt, g_rr, g_thth, g_phiphi])
    
    def get_dr_metric(self, r: float, theta: float) -> NDArray:

        r_eff = sqrt(r**2 + self.Param**2)

        g_tt = -(1 - 2 * self.M / r_eff)

        dr_g_tt = -2 * self.M * r / r_eff**3

        """ Handles division by zero that happens at the horizons. """
        dr_g_rr = -1 / g_tt**2 * dr_g_tt if g_tt else infty

        dr_g_thth   = 2 * r
        dr_g_phiphi = 2 * r * sin(theta)**2

        return diag([dr_g_tt, dr_g_rr, dr_g_thth, dr_g_phiphi])
        
    def get_photon_orbits(self) -> tuple[float, float]:

        r_prograde = sqrt((3 * self.M)**2 - self.Param**2)
        r_retrograde = r_prograde

        return r_prograde, r_retrograde

    def get_ISCO(self) -> tuple[float, float]:

        r_prograde = sqrt((6 * self.M)**2 - self.Param**2)
        r_retrograde = r_prograde

        return r_prograde, r_retrograde
    
        
    def Identify(self) -> int:
        
        return self.Spacetime_enums.Regular_Black_Hole.value

class Wormhole(Spacetime):

    def __init__(self, r_throat: float, parameter: float, use_global_coords: bool = True) -> None:
        self.M = 1
        self.R_throat = r_throat
        self.Param = parameter
        self.Use_global_coords = use_global_coords

    def __convert_global_to_local_coords__(self, Ell, p) -> tuple[float, float]:

        dEll_dr = sqrt(Ell**2 + self.R_throat**2) / Ell
        
        p[0] = p[0] * dEll_dr
        
        return p
        
    def get_metric(self, r: float, theta: float) -> NDArray:

        if self.Use_global_coords:
                
            """  The metric uses global coordinates - the range is [-inf, + inf]  """      
            Ell = r
            r = sqrt(Ell**2 + 1)

        g_tt = -exp(-2 * self.M / r - 2 * self.Param * (self.M / r)**2)

        """ Handles division by zero that happens at the throat. """
        g_rr = 1 / (1 - self.R_throat / r) if  1 / (1 - self.R_throat / r) else infty
        g_thth   = r**2
        g_phiphi = g_thth * sin(theta)**2

        return diag([g_tt, g_rr, g_thth, g_phiphi])
    
    def get_dr_metric(self, r: float, theta: float) -> NDArray:
        
        if self.Use_global_coords:
                
            """  The metric uses global coordinates - the range is [-inf, + inf]  """      
            Ell = r
            r = sqrt(Ell**2 + 1)

        dr_g_tt = -2 * exp(-2 * self.M / r - 2 * self.Param * (self.M / r)**2) * (self.M / r**2 + 2 * self.Param * self.M**2 / r**3)

        """ Handles division by zero that happens at the throat. """
        dr_g_rr = -1 / (1 - self.R_throat / r)**2 * (self.R_throat / r**2) if  1 / (1 - self.R_throat / r) else infty
        dr_g_thth   = 2 * r
        dr_g_phiphi = dr_g_thth * sin(theta)**2

        return diag([dr_g_tt, dr_g_rr, dr_g_thth, dr_g_phiphi])

    def get_photon_orbits(self):

        r_prograde = self.M / 2 * (1 + sqrt(1 + 8 * self.Param))

        if self.Use_global_coords:
            r_prograde = sqrt(r_prograde**2 - self.R_throat**2)

        r_retrograde = r_prograde

        return r_prograde, r_retrograde

    def get_ISCO(self):

        r_prograde = 2 * self.M * (sqrt(4 / 9 * (6 * self.Param + 1) ) * cosh(1 / 3 * arccosh( (1 + 9 * self.Param + 27 / 2 * self.Param**2) / pow(6 * self.Param + 1, 3 / 2))) + 1 / 3)
        
        if self.Use_global_coords:
            r_prograde = sqrt(r_prograde**2 - self.R_throat**2)

        r_retrograde = r_prograde

        return r_prograde, r_retrograde
    
    def Identify(self) -> int:
        
        return self.Spacetime_enums.Wormhole.value

class Janis_Newman_Winicour(Spacetime):

    def __init__(self, parameter):
        self.M = 1
        self.Param = parameter

    def get_metric(self, r: float, theta: float) -> NDArray:

        r_singularity = 2 * self.M / self.Param

        g_tt = -pow(1 - r_singularity / r, self.Param)

        """ Handles division by zero that happens at the horizons. """
        g_rr = -1 / g_tt if g_tt else infty

        g_thth   = r**2 * pow(1 - r_singularity / r, 1 - self.Param) 
        g_phiphi = g_thth * sin(theta)**2

        return diag([g_tt, g_rr, g_thth, g_phiphi])
    
    def get_dr_metric(self, r: float, theta: float) -> NDArray:

        r_singularity = 2 * self.M / self.Param

        g_tt = -pow(1 - r_singularity / r, self.Param)

        dr_g_tt = -self.Param * pow(1 - r_singularity / r, self.Param - 1) * (r_singularity / r**2)

        """ Handles division by zero that happens at the horizons. """
        dr_g_rr = - 1 / g_tt**2 * dr_g_tt if 1 / g_tt else infty

        dr_g_thth   = 2 * r * pow(1 - r_singularity / r, 1 - self.Param) + pow(1 - r_singularity / r, -self.Param) * r_singularity * (1 - self.Param)
        dr_g_phiphi = dr_g_thth * sin(theta)**2

        return diag([dr_g_tt, dr_g_rr, dr_g_thth, dr_g_phiphi])

    def get_photon_orbits(self) -> tuple[float, float]:

        if self.Param > 0.5:

            r_prograde = self.M / self.Param * (2 * self.Param + 1)
            r_retrograde = r_prograde

            return r_prograde, r_retrograde

        else:

            print("Photon sphere does not exist for values of Param < 0.5!")
            exit(1)

    def get_ISCO(self) -> tuple[float, float]:

        if self.Param > 1 / sqrt(5):

            r_ISCO_outer = 1 / self.Param * (3 * self.Param + 1 + sqrt(5 * self.Param**2 - 1))
            r_ISCO_inner = 1 / self.Param * (3 * self.Param + 1 - sqrt(5 * self.Param**2 - 1))

        else:

            r_ISCO_inner = 2 * self.M / self.Param
            r_ISCO_outer = r_ISCO_inner

        return r_ISCO_outer, r_ISCO_inner
    
    def Identify(self) -> int:
        
        return self.Spacetime_enums.Janis_Newman_Winicour.value
        
class Gauss_Bonnet(Spacetime):

    def __init__(self, Param: float):

        self.M = 1
        self.Param = Param
 
        if self.Param < 3 / 4 * sqrt(3):

            self.HAS_PHOTON_SPHERE = True

        else:

            self.HAS_PHOTON_SPHERE = False

        self.HAS_R_OF_B = True

    def get_metric(self, r: float, theta: float) -> NDArray:

        if self.Param != 0:
            f = 1 + r**2 / 2 / self.Param * (1 - sqrt(1 + 8 * self.Param * self.M / r**3))
        else:
            f = 1 - 2 * self.M / r

        g_tt = -f

        """ Handles division by zero that happens at the horizons. """
        g_rr = - 1 / g_tt if g_tt else infty

        g_thth = r**2
        g_phiphi = g_thth * sin(theta)**2

        return diag([g_tt, g_rr, g_thth, g_phiphi])
    
    def get_dr_metric(self, r: float, theta: float) -> NDArray:

        if self.Param != 0:
            f = 1 + r**2 / 2 / self.Param * (1 - sqrt(1 + 8 * self.Param * self.M / r**3))
            dr_f = 2. / r * (f - 1) + 6. * self.M / sqrt(r**4 + 8 * self.Param * self.M * r)
        else:
            f = 1 - 2 * self.M / r
            dr_f = 2 * self.M / r**2

        g_tt = -f
        dr_g_tt = -dr_f

        """ Handles division by zero that happens at the horizons. """
        dr_g_rr = 1 / g_tt**2 * dr_g_tt if g_tt else infty

        dr_g_thth = 2 * r
        dr_g_phiphi = dr_g_thth * sin(theta)**2

        return diag([dr_g_tt, dr_g_rr, dr_g_thth, dr_g_phiphi])
    
    def get_photon_orbits(self) -> tuple[float, float]:

        if self.Param > 3 / 4 * sqrt(3):

            print("No photon sphere exists for Param > 3 / 4 * sqrt(3)!")
            exit(1)

        polynomial_coefs = [1, 0, -9 * self.M**2, 8 * self.M * self.Param]

        pol_roots = roots(polynomial_coefs)
        pol_roots = pol_roots[pol_roots > 0]

        r_outer = max(pol_roots)
        r_inner = min(pol_roots)

        return r_outer, r_inner
    
    def Identify(self) -> int:
        
        return self.Spacetime_enums.Gauss_Bonnet.value