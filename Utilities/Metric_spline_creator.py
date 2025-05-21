from numpy import array, flip, append, pi, exp, sin, sqrt, log, roots

from Support_functions.Surface_Cubic_B_spline import Surface_Cubic_B_spline
import matplotlib.pyplot as plt

import xml.etree.cElementTree as ET
import xml.dom.minidom

class Numerical_metric_parser_class():

    def __init__(self, File_path: str, M_ADM, a_ADM, r_H) -> tuple[array, array, array, array, array, array]:
        
        self.M_ADM = M_ADM
        self.a_ADM = a_ADM
        self.r_H = r_H
        
        """ ====================== Initialize the arrays that hold the metric functions ====================== """
        
        GRID_R_SIZE = 251
        GRID_THETA_SIZE = 30
        
        x_coord = []
        theta_coord = []
        F_0 = []
        F_1 = []
        F_2 = []
        W   = []
        
        with open(File_path, "r") as file:
            
            for line in file:
                
                Line_contents = line.strip().split(" ")
                Line_contents = [x for x in Line_contents if "" != x]
                
                """ Different theta values on the grid are seperated by a "\n" character - parsing it results in an empty list. """
                if len(Line_contents) == 0:
                    continue
                
                # """ Skip parsing the compactified radial coordiante grid point at the horizon (a.e. at x = 0) -> its g_rr diverges there so I can't use it. """
                # if float(Line_contents[0]) == 0:
                #     continue
                
                # """ Skip parsing the compactified radial coordiante grid point at infintiy (a.e. at x = 1) -> its not a useful value. """
                # if float(Line_contents[0]) == 1:
                #     continue
                
                x_coord.append(float(Line_contents[0]))
                theta_coord.append(float(Line_contents[1]))
                F_1.append(float(Line_contents[2]))
                F_2.append(float(Line_contents[3]))
                F_0.append(float(Line_contents[4]))
                W.append(float(Line_contents[5]))
                
        """ The metric is calculated only for theta values in the range [0, pi / 2]. We use the reflection symmetry of the problem to get the rest of the grid. """
        x_coord = array(x_coord).reshape(GRID_THETA_SIZE, GRID_R_SIZE)
        self.x_coord = append(x_coord, x_coord[1:], axis = 0)
        
        theta_coord = array(theta_coord).reshape(GRID_THETA_SIZE, GRID_R_SIZE)
        self.theta_coord = append(theta_coord, theta_coord[1:] + pi / 2, axis = 0)
        
        F_0 = array(F_0).reshape(GRID_THETA_SIZE, GRID_R_SIZE)
        self.F_0 = append(F_0, flip(F_0[1:], axis = 0), axis = 0)
        
        F_1 = array(F_1).reshape(GRID_THETA_SIZE, GRID_R_SIZE)  
        self.F_1 = append(F_1, flip(F_1[1:], axis = 0), axis = 0)  
        
        F_2 = array(F_2).reshape(GRID_THETA_SIZE, GRID_R_SIZE) 
        self.F_2 = append(F_2, flip(F_2[1:], axis = 0), axis = 0)  
        
        W = array(W).reshape(GRID_THETA_SIZE, GRID_R_SIZE)
        self.W = append(W, flip(W[1:], axis = 0), axis = 0)
        
        """ Convert the compactified coordinate x to the (mass normalized) unbounded radial coordinate.
            NOTE: This is NOT in Boyer-Linguist coordinates. """
        self.r_coord = sqrt((self.x_coord / (1 - self.x_coord))**2 + self.r_H**2) / self.M_ADM
        
        """ Convert the uncompactified radial coordinate to the Boyer-Linguist radial coordinate. """
        self.r_BL_coord = self.r_coord + self.a_ADM**2 / (1 + sqrt(1 - self.a_ADM**2))

    def get_parsed_results(self) -> tuple[array, array, array, array, array, array]:
        
        return self.x_coord, self.r_coord, self.r_BL_coord, self.theta_coord, self.F_0, self.F_1, self.F_2, self.W
        
    def get_metric_functions(self):
                
        """ This is how the metric function W(r) scales with the BH mass - stems from the fact that rW should be dimensionless. """
        W = self.W * self.M_ADM
        N = 1 - (self.r_H / self.M_ADM) / self.r_coord

        g_tt     = -exp(2 * self.F_0) * N + exp(2 * self.F_2) * (W * self.r_coord * sin(self.theta_coord))**2
        g_tphi   = -exp(2 * self.F_2) * W * (self.r_coord * sin(self.theta_coord))**2
        g_rr     =  exp(2 * self.F_1) / N
        g_thth   =  exp(2 * self.F_1) * self.r_coord**2 
        g_phiphi =  exp(2 * self.F_2) * self.r_coord**2 * sin(self.theta_coord)**2
        
        return g_tt, g_tphi, g_rr, g_thth, g_phiphi
    
    def plot_metric_functions(self, radial_coord_to_plot: str = "BL", radial_coord_cutoff: float = 3.5) -> None:
        
        g_tt, g_tphi, g_rr, g_thth, g_phiphi = self.get_metric_functions()
        
        match radial_coord_to_plot:
            case "BL":
                radial_coord = self.r_BL_coord
            case _:
                radial_coord = self.x_coord
            
        idx = radial_coord < radial_coord_cutoff
        
        _, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(ncols = 5, nrows = 1, subplot_kw = dict(projection='3d'))
        
        ax1.plot_surface(*(x.reshape(60, int(len(g_tt[idx].flatten()) / 60)) for x in (radial_coord[idx], self.theta_coord[idx], g_tt[idx])))
        ax1.set_xlabel(r'r')
        ax1.set_ylabel(r'$\theta$')
        ax1.set_zlabel(r"$g_{tt}$")
        
        ax2.plot_surface(*(x.reshape(60, int(len(g_tphi[idx].flatten()) / 60)) for x in (radial_coord[idx], self.theta_coord[idx], g_tphi[idx])))
        ax2.set_xlabel(r'r')
        ax2.set_ylabel(r'$\theta$')
        ax2.set_zlabel(r"$g_{t\phi}$")
     
        ax3.plot_surface(*(x.reshape(60, int(len(g_rr[idx].flatten()) / 60)) for x in (radial_coord[idx], self.theta_coord[idx], log(g_rr[idx]))))
        ax3.set_xlabel(r'r')
        ax3.set_ylabel(r'$\theta$')
        ax3.set_zlabel(r"ln$g_{rr}$")
        
        ax4.plot_surface(*(x.reshape(60, int(len(g_thth[idx].flatten()) / 60)) for x in (radial_coord[idx], self.theta_coord[idx], g_thth[idx])))
        ax4.set_xlabel(r'r')
        ax4.set_ylabel(r'$\theta$')
        ax4.set_zlabel(r"$g_{\theta\theta}$")
        
        ax5.plot_surface(*(x.reshape(60, int(len(g_phiphi[idx].flatten()) / 60)) for x in (radial_coord[idx], self.theta_coord[idx], g_phiphi[idx])))
        ax5.set_xlabel(r'r')
        ax5.set_ylabel(r'$\theta$')
        ax5.set_zlabel(r"$g_{\phi\phi}$")

        plt.show()
        
    def export_spline_to_XML(self, Metric_name: str, Radial_control_vectors: array, Theta_control_vectors: array, Metric_control_vectors: array, Control_vector_order: list):
        
        X_knot_points, _, _, Theta_knot_points, _, _, _, _ = Numerical_metric_parser.get_parsed_results()
        
        Encoding = 'UTF-8'
        XML_root_node = ET.Element("Metric_spline_coefficients", {"Metric_name": Metric_name})
        
        Grid_subelement = ET.SubElement(XML_root_node, "Coordinate_grid")
        
        """ =============================== The radial grid knot points and control vector =============================== """
        """ All the coorinate control vectors SHOULD be identical, so I just intdex the one for the first provided metric component """
        
        X_coordinate_subelement = ET.SubElement(Grid_subelement, "Compactified_radial_coordinate_grid", units = "[-]")
        X_grid_knots_subelement = ET.SubElement(X_coordinate_subelement, "Grid_knots", Grid_size = "{}".format(len(X_knot_points[0])))
        for Grid_point_idx, X_point in enumerate(X_knot_points[0]):
            ET.SubElement(X_grid_knots_subelement, "Grid_point_idx_{}".format(Grid_point_idx)).text = "{}".format(X_point)
            
        X_grid_control_vector_subelement = ET.SubElement(X_coordinate_subelement, "Control_vector", Component_number = "{}".format(len(Radial_control_vectors[0])))
        for Component_idx, X_point in enumerate(Radial_control_vectors[0]):
            ET.SubElement(X_grid_control_vector_subelement, "Component_idx_{}".format(Component_idx)).text = "{}".format(X_point)
        
        """ =============================== The theta grid knot points and control vector =============================== """
        """ All the coorinate control vectors SHOULD be identical, so I just intdex the one for the first provided metric component """
        
        Theta_coordinate_subelement = ET.SubElement(Grid_subelement, "Theta_coordinate_grid", units = "[Rad]")
        Theta_grid_knots_subelement = ET.SubElement(Theta_coordinate_subelement, "Grid_knots", Grid_size = "{}".format(len(Theta_knot_points.T[0])))
        for Knot_idx, Theta_point in enumerate(Theta_knot_points.T[0]):
            ET.SubElement(Theta_grid_knots_subelement, "Grid_point_idx_{}".format(Knot_idx)).text = "{}".format(Theta_point)
               
        Theta_grid_control_vector_subelement = ET.SubElement(Theta_coordinate_subelement, "Control_vector", Component_number = "{}".format(len(Theta_control_vectors[0])))
        for Component_idx, Y_point in enumerate(Theta_control_vectors[0]):
            ET.SubElement(Theta_grid_control_vector_subelement, "Component_idx_{}".format(Component_idx)).text = "{}".format(Y_point)
        
        """ =============================== The metric component control vectors =============================== """

        for Metric_component_idx, Metric_component in enumerate(Control_vector_order):
            
            Metric_component_subelement = ET.SubElement(XML_root_node, Metric_component)
            
            Control_vector_z_element = ET.SubElement(Metric_component_subelement, "Control_vector", Component_number = "{}".format(len(Metric_control_vectors[Metric_component_idx])))
            for Vector_component_idx, Control_vec_component in enumerate(Metric_control_vectors[Metric_component_idx]):
                ET.SubElement(Control_vector_z_element, "Component_idx_{}".format(Vector_component_idx)).text = "{}".format(Control_vec_component)
                
        XML_struct = xml.dom.minidom.parseString(ET.tostring(XML_root_node))
        formatted_XML_string = XML_struct.toprettyxml()
        Header, Body = formatted_XML_string.split('?>')

        with open("test.XML", 'w') as xfile:
            xfile.write(Header + 'encoding=\"{}\"?>\n'.format(Encoding) + Body)
            xfile.close()
            
if __name__ == "__main__":
    
    """ =================== Some post-evolution calculated metric parameters =================== """
    
    """ This is a rounded value for the Black Hole mass - we take is for granted and from it, calculate what the spin parameter SHOULD be to give their reported event horizon raius. """
    M_ADM = 0.415
    
    """ Their reported event horizon radius - this is an input to the simulation, and they gave a decent amount of digits, so it should be fine.
        NOTE: This is NOT in Boyer-Linguist coordinates. """
    r_H = 0.0662902
    
    """ This is the (mass normalized) calculated spin parameter that gives their event horizon radius (its one of the real solutions to a quintic equation).
        Let alpha = 2M^2 - Mr_H and beta = 2M + r_h, then the quintic is 4 a^4 + (beta^2 - 4 alpha) a^2 + alpha^2 - M^2 beta^2 = 0. """
    a_ADM = 0.41399683 / M_ADM
    # a_ADM = 0.172 / M_ADM**2

    Numerical_metric_parser = Numerical_metric_parser_class("Numerical_metrics/configuration-II.dat", M_ADM = M_ADM, a_ADM = a_ADM, r_H = r_H)
    x_coord, _, r_BL_coord, theta_coord, F_0, F_1, F_2, W = Numerical_metric_parser.get_parsed_results()

    F_0_spline_instance = Surface_Cubic_B_spline(x_grid = theta_coord, y_grid = x_coord, z_grid = F_0, X_patch_number = 59, Y_patch_number = 251)
    F_1_spline_instance = Surface_Cubic_B_spline(x_grid = theta_coord, y_grid = x_coord, z_grid = F_1, X_patch_number = 59, Y_patch_number = 251)
    F_2_spline_instance = Surface_Cubic_B_spline(x_grid = theta_coord, y_grid = x_coord, z_grid = F_2, X_patch_number = 59, Y_patch_number = 251)
    W_spline_instance = Surface_Cubic_B_spline(x_grid = theta_coord, y_grid = x_coord, z_grid = W, X_patch_number = 59, Y_patch_number = 251)

    Numerical_metric_parser.export_spline_to_XML(Metric_name = "Numerical_Kerr_Config_II", 
                                                 Theta_control_vectors  = [F_0_spline_instance.Control_vector_X, 
                                                                           F_1_spline_instance.Control_vector_X, 
                                                                           F_2_spline_instance.Control_vector_X, 
                                                                           W_spline_instance.Control_vector_X], 
                                                 Radial_control_vectors = [F_0_spline_instance.Control_vector_Y, 
                                                                           F_1_spline_instance.Control_vector_Y, 
                                                                           F_2_spline_instance.Control_vector_Y, 
                                                                           W_spline_instance.Control_vector_Y], 
                                                 Metric_control_vectors = [F_0_spline_instance.Control_vector_Z, 
                                                                           F_1_spline_instance.Control_vector_Z, 
                                                                           F_2_spline_instance.Control_vector_Z, 
                                                                           W_spline_instance.Control_vector_Z], 
                                                 Control_vector_order   = ["F_0", "F_1", "F_2", "W"])
    
    Theta_surface, Radial_surface, F_0_surface = F_0_spline_instance.evaluate_spline(Patch_discretization = 5)
    Theta_surface, Radial_surface, F_1_surface = F_1_spline_instance.evaluate_spline(Patch_discretization = 5)
    Theta_surface, Radial_surface, F_2_surface = F_2_spline_instance.evaluate_spline(Patch_discretization = 5)
    Theta_surface, Radial_surface, W_surface = W_spline_instance.evaluate_spline(Patch_discretization = 5)
    
    _, (ax1, ax2, ax3, ax4) = plt.subplots(ncols = 4, nrows = 1, subplot_kw = dict(projection = '3d'))
    
    ax1.plot_surface(Radial_surface, Theta_surface, F_0_surface, color = 'orange') 
    ax1.plot_surface(x_coord, theta_coord, F_0, color = 'blue') 
    
    ax2.plot_surface(Radial_surface, Theta_surface, F_1_surface, color = 'orange') 
    ax2.plot_surface(x_coord, theta_coord, F_1, color = 'blue')
    
    ax3.plot_surface(Radial_surface, Theta_surface, F_2_surface, color = 'orange') 
    ax3.plot_surface(x_coord, theta_coord, F_2, color = 'blue')
    
    ax4.plot_surface(Radial_surface, Theta_surface, W_surface, color = 'orange') 
    ax4.plot_surface(x_coord, theta_coord, W, color = 'blue')
    
    plt.show()
    
    # Numerical_metric_parser.plot_metric_functions("x", radial_coord_cutoff = 10)