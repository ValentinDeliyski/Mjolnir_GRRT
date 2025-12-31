from numpy import array, flip, append, pi, exp, sin, sqrt, log, roots, roll
import itertools

from Support_functions.Surface_Cubic_B_spline import Surface_Cubic_B_spline
import matplotlib.pyplot as plt

import xml.etree.cElementTree as ET
import xml.dom.minidom

from numpy.typing import NDArray

class Numerical_metric_parser_class():

    def __init__(self, File_path: str) -> None:
        
        """ ====================== Initialize the arrays that hold the metric functions ====================== """
        
        GRID_R_SIZE = 91
        GRID_THETA_SIZE = 51
        
        x_coord = []
        theta_coord = []
        F_0 = []
        F_1 = []
        F_2 = []
        W = []

        with open(File_path, "r") as file:
            
            for line in file:
                
                Line_contents = line.strip().split(" ")
                Line_contents = [x for x in Line_contents if "" != x]
                
                """ Different theta values on the grid are seperated by a "\n" character - parsing it results in an empty list. """
                if len(Line_contents) == 0:
                    continue

                x_coord.append(float(Line_contents[0]))
                theta_coord.append(float(Line_contents[1]))
                F_0.append(exp(2 * float(Line_contents[2])))
                F_1.append(exp(2 * float(Line_contents[3])))
                F_2.append(exp(2 * float(Line_contents[4])))
                W.append(-float(Line_contents[5]))
                
        """ The metric is calculated only for theta values in the range [0, pi / 2]. We use the reflection symmetry of the problem to get the rest of the grid. """
        x_coord = array(x_coord).reshape(GRID_THETA_SIZE, GRID_R_SIZE)
        self.x_coord = append(x_coord, x_coord[1:], axis = 0)
        
        theta_coord = array(theta_coord).reshape(GRID_THETA_SIZE, GRID_R_SIZE)
        self.theta_coord = append(theta_coord, flip(pi - theta_coord)[1:], axis = 0)
         
        F_0 = array(F_0).reshape(GRID_THETA_SIZE, GRID_R_SIZE)  
        self.F_0 = append(F_0, flip(F_0, axis = 0)[1:], axis = 0)  
        
        F_1 = array(F_1).reshape(GRID_THETA_SIZE, GRID_R_SIZE)  
        self.F_1 = append(F_1, flip(F_1, axis = 0)[1:], axis = 0)  
        
        F_2 = array(F_2).reshape(GRID_THETA_SIZE, GRID_R_SIZE) 
        self.F_2 = append(F_2, flip(F_2, axis = 0)[1:], axis = 0)  
        
        W = array(W).reshape(GRID_THETA_SIZE, GRID_R_SIZE)
        self.W = append(W, flip(W, axis = 0)[1:], axis = 0)

    def get_parsed_results(self) -> tuple[NDArray, NDArray, NDArray, NDArray, NDArray, NDArray]:
        
        return self.x_coord, self.theta_coord, self.F_0, self.F_1, self.F_2, self.W
        
    def export_spline_to_XML(self, Metric_name: str, Radial_control_vectors: list[NDArray], Theta_control_vectors: list[NDArray], Metric_control_vectors: list[NDArray], Control_vector_order: list[str]):
        
        X_knot_points, Theta_knot_points, _, _, _, _ = Numerical_metric_parser.get_parsed_results()
        
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

        with open("{}.XML".format(Metric_name), 'w') as xfile:
            xfile.write(Header + 'encoding=\"{}\"?>\n'.format(Encoding) + Body)
            xfile.close()
            
if __name__ == "__main__":
    
    Numerical_metric_parser = Numerical_metric_parser_class("Numerical_metrics/Zero_curvature/rh=0.1_om=0.738499261269268_h0=0.0594286954018100.dat")
    x_coord, theta_coord, F_0, F_1, F_2, W = Numerical_metric_parser.get_parsed_results()

    F_0_spline_instance = Surface_Cubic_B_spline(x_grid = theta_coord, y_grid = x_coord, z_grid = F_0, X_patch_number = 101, Y_patch_number = 91)
    F_1_spline_instance = Surface_Cubic_B_spline(x_grid = theta_coord, y_grid = x_coord, z_grid = F_1, X_patch_number = 101, Y_patch_number = 91)
    F_2_spline_instance = Surface_Cubic_B_spline(x_grid = theta_coord, y_grid = x_coord, z_grid = F_2, X_patch_number = 101, Y_patch_number = 91)
    W_spline_instance = Surface_Cubic_B_spline(x_grid = theta_coord, y_grid = x_coord, z_grid = W, X_patch_number = 101, Y_patch_number = 91)

    # exit(0)
    Numerical_metric_parser.export_spline_to_XML(Metric_name = "Galin_zero_curvature_config_IV", 
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
    
    Theta_surface, Radial_surface, F_0_surface = F_0_spline_instance.evaluate_spline(Patch_discretization = 15)
    Theta_surface, Radial_surface, F_1_surface = F_1_spline_instance.evaluate_spline(Patch_discretization = 15)
    Theta_surface, Radial_surface, F_2_surface = F_2_spline_instance.evaluate_spline(Patch_discretization = 15)
    Theta_surface, Radial_surface, W_surface = W_spline_instance.evaluate_spline(Patch_discretization = 15)
    
    _, (ax1, ax2, ax3, ax4) = plt.subplots(ncols = 4, nrows = 1, subplot_kw = dict(projection = '3d'))
    
    ax1.plot_surface(Radial_surface, Theta_surface, F_0_surface, color = 'orange') 
    ax1.plot_surface(x_coord, theta_coord, F_0, color = 'blue') 
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z') # type: ignore
    # ax1.set_xlim(0,1)  # type: ignore
    # ax1.set_ylim(0,3.14)  # type: ignore
    
    ax2.plot_surface(Radial_surface, Theta_surface, F_1_surface, color = 'orange') 
    ax2.plot_surface(x_coord, theta_coord, F_1, color = 'blue')
    
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('z') # type: ignore
    # ax2.set_xlim(0,1)  # type: ignore
    # ax2.set_ylim(0,3.14)  # type: ignore
    
    ax3.plot_surface(Radial_surface, Theta_surface, F_2_surface, color = 'orange') 
    ax3.plot_surface(x_coord, theta_coord, F_2, color = 'blue')
    
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_zlabel('z') # type: ignore
    # ax3.set_xlim(0,1)  # type: ignore
    # ax3.set_ylim(0,3.14)  # type: ignore
    
    ax4.plot_surface(Radial_surface, Theta_surface, W_surface, color = 'orange') 
    ax4.plot_surface(x_coord, theta_coord, W, color = 'blue') 
    
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    ax4.set_zlabel('z') # type: ignore
    # ax4.set_xlim(0,1)  # type: ignore
    # ax4.set_ylim(0,3.14)  # type: ignore
    
    plt.show()
    