from Support_functions.Spacetimes_new import Kerr

from enum import Enum
from numpy.typing import NDArray
from numpy import pi, array, flip, append, linspace, meshgrid, column_stack, isnan, log10, cos, sin, argsort, stack, sqrt, mean, zeros_like
from scipy.interpolate import CloughTocher2DInterpolator, LinearNDInterpolator, interp1d, RegularGridInterpolator
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt

import xml.etree.cElementTree as ET
import xml.dom.minidom

class Coords(Enum):

    e_t = 0
    e_r = 1
    e_theta = 2
    e_phi = 3

    e_x = 0
    e_y = 1
    e_z = 2

class Disk_model():
    
    def Parse_raw_density_file(self):
        
        raw_r_coord = []
        raw_theta_coord = []
        raw_Density = []

        with open(self.Density_file_path, "r") as file:

            file.__next__()

            line = file.__next__()
            Line_contents = line.strip().split(" ")

            self.r_grid_len = int(Line_contents[1])

            """ Half the theta grid, because I take only one side of the equator. """
            self.theta_grid_len = int(int(Line_contents[2]) / 2) + 1

            file.__next__()
            
            for line in file:
                
                Line_contents = line.strip().split(" ")
                Line_contents = [x for x in Line_contents if "" != x]
                
                """ Different theta values on the grid are seperated by a "\n" character - parsing it results in an empty list. """
                if (len(Line_contents) == 0 or cos(float(Line_contents[1])) < 0):
                    continue

                raw_r_coord.append(float(Line_contents[0]))
                raw_theta_coord.append(cos(float(Line_contents[1])))
                raw_Density.append(float(Line_contents[2]))    

            self.r_grid_len = int(raw_r_coord.__len__() / self.theta_grid_len)

            self.raw_r_coord: NDArray = array(raw_r_coord)
            self.raw_theta_coord: NDArray = array(raw_theta_coord)
            self.raw_Density: NDArray = array(raw_Density).reshape((self.theta_grid_len, self.r_grid_len))
            
            self.rho_coord = (self.raw_r_coord * sqrt(1 - self.raw_theta_coord**2)).reshape((self.theta_grid_len, self.r_grid_len))
            self.z_coord = (self.raw_r_coord * self.raw_theta_coord).reshape((self.theta_grid_len, self.r_grid_len))
        
    def interpolate_raw_data(self):
        
        self.r_coord_range = linspace(min(self.raw_r_coord), max(self.raw_r_coord), self.GRID_R_SIZE)
        self.theta_coord_range = linspace(min(self.raw_theta_coord), max(self.raw_theta_coord), self.GRID_THETA_SIZE)

        self.r_coord, self.theta_coord = meshgrid(self.r_coord_range, self.theta_coord_range, indexing = "ij")
        
        interp = RegularGridInterpolator((self.raw_theta_coord[::self.r_grid_len], self.raw_r_coord[:self.r_grid_len]), self.raw_Density, "nearest")
        
        self.Density = interp(stack([self.theta_coord, self.r_coord], axis = -1), method = "nearest")
        self.Density[isnan(self.Density)] = 0.0
        self.Density[self.Density < 0.0] = 0.0
        
        plt.figure(figsize = (8, 6)) 
        plt.contourf(self.rho_coord, self.z_coord, self.raw_Density, levels = 450) 
        plt.contour(self.rho_coord, self.z_coord, self.raw_Density, levels = [1e-5 * 5.211680725829852e-17], colors = "r")
        
        plt.show()
             
    def Export_interpolated_data_to_XML(self, Model_name: str):
        
        Encoding = 'UTF-8'
        XML_root_node = ET.Element("Numerical_disk_data", {"Disk_model_name": Model_name})
             
        Parameters_subelement = ET.SubElement(XML_root_node, "Parameters")
        ET.SubElement(Parameters_subelement, "Geometric_Central_Density").text = "{}".format(self.Density_at_center)
        ET.SubElement(Parameters_subelement, "Geometric_K_rho").text = "{}".format(self.K_rho_geom)
        ET.SubElement(Parameters_subelement, "Geometric_K_mag").text = "{}".format(self.K_mag_geom)
        ET.SubElement(Parameters_subelement, "Gamma_rho").text = "{}".format(self.Gamma_rho)
        ET.SubElement(Parameters_subelement, "Gamma_mag").text = "{}".format(self.Gamma_mag)
        
        Rho_coord_subelement = ET.SubElement(XML_root_node, "Rho_coord_grid", units = "[M]", Grid_size = "{}".format(len(self.r_coord_range)))
        
        for Grid_point_idx, Rho_coord in enumerate(self.r_coord_range):
            ET.SubElement(Rho_coord_subelement, "Grid_point_idx_{}".format(Grid_point_idx)).text = "{}".format(Rho_coord)
           
        Z_coord_subelement = ET.SubElement(XML_root_node, "Z_coord_grid", units = "[M]", Grid_size = "{}".format(len(self.theta_coord_range)))
        
        for Grid_point_idx, Z_coord in enumerate(self.theta_coord_range):
            ET.SubElement(Z_coord_subelement, "Grid_point_idx_{}".format(Grid_point_idx)).text = "{}".format(Z_coord)
            
        Density_subelement = ET.SubElement(XML_root_node, "Density", units = "[-]", Component_number = "{}".format(len(self.Density.flatten())))
        
        for Grid_point_idx, Density in enumerate(self.Density.T.flatten()):
            ET.SubElement(Density_subelement, "Component_idx_{}".format(Grid_point_idx)).text = "{}".format(Density)
                   
        XML_struct = xml.dom.minidom.parseString(ET.tostring(XML_root_node))
        formatted_XML_string = XML_struct.toprettyxml()
        Header, Body = formatted_XML_string.split('?>')

        with open("{}.XML".format(Model_name), 'w') as xfile:
            xfile.write(Header + 'encoding=\"{}\"?>\n'.format(Encoding) + Body)
            xfile.close()
        
    def __init__(self, Density_file_path: str, Grid_R_size: int, Grid_Theta_size: int) -> None:
        
        """ === The reference for this is https://journals.aps.org/prd/pdf/10.1103/PhysRevD.99.043002 === """
        
        self.Kerr_instance = Kerr(mass = 1, spin_param = 0.9)
        
        """ The radial coordinate of the disk center - from "physical_formatted_output_data.dat" """
        r_center = 3.3342120001288786
        
        """ This is the metric factor appearing in the magnetic pressure EOS - expression (14) in the reference """
        Metric_at_center = self.Kerr_instance.get_metric(r = r_center, theta = pi / 2)
        L_at_center = Metric_at_center[Coords.e_t.value][Coords.e_phi.value]**2 - Metric_at_center[Coords.e_t.value][Coords.e_t.value] * Metric_at_center[Coords.e_phi.value][Coords.e_phi.value]
        
        """ This is mentioned under section III - Methodology """
        self.Density_at_center = 5.2e-17
        
        """ The density polytropic index - mentioned under section III - Methodology """
        self.Gamma_rho = 4 / 3
        self.Gamma_mag = 4 / 3
        
        """ Gotten from "physical_formatted_output_data.dat """
        Gas_Pressure_at_center = 8.224680416168014e-19
        Mag_Pressure_at_center = 8.224680341623081e-20
        
        """ Compute the geometric polytropic coeff """
        self.K_rho_geom = Gas_Pressure_at_center / self.Density_at_center**self.Gamma_rho
        
        """ Compute the specific enthalpy at the disk center - Rizzolla (2.229) + (2.248), but also mentioned under expression (23) in the reference """
        Specific_int_energy_at_center = 0.04745009563391138
        Specific_enthalpy_at_center = 1 + self.Gamma_rho * Specific_int_energy_at_center
        
        self.K_mag_geom = Mag_Pressure_at_center / L_at_center**(self.Gamma_mag - 1) / (self.Density_at_center * Specific_enthalpy_at_center)**(self.Gamma_mag)
        
        self.GRID_R_SIZE = Grid_R_size
        self.GRID_THETA_SIZE = Grid_Theta_size
        
        self.Density_file_path = Density_file_path

Disk_model_instance = Disk_model("Numerical_disks/a0.9_alpha0.5_dW0.9_mag1e1_sub/BHAC_out_rescaled.dat", 150, 150)
Disk_model_instance.Parse_raw_density_file()
Disk_model_instance.interpolate_raw_data()
Disk_model_instance.Export_interpolated_data_to_XML("Test")
        
        