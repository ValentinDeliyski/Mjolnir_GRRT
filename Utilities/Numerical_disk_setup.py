from Support_functions.Spacetimes_new import Kerr

from enum import Enum
from numpy import pi, array, flip, append, linspace, meshgrid, column_stack, isnan
from scipy.interpolate import CloughTocher2DInterpolator
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
        
        self.raw_rho_coord = []
        self.raw_z_coord = []
        self.raw_Density = []

        with open(self.Density_file_path, "r") as file:
            
            for line in file:
                
                Line_contents = line.strip().split(" ")
                Line_contents = [x for x in Line_contents if "" != x]
                
                """ Different theta values on the grid are seperated by a "\n" character - parsing it results in an empty list. """
                if len(Line_contents) == 0:
                    continue

                self.raw_rho_coord.append(float(Line_contents[0]))
                self.raw_z_coord.append(float(Line_contents[1]))
                self.raw_Density.append(float(Line_contents[2]))    
                
            self.raw_Density = array(self.raw_Density)
        
    def interpolate_raw_data(self):
        
        self.rho_coord_range = linspace(min(self.raw_rho_coord), max(self.raw_rho_coord), self.GRID_RHO_SIZE)
        self.z_coord_range = linspace(min(self.raw_z_coord), max(self.raw_z_coord), self.GRID_Z_SIZE)

        self.rho_coord, self.z_coord = meshgrid(self.rho_coord_range, self.z_coord_range)
        
        interp = CloughTocher2DInterpolator(column_stack((self.raw_rho_coord, self.raw_z_coord)), self.raw_Density)
        
        self.Density = interp(self.rho_coord, self.z_coord)
        self.Density[isnan(self.Density)] = 0.0
        self.Density[self.Density < 0.0] = 0.0
        
        # self.rho_coord = self.rho_coord.reshape(self.GRID_RHO_SIZE, self.GRID_Z_SIZE)
        # self.rho_coord = append(self.rho_coord[:-1], self.rho_coord, axis = 0)
        
        # self.z_coord = self.z_coord.reshape(self.GRID_RHO_SIZE, self.GRID_Z_SIZE)
        # self.z_coord = append(-flip(self.z_coord)[:-1], self.z_coord, axis = 0)
         
        # self.Density = self.Density.reshape(self.GRID_RHO_SIZE, self.GRID_Z_SIZE)  
        # self.Density = append(flip(self.Density, axis = 0)[:-1], self.Density, axis = 0)
        
        plt.figure(figsize=(8, 6))
        cf = plt.contourf( self.rho_coord, self.z_coord, self.Density, levels=50)
        # plt.scatter(r, z, c='k', s=3, alpha=0.3)
        plt.xlabel("r")
        plt.ylabel("z")
        plt.colorbar(cf, label="rho")
        plt.tight_layout()
        plt.show()
        
    def Export_interpolated_data_to_XML(self, Model_name: str):
        
        Encoding = 'UTF-8'
        XML_root_node = ET.Element("Numerical_disk_data", {"Disk_model_name": Model_name})
             
        Rho_coord_subelement = ET.SubElement(XML_root_node, "Rho_coord_grid", units = "[M]", Grid_size = "{}".format(len(self.rho_coord_range)))
        
        for Grid_point_idx, Rho_coord in enumerate(self.rho_coord_range):
            ET.SubElement(Rho_coord_subelement, "Grid_point_idx_{}".format(Grid_point_idx)).text = "{}".format(Rho_coord)
           
        Z_coord_subelement = ET.SubElement(XML_root_node, "Z_coord_grid", units = "[M]", Grid_size = "{}".format(len(self.z_coord_range)))
        
        for Grid_point_idx, Z_coord in enumerate(self.z_coord_range):
            ET.SubElement(Z_coord_subelement, "Grid_point_idx_{}".format(Grid_point_idx)).text = "{}".format(Z_coord)
            
        Density_subelement = ET.SubElement(XML_root_node, "Density", units = "[-]", Component_number = "{}".format(len(self.Density.flatten())))
        
        for Grid_point_idx, Density in enumerate(self.Density.flatten()):
            ET.SubElement(Density_subelement, "Component_idx_{}".format(Grid_point_idx)).text = "{}".format(Density)
                   
        XML_struct = xml.dom.minidom.parseString(ET.tostring(XML_root_node))
        formatted_XML_string = XML_struct.toprettyxml()
        Header, Body = formatted_XML_string.split('?>')

        with open("{}.XML".format(Model_name), 'w') as xfile:
            xfile.write(Header + 'encoding=\"{}\"?>\n'.format(Encoding) + Body)
            xfile.close()
        
    def __init__(self, Density_file_path: str, Grid_Rho_size: int, Grid_Z_size: int) -> None:
        
        """ === The reference for this is https://journals.aps.org/prd/pdf/10.1103/PhysRevD.99.043002 === """
        
        self.Kerr_instance = Kerr(mass = 1, spin_param = 0.5)
        
        """ The radial coordinate of the disk center - from "physical_formatted_output_data.dat" """
        r_center = 6.7016292925899705
        
        """ This is the metric factor appearing in the magnetic pressure EOS - expression (14) in the reference """
        Metric_at_center = self.Kerr_instance.get_metric(r = r_center, theta = pi / 2)
        L_at_center = Metric_at_center[Coords.e_t.value][Coords.e_phi.value]**2 - Metric_at_center[Coords.e_t.value][Coords.e_t.value] * Metric_at_center[Coords.e_phi.value][Coords.e_phi.value]
        
        """ This is mentioned under section III - Methodology """
        Density_at_center = 1
        
        """ The density polytropic index - mentioned under section III - Methodology """
        self.Gamma_rho = 4 / 3
        self.Gamma_mag = 4 / 3
        
        """ Gotten from "physical_formatted_output_data.dat """
        Gas_Pressure_at_center = 0.004029623980013597
        Mag_Pressure_at_center = 0.0004029623864019536
        
        """ Compute the geometric polytropic coeff """
        self.K_rho_geom = Gas_Pressure_at_center / Density_at_center**self.Gamma_rho
        
        """ Compute the specific enthalpy at the disk center - Rizzolla (2.229) + (2.248), but also mentioned under expression (23) in the reference """
        Specific_int_energy_at_center = 0.012088921310206105
        Specific_enthalpy_at_center = 1 + self.Gamma_rho * Specific_int_energy_at_center
        
        self.K_mag_geom = Mag_Pressure_at_center / L_at_center**(self.Gamma_mag - 1) / (Density_at_center * Specific_enthalpy_at_center)**(self.Gamma_mag)
        
        self.GRID_RHO_SIZE = Grid_Rho_size
        self.GRID_Z_SIZE = Grid_Z_size
        
        self.Density_file_path = Density_file_path

Disk_model_instance = Disk_model("Numerical_disks/BL_density.dat", 150, 150)
Disk_model_instance.Parse_raw_density_file()
Disk_model_instance.interpolate_raw_data()
Disk_model_instance.Export_interpolated_data_to_XML("Test")
        
        