from Support_functions.Parsers import Units_class, Simulation_Parser
import matplotlib.pyplot as plt
import numpy as np
import os 
 
class Hotspot_analyzer():
    
    def __init__(self, sim_path: str, metric: str) -> None:
        
        parent_directory = os.path.abspath('...')
        
        self.Sim_parser_n0 = Simulation_Parser(parent_directory + sim_path + metric + "_n0")
        I_Intensity_0, Q_Intensity_0, _, V_Intensity_0, _, _, _ = self.Sim_parser_n0.get_plottable_sim_data()
                
        self.Sim_parser_n1 = Simulation_Parser(parent_directory + sim_path + metric + "_n1")
        I_Intensity_1, Q_Intensity_1, _, V_Intensity_1, _, _, _ = self.Sim_parser_n1.get_plottable_sim_data()
                
        self.Sim_parser_n2 = Simulation_Parser(parent_directory + sim_path + metric + "_n2")
        I_Intensity_2, Q_Intensity_2, _, V_Intensity_2, _, _, _ = self.Sim_parser_n2.get_plottable_sim_data()
                
        self.Sim_parser_n3 = Simulation_Parser(parent_directory + sim_path + metric + "_n3")
        I_Intensity_3, Q_Intensity_3, _, V_Intensity_3, _, _, _ = self.Sim_parser_n3.get_plottable_sim_data()
        
        self.Total_I_intensity = I_Intensity_0 + I_Intensity_1 + I_Intensity_2 + I_Intensity_3
        self.Total_Q_intensity = Q_Intensity_0 + Q_Intensity_1 + Q_Intensity_2 + Q_Intensity_3
        self.Total_V_intensity = V_Intensity_0 + V_Intensity_1 + V_Intensity_2 + V_Intensity_3
        
    def compute_centroid(self) -> list:
        
        _1D_I_intensity_array = self.Sim_parser_n0.I_Intensity + self.Sim_parser_n1.I_Intensity + self.Sim_parser_n2.I_Intensity + self.Sim_parser_n3.I_Intensity
        Image_X_coords = self.Sim_parser_n0.X_coords
        Image_Y_coords = self.Sim_parser_n0.Y_coords
        
        Intensity_sum = sum(_1D_I_intensity_array)
        X_centroid = 0
        Y_centroid = 0
        
        for Intensity, X_coords, Y_coords in zip(_1D_I_intensity_array, Image_X_coords, Image_Y_coords):
            
            X_centroid = X_centroid + X_coords * Intensity / Intensity_sum
            Y_centroid = Y_centroid + Y_coords * Intensity / Intensity_sum
        
        return [X_centroid, Y_centroid]
    
    
if __name__ == "__main__":
    
    Centroid_X = []
    Centroid_Y = []
    
    for run in range(10):
    
        Hotspot_analyzer_instance = Hotspot_analyzer("Reference_simulations\\Hotspot_Reference_Simulation_{}\\".format(run), "Kerr")
        Hotspot_centroid = Hotspot_analyzer_instance.compute_centroid()
        
        Centroid_X.append(Hotspot_centroid[0])
        Centroid_Y.append(Hotspot_centroid[1])
        
    plt.plot(Centroid_X, Centroid_Y)
    plt.show()