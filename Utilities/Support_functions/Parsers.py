import csv
import numpy as np
from itertools import islice

class Simulation_Parser():

    def __init__(self, File_name: str) -> None:

        with open(File_name + ".txt", 'r') as file:

            csvreader = csv.reader(file, delimiter = ":")

            _ = csvreader.__next__()

            self.metric = str(csvreader.__next__()[1]).split(" ")[1]

            _ = csvreader.__next__()

            self.OBS_DISTANCE    = float(csvreader.__next__()[1])
            self.OBS_INCLICATION = float(csvreader.__next__()[1])
            self.OBS_FREQUENCY   = float(csvreader.__next__()[1])

            for i in range(1):
                _ = csvreader.__next__()

            self.disk_profile = str(csvreader.__next__()[1])
            self.disk_profile = self.disk_profile.split(" ")[1] + " " + self.disk_profile.split(" ")[2]
            
            if self.disk_profile == "Power law":
                self.disk_opening_angle = float(csvreader.__next__()[1])
                self.R_0 = float(csvreader.__next__()[1])
                self.R_Cutoff = float(csvreader.__next__()[1])
                self.R_Cutoff_Scale = float(csvreader.__next__()[1])

            else:
                self.height_scale   = float(csvreader.__next__()[1])
                self.radial_scale   = float(csvreader.__next__()[1])


            self.emission_model = str(csvreader.__next__()[1])

            if self.emission_model == " Phenomenological":

                for i in range(3):
                    _ = csvreader.__next__()

                self.Emission_Scale = csvreader.__next__()[1]

                for i in range(7):
                    _ = csvreader.__next__()

            else:

                for i in range(11):
                    _ = csvreader.__next__()

            self.WINDOW_LIMITS = [float(limit) for limit in csvreader.__next__()[1].split(',')]
            
            Resolution_list = csvreader.__next__()[1].split(' ')

            self.X_PIXEL_COUNT   = int(Resolution_list[1])
            self.Y_PIXEL_COUNT   = int(Resolution_list[3])

            self.params = csvreader.__next__()

            self.Legend = csvreader.__next__()

            csvreader = csv.reader(file, delimiter = " ")

            self.X_coords        = np.zeros(self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT)
            self.Y_coords        = np.zeros(self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT)
            self.NT_Flux         = np.zeros(self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT)
            self.Intensity       = np.zeros(self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT)
            self.NT_Redshift     = np.zeros(self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT)
            self.NT_Flux_Shifted = np.zeros(self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT)

            index = 0

            for row in csvreader:

                try:

                    self.X_coords[index] = row[0]
                    self.Y_coords[index] = row[1]

                    self.NT_Redshift[index]  = row[2]
                    self.NT_Flux[index]      = row[3]
                    self.Intensity[index]    = row[4]

                    self.NT_Flux_Shifted[index] = self.NT_Redshift[index]**4*self.NT_Flux[index]

                    index += 1

                except:
                    break


    def get_total_flux(self, obs_pos):

        dx = np.abs(self.X_coords[0] - self.X_coords[1])
        dy = np.abs(self.Y_coords[0] - self.Y_coords[self.X_PIXEL_COUNT + 1])

        return np.sum(self.Intensity) * dx * dy / obs_pos**2

    def get_plottable_sim_data(self) -> tuple:

        # Arrays need to be flipped, because mpl treats y = 0 as the top, 
        # and the simulator (aka openGL) treats it as the bottom

        Intensity = self.Intensity.reshape(self.X_PIXEL_COUNT, self.Y_PIXEL_COUNT)
        Intensity = np.flip(Intensity, 0)

        NT_Flux         = self.NT_Flux.reshape(self.X_PIXEL_COUNT,self.Y_PIXEL_COUNT)
        NT_Flux         = np.flip(NT_Flux, 0)

        NT_Redshift     = self.NT_Redshift.reshape(self.X_PIXEL_COUNT,self.Y_PIXEL_COUNT)
        NT_Redshift     = np.flip(NT_Redshift, 0)

        NT_Flux_Shifted = self.NT_Flux_Shifted.reshape(self.X_PIXEL_COUNT,self.Y_PIXEL_COUNT)
        NT_Flux_Shifted = np.flip(NT_Flux_Shifted, 0)

        return Intensity, NT_Flux, NT_Redshift, NT_Flux_Shifted
    
    def export_ehtim_data(self, Spacetime: str, data: np.array, path: str):

        ehtim_x_fov = 2 * 5.500000e-05
        ehtim_y_fov = 2 * 5.500000e-05

        Units = Units_class()
        
        dx = np.abs(self.X_coords[0] - self.X_coords[1])
        dy = np.abs(self.Y_coords[0] - self.Y_coords[self.X_PIXEL_COUNT + 1])

        stupid_eht_scaling = 1

        formatted_sim_data = data.reshape(1, self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT).flatten()

        formatted_sim_data = formatted_sim_data

        array_to_export = np.array([self.X_coords / max(self.X_coords) * ehtim_x_fov / 2, 
                                    self.Y_coords / max(self.Y_coords) * ehtim_y_fov / 2, 
                                    formatted_sim_data * dx * dy / Units.M87_DISTANCE_GEOMETRICAL**2 * stupid_eht_scaling]).T

        header = ("SRC: M87 \n"                   + 
                  "RA: 12 h 30 m 49.3920 s \n"    +
                  "DEC: 12 deg 23 m 27.9600 s \n" +
                  "MJD: 58211.000000 \n"          + 
                  "RF: {} GHz \n".format(self.OBS_FREQUENCY / 1e9)    +
                  "FOVX: {} pix 0.000110 as \n".format(self.X_PIXEL_COUNT) +
                  "FOVY: {} pix 0.000110 as \n".format(self.Y_PIXEL_COUNT) +
                  "------------------------------------ \n" +
                  "x (as)     y (as)       I (Jy/pixel)")

        with open(path + '{}_data_for_ehtim_{}.csv'.format(Spacetime, int(self.OBS_FREQUENCY / 1e9)), 'w') as my_file:
            np.savetxt(my_file, array_to_export, fmt = '%0.4e', header = header)

        print('Array exported to file!')

class ehtim_Parser():

    def __init__(self, File_name: str) -> None:

        with open(File_name + ".txt", 'r') as file:

            self.HEADER_ROW_COUNT = 9

            csvreader = csv.reader(file, delimiter = " ")

            for i in range(4):
                    _ = csvreader.__next__()

            self.OBS_FREQUENCY = float(csvreader.__next__()[2])

            X_data_line = csvreader.__next__()
    
            self.X_PIXEL_COUNT = int(X_data_line[2])
            X_range            = float(X_data_line[4]) / 2

            Y_data_line = csvreader.__next__()

            self.Y_PIXEL_COUNT = int(Y_data_line[2])
            Y_range            = float(Y_data_line[4]) / 2

            self.WINDOW_LIMITS = [-X_range, X_range, -Y_range, Y_range]

            for i in range(2):
                    _ = csvreader.__next__()


            self.X_coords        = np.zeros(self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT)
            self.Y_coords        = np.zeros(self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT)
            self.Intensity       = np.zeros(self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT)

            index = 0

            for row in csvreader:

                self.X_coords[index]  = float(row[0])
                self.Y_coords[index]  = float(row[1])
                self.Intensity[index] = float(row[2])

                index += 1

    def get_plottable_ehtim_data(self) -> tuple:

        Intensity = self.Intensity.reshape(self.X_PIXEL_COUNT, self.Y_PIXEL_COUNT)

        Metadata = self.WINDOW_LIMITS
        
        return Intensity, Metadata
    
    def get_total_flux(self):

        return np.sum(self.Intensity)

class VIDA_params_Parser():

    def __init__(self, File_name):

        with open(File_name + ".csv", 'r') as file:

            csvreader = csv.reader(file, delimiter = " ")

            self.template_params = {"Gaussian_1":{}, 
                                    "Gaussian_2":{}}
            
            self.template_params["Gaussian_1"]["d0"]          = 2 * float(csvreader.__next__()[0])
            self.template_params["Gaussian_1"]["sigma"]       = float(csvreader.__next__()[0]) 
            self.template_params["Gaussian_1"]["tau"]         = float(csvreader.__next__()[0])
            self.template_params["Gaussian_1"]["rot_angle"]   = float(csvreader.__next__()[0])
            self.template_params["Gaussian_1"]["slash"]       = float(csvreader.__next__()[0])
            self.template_params["Gaussian_1"]["slash_angle"] = float(csvreader.__next__()[0])
            self.template_params["Gaussian_1"]["x0"]          = float(csvreader.__next__()[0])
            self.template_params["Gaussian_1"]["y0"]          = float(csvreader.__next__()[0])

            try:

                self.template_params["Gaussian_2"]["d0"]          = 2 * float(csvreader.__next__()[0])
                self.template_params["Gaussian_2"]["sigma"]       = float(csvreader.__next__()[0]) 
                self.template_params["Gaussian_2"]["tau"]         = float(csvreader.__next__()[0])
                self.template_params["Gaussian_2"]["rot_angle"]   = float(csvreader.__next__()[0])
                self.template_params["Gaussian_2"]["slash"]       = float(csvreader.__next__()[0])
                self.template_params["Gaussian_2"]["slash_angle"] = float(csvreader.__next__()[0])
                self.template_params["Gaussian_2"]["x0"]          = float(csvreader.__next__()[0])
                self.template_params["Gaussian_2"]["y0"]          = float(csvreader.__next__()[0])

            except:

                self.template_params["Gaussian_2"] = None

class Units_class():

    def __init__(self) -> None:

        #============== Useful scaling constants ==============#

        self.KILO = 1e3
        self.MEGA = 1e6
        self.GIGA = 1e9

        #============== Physical constants ==============#
        
        self.C_LIGHT_SI   = 299792458
        self.G_NEWTON_SI  = 6.6743e-11
        self.BOLTZMANN_SI = 1.380649e-23
        self.PLANCK_SI    = 6.62607015e-34

        self.M_SUN_SI     = 1.989e30
        self.M_M87_BH_SI  = 6.2e9 * self.M_SUN_SI
        self.M_SGRA_BH_SI = 4.154e6 * self.M_SUN_SI

        #==============  Time conversions  ==============#

        self.YEAR_TO_SEC = 31556952

        #============== Angular conversions ==============#

        self.RAD_TO_DEG      = 180 / np.pi
        self.ARCSEC_TO_RAD   = np.pi / 180 / 3600
        self.DEG_TO_AS       = 3600
        self.RAD_TO_MICRO_AS = self.RAD_TO_DEG * self.DEG_TO_AS * 1e6

        #============== Distance conversions ==============#

        self.LY_TO_METER      = self.C_LIGHT_SI * self.YEAR_TO_SEC
        self.GR_MASS_TO_METER = self.G_NEWTON_SI / self.C_LIGHT_SI**2

        self.M87_DISTANCE_LY = 53.49e6
        self.M87_DISTANCE_GEOMETRICAL = self.M87_DISTANCE_LY * self.LY_TO_METER / self.GR_MASS_TO_METER / self.M_M87_BH_SI

        self.SGRA_DISTANCE_LY = 26673
        self.SGRA_DISTANCE_GEOMETRICAL = self.SGRA_DISTANCE_LY * self.LY_TO_METER / self.GR_MASS_TO_METER / self.M_SGRA_BH_SI

        #============== Flux conversions ==============#

        self.W_M2_TO_JY = 1e26
        self.J_TO_ERG   = 1e7
    
    def Spectral_density_to_T(self, I_nu, f):

        I_nu += 1e-10 # To avoid division by 0 errors

        return self.PLANCK_SI * f / self.BOLTZMANN_SI / np.log(1 + 2 * self.PLANCK_SI * f**3 / self.C_LIGHT_SI**2 / I_nu)