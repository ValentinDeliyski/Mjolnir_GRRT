import csv
import numpy as np
from itertools import islice

class Simulation_Parser():

    def __init__(self, File_name: str) -> None:

        with open(File_name + ".txt", 'r') as file:


            csvreader = csv.reader(file, delimiter = ":")

            _ = csvreader.__next__()

            self.metric = str(csvreader.__next__()[1])

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

        Metadata = self.OBS_DISTANCE, self.OBS_INCLICATION, self.WINDOW_LIMITS, self.Legend
        
        return Intensity, NT_Flux, NT_Redshift, NT_Flux_Shifted, Metadata
    
    def export_ehtim_data(self, data: np.array):

        ehtim_x_fov = 2 * 5.500000e-05
        ehtim_y_fov = 2 * 5.500000e-05

        Units = Units_class()
        
        dx = np.abs(self.X_coords[0] - self.X_coords[1])
        dy = np.abs(self.Y_coords[0] - self.Y_coords[self.X_PIXEL_COUNT + 1])

        stupid_eht_scaling = 1

        formatted_sim_data = data.reshape(1, self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT).flatten()

        array_to_export = np.array([self.X_coords / max(self.X_coords) * ehtim_x_fov / 2, 
                                    self.Y_coords / max(self.Y_coords) * ehtim_y_fov / 2, 
                                    formatted_sim_data * dx * dy / Units.M87_DISTANCE_GEOMETRICAL**2 * stupid_eht_scaling]).T

        print(sum(formatted_sim_data)* dx * dy / Units.M87_DISTANCE_GEOMETRICAL**2 * stupid_eht_scaling)

        header = ("SRC: M87 \n"                   + 
                  "RA: 12 h 30 m 49.3920 s \n"    +
                  "DEC: 12 deg 23 m 27.9600 s \n" +
                  "MJD: 58211.000000 \n"          + 
                  "RF: 230.0000 GHz \n"           +
                  "FOVX: {} pix 0.000150 as \n".format(self.X_PIXEL_COUNT) +
                  "FOVY: {} pix 0.000150 as \n".format(self.Y_PIXEL_COUNT) +
                  "------------------------------------ \n" +
                  "x (as)     y (as)       I (Jy/pixel)")

        with open('data_for_ehtim.csv', 'w') as my_file:
            np.savetxt(my_file, array_to_export, fmt = '%0.4e', header = header)

        print('Array exported to file!')

class ehtim_Parser():

    def __init__(self, File_name: str) -> None:

        def get_csv_line(path, line_number) -> str:
            with open(path) as file:
                return next(islice(csv.reader(file), line_number, None))[0]

        with open("..\\Sim_Results\\" + File_name + ".txt", 'r') as file:

            self.HEADER_ROW_COUNT = 9

            csvreader = csv.reader(file, delimiter = " ")
            row_count = sum(1 for row in file) - self.HEADER_ROW_COUNT
            file.seek(0)
    
            self.X_PIXEL_COUNT = int(get_csv_line("..\\Sim_Results\\" + File_name + ".txt", 5).split(" ")[2])
            self.Y_PIXEL_COUNT = int(get_csv_line("..\\Sim_Results\\" + File_name + ".txt", 6).split(" ")[2])

            X_range = float(get_csv_line("..\\Sim_Results\\" + File_name + ".txt", 5).split(" ")[4])
            Y_range = float(get_csv_line("..\\Sim_Results\\" + File_name + ".txt", 6).split(" ")[4])

            self.WINDOW_LIMITS = [-X_range, X_range, -Y_range, Y_range]

            self.Legend = get_csv_line("..\\Sim_Results\\" + File_name + ".txt", 8).split("  ")
            self.Legend = [self.Legend[0][2:8], self.Legend[2][1:7], self.Legend[5][1:13]]
            
            next(islice(csvreader, 10, None))

            self.X_coords        = np.zeros(row_count)
            self.Y_coords        = np.zeros(row_count)
            self.Intensity       = np.zeros(row_count)

            index = 0

            for row in csvreader:

                self.X_coords[index]  = row[0]
                self.Y_coords[index]  = row[1]
                self.Intensity[index] = row[2]

                index += 1

    def get_plottable_ehtim_data(self) -> tuple:

        # Arrays need to be flipped, because mpl treats y = 0 as the top, 
        # and the simulator (aka openGL) treats it as the bottom

        Intensity = self.Intensity.reshape(self.X_PIXEL_COUNT, self.Y_PIXEL_COUNT)
        # Intensity = np.flip(Intensity, 0)

        Metadata = self.WINDOW_LIMITS, self.Legend
        
        return Intensity, Metadata


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

        return self.PLANCK_SI * f / self.BOLTZMANN_SI / np.log(1 + 2 * self.PLANCK_SI * f**3 / self.C_LIGHT_SI**2 / I_nu)