import csv
import numpy as np
from itertools import islice

class Simulation_Parser():

    def __init__(self, File_name: str) -> None:

        with open("..\\Sim_Results\\" + File_name + ".txt", 'r') as file:

            self.HEADER_ROW_COUNT = 5

            csvreader = csv.reader(file, delimiter = ":")
            row_count = sum(1 for row in file) - self.HEADER_ROW_COUNT
            file.seek(0)

            self.OBS_DISTANCE    = float(csvreader.__next__()[1])
            self.OBS_INCLICATION = float(csvreader.__next__()[1])

            self.WINDOW_LIMITS = [float(limit) for limit in csvreader.__next__()[1].split(',')]
            
            Resolution_list = csvreader.__next__()[1].split(' ')

            self.X_PIXEL_COUNT   = int(Resolution_list[1])
            self.Y_PIXEL_COUNT   = int(Resolution_list[3])

            self.Legend = csvreader.__next__()

            csvreader = csv.reader(file, delimiter = " ")

            self.X_coords        = np.zeros(row_count)
            self.Y_coords        = np.zeros(row_count)
            self.NT_Flux         = np.zeros(row_count)
            self.Intensity       = np.zeros(row_count)
            self.NT_Redshift     = np.zeros(row_count)
            self.NT_Flux_Shifted = np.zeros(row_count)

            index = 0

            for row in csvreader:

                self.X_coords[index] = row[0]
                self.Y_coords[index] = row[1]

                self.NT_Redshift[index]  = row[2]
                self.NT_Flux[index]      = row[3]
                self.Intensity[index]    = row[4]

                self.NT_Flux_Shifted[index] = self.NT_Redshift[index]**4*self.NT_Flux[index]

                index += 1

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

        ehtim_x_fov = 2 * 7.500000e-05
        ehtim_y_fov = 2 * 7.500000e-05

        R_M87_LY = 53490000
        LY_TO_M  = 9.461e+15
        R_M87_M  = R_M87_LY * LY_TO_M

        M_SUN_KG = 2e30
        M_M87_KG = 6.5e9 * M_SUN_KG

        G = 6.67e-11
        c = 3e8

        R_SCH_M87 = M_M87_KG * G / c**2

        R_OBJ_DIMENTIONLESS = R_M87_M / R_SCH_M87

        RAD_TO_ARCS     = 206265
        CARTESIAN_TO_AS = RAD_TO_ARCS / R_OBJ_DIMENTIONLESS
        
        dx = (max(self.X_coords) - min(self.X_coords)) * CARTESIAN_TO_AS
        dy = (max(self.Y_coords) - min(self.Y_coords)) * CARTESIAN_TO_AS

        stupid_eht_scaling = 1e9

        formatted_sim_data = data.reshape(1, self.X_PIXEL_COUNT * self.Y_PIXEL_COUNT).flatten()

        array_to_export = np.array([self.X_coords / max(self.X_coords) * ehtim_x_fov / 2, 
                                    self.Y_coords / max(self.Y_coords) * ehtim_y_fov / 2, 
                                    formatted_sim_data * dx * dy / R_OBJ_DIMENTIONLESS**2 * stupid_eht_scaling]).T

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
