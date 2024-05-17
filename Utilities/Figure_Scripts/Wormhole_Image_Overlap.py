import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import csv

def parse_file(file_path: str):

    with open(file_path, 'r') as file:

        csvreader = csv.reader(file, delimiter = " ")

        r_image = []

        for _, row in enumerate(csvreader):

            r_image.append(float(row[0]))

    return r_image

if __name__  == "__main__":

    params = {"ytick.color" : "black",
              "xtick.color" : "black",
              "axes.labelcolor" : "black",
              "axes.edgecolor" : "black",
              "text.usetex" : True,
              "font.family" : "serif",
              "font.serif" : ["Computer Modern Serif"]}
    plt.rcParams.update(params)

    r_image_in_SCHW  = parse_file("Scwarzchild_r_4.5_20.0_deg_indirect.csv")
    r_image_out_SCHW = parse_file("Scwarzchild_r_500_20.0_deg_indirect.csv")
    r_image_in_WH    = []
    r_image_out_WH   = []

    Figure_Overlap, (Inner_plot) = plt.subplots(1, 1, gridspec_kw = {'width_ratios': [1]}, constrained_layout = True)

    for gamma in np.linspace(0.52, 1, 12):

        r_image_in_WH.append(parse_file("JNW_r_4.5_gamma_{}_20.0_deg_indirect.csv".format(round(gamma, 2))))
        r_image_out_WH.append(parse_file("JNW_r_500_gamma_{}_20.0_deg_indirect.csv".format(round(gamma, 2))))
    
    x_axis = np.linspace(0, 2, 500) 

    Fontsize = 32

    # =================================== Colormap =================================== #

    cmap = matplotlib.colormaps['plasma']
    
    colorbar_map = matplotlib.cm.ScalarMappable(cmap = matplotlib.colormaps['plasma'])
    colorbar_map.set_clim([0.52,1])

    Colorbar = Figure_Overlap.colorbar(colorbar_map, ax = Inner_plot)
    Colorbar.set_label(r"$\gamma$", fontsize = Fontsize)
    Colorbar.ax.tick_params(labelsize = Fontsize)     

    # ================================================================================ #

    Inner_plot.plot(x_axis, r_image_out_SCHW, "k--")
    Inner_plot.plot(x_axis, r_image_in_SCHW, "k--")
       
    # Outer_Plot.plot(x_axis, r_image_out_SCHW, "k--")
    # Outer_Plot.plot(x_axis, r_image_in_SCHW, "k--")


    for gamma_value in range(12):

        if gamma_value % 1 == 0:

            Color_range = cmap(gamma_value / 11) 

            Inner_plot.plot(x_axis, r_image_in_WH[gamma_value], color = Color_range)
            Inner_plot.plot(x_axis, r_image_out_WH[gamma_value], color = Color_range, label = r"$\gamma = ${}".format(round(gamma_value / 11 + (1 - gamma_value / 11) * 0.52, 2)))

            # Inner_plot.fill_between(x_axis, y1 = r_image_in_WH[gamma_value], y2 = r_image_out_WH[gamma_value], where = x_axis, alpha = 0.5, color = Color_range)

    Inner_plot.set_ylim([4.5, 7])
    Inner_plot.set_ylabel(r"$\xi$ [M]", fontsize = Fontsize)
    Inner_plot.set_xlabel(r"$\phi$ [$\pi$]", fontsize = Fontsize)
    Inner_plot.set_xlim([0, 2])
    Inner_plot.tick_params(axis='x', labelsize=Fontsize)
    Inner_plot.tick_params(axis='y', labelsize=Fontsize)

    # Inner_plot.legend(loc = "upper left", fontsize = 22)

    plt.show()