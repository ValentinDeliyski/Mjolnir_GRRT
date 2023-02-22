import numpy as np
from matplotlib import pyplot as plt

from Support_functions.Parsers import*
from Support_functions.Shadows import*

Sim_Parser_0   = Simulation_Parser("JNW_Data0")
Sim_Parser_1   = Simulation_Parser("JNW_Data1")
Sim_Parser_2   = Simulation_Parser("JNW_Data2")
Sim_Parser_3   = Simulation_Parser("JNW_Data3")

Intensity_0, NT_Flux_0, NT_Redshift_0, NT_Flux_Shifted_0, Metadata_0 = Sim_Parser_0.get_plottable_sim_data()
Intensity_1, NT_Flux_1, NT_Redshift_1, NT_Flux_Shifted_1, Metadata_1 = Sim_Parser_1.get_plottable_sim_data()
Intensity_2, NT_Flux_2, NT_Redshift_2, NT_Flux_Shifted_2, Metadata_2 = Sim_Parser_2.get_plottable_sim_data()
Intensity_3, NT_Flux_3, NT_Redshift_3, NT_Flux_Shifted_3, Metadata_3 = Sim_Parser_3.get_plottable_sim_data()

# Sim_Parser.export_ehtim_data(Intensity)
axes_limits       = np.array([(limit) for limit in Metadata_0[2]])

Main_Figure = plt.figure()
Subplot = Main_Figure.add_subplot(111)

Data_to_plot = Intensity_1
Subplot.imshow(Data_to_plot, interpolation = 'bilinear', cmap = 'hot', extent = axes_limits, vmin = 0, vmax = np.max(Data_to_plot))

# figure = plt.figure()

# subfig_sim = figure.add_subplot(1, 2, 1)
# imgplot = plt.imshow(Intensity, interpolation = 'bilinear', cmap = 'hot', extent = axes_limits)
# subfig_sim.set_title('Ray Tracer')
# subfig_sim.set_xlabel(r'$\alpha$ [rad]')
# subfig_sim.set_ylabel(r'$\beta$ [rad]')

# subfig_eht = figure.add_subplot(1, 2, 2)
# imgplot = plt.imshow(ehtim_Intensity, interpolation = 'bilinear', cmap = 'hot', extent = ehtim_axes_limits*1e6)
# subfig_eht.set_title('Eht Imager')
# subfig_eht.set_xlabel(r'$\alpha$ [$\mu$as]')
# subfig_eht.set_ylabel(r'$\beta$  [$\mu$as]')

# # add_Kerr_Shadow(0.998, obs_distance = Metadata[0], obs_inclanation = np.deg2rad(Metadata[1]))
# add_Wormhole_Shadow(spin = -0.98, alpha = 2, obs_distance = Metadata_0[0], obs_inclanation = np.deg2rad(Metadata_0[1]), figure = Subplot)

# #---- Figure Labels -----#

print(np.max(Intensity_0), np.max(Intensity_1), np.max(Intensity_2))

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
plt.show()

"""
    Image cases:

        Case 1: 
            @ Thin disk - DISK_HEIGHT_SCALE = 100. / 3
            @ Radial Scale = 10 [M]
            @ Inclination = 80 [deg]
            @ Obs Window - X: [-30 30] [M], Y: [-30 30] [M]
            @ Absorbtion Coeff = 0
            @ Spin = 0.9
            @ Wormhole Alpha = 2

        Case 2: 
            @ Thick disk - DISK_HEIGHT_SCALE = 10. / 3
            @ Radial Scale = 10 [M]
            @ Inclination = 80 [deg]
            @ Obs Window - X: [-30 30] [M], Y: [-30 30] [M]
            @ Absorbtion Coeff = 0
            @ Spin = 0.9
            @ Wormhole Alpha = 2

        Case 3: 
            @ Thick disk - DISK_HEIGHT_SCALE = 10. / 3
            @ Radial Scale = 10 [M]
            @ Inclination = 20 [deg]
            @ Obs Window - X: [-30 30] [M], Y: [-30 30] [M]
            @ Absorbtion Coeff = 0
            @ Spin = 0.9
            @ Wormhole Alpha = 2

        Case 4: 
            @ Thick disk - DISK_HEIGHT_SCALE = 10. / 3
            @ Radial Scale = 5 [M]
            @ Inclination = 20 [deg]
            @ Obs Window - X: [-30 30] [M], Y: [-30 30] [M]
            @ Absorbtion Coeff = 10e5
            @ Spin = 0.9
            @ Wormhole Alpha = 2

        Case 5: 
            @ Thick disk - DISK_HEIGHT_SCALE = 10. / 3
            @ Radial Scale = 5 [M]
            @ Inclination = 20 [deg]
            @ Obs Window - X: [-30 30] [M], Y: [-30 30] [M]
            @ Absorbtion Coeff = 1e3
            @ Spin = 0.9
            @ Wormhole Alpha = 2
            
"""
