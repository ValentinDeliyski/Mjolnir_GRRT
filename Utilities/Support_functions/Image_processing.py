from numpy import sin, cos, sqrt, exp, arctan2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cl
from Support_functions.Parsers import*

def generate_general_gaussian_template(N_pixels: int, template_params: dict, FOV: float) -> np.ndarray:

    elipse_x_offset = template_params["x0"] + FOV / 2
    elipse_y_offset = template_params["y0"] + FOV / 2
    rot_angle       = template_params["rot_angle"]
    template_std    = template_params["sigma"]
    d0              = template_params["d0"]
    tau             = template_params["tau"]
    slash           = template_params["slash"]
    slash_angle     = template_params["slash_angle"]

    a_elipse = d0 / 2 / sqrt(1 - tau)
    b_elipse = d0 / 2 * sqrt(1 - tau)

    elipse_samples = np.linspace(0, 2 * np.pi, 1000)
    x_elipse = a_elipse * cos(elipse_samples) * cos(rot_angle) - b_elipse * sin(elipse_samples) * sin(rot_angle) + elipse_x_offset 
    y_elipse = b_elipse * sin(elipse_samples) * cos(rot_angle) + a_elipse * cos(elipse_samples) * sin(rot_angle) + elipse_y_offset 
 
    template_image = np.zeros((N_pixels, N_pixels))

    for px in range(N_pixels):

        for py in range(N_pixels):

            x_pixel = px / (N_pixels - 1) * FOV
            y_pixel = py / (N_pixels - 1) * FOV

            d_squared = (x_pixel - x_elipse)**2 + (y_pixel - y_elipse)**2 

            azimuth_rel_to_center = arctan2(y_pixel - elipse_y_offset, x_pixel - elipse_x_offset)

            template_image[N_pixels - 1 - py, N_pixels - 1 - px] = (1 + slash * cos(azimuth_rel_to_center - slash_angle)) * exp(-np.min(d_squared) / 2 / template_std**2)
            
    template_image = template_image / np.max(template_image)
    
    return template_image

def get_template_pixel_mask(template_params, FOV, N_pixels):

    if template_params["Gaussian_2"] != None:

        elipse_x_offset =     (template_params["Gaussian_1"]["x0"] + template_params["Gaussian_2"]["x0"]) / 2 + FOV / 2
        elipse_y_offset =     (template_params["Gaussian_1"]["y0"] + template_params["Gaussian_2"]["y0"]) / 2 + FOV / 2
        rot_angle       =     (template_params["Gaussian_1"]["rot_angle"] + template_params["Gaussian_2"]["rot_angle"]) / 2 
        template_std    = sqrt(template_params["Gaussian_1"]["sigma"]**2 + template_params["Gaussian_2"]["sigma"]**2)
        d0              =     (template_params["Gaussian_1"]["d0"] + template_params["Gaussian_2"]["d0"]) / 2
        tau             =     (template_params["Gaussian_1"]["tau"] + template_params["Gaussian_2"]["tau"]) / 2

    else:

        elipse_x_offset = template_params["Gaussian_1"]["x0"] + FOV / 2
        elipse_y_offset = template_params["Gaussian_1"]["y0"] + FOV / 2
        rot_angle       = template_params["Gaussian_1"]["rot_angle"] 
        template_std    = template_params["Gaussian_1"]["sigma"]
        d0              = template_params["Gaussian_1"]["d0"]
        tau             = template_params["Gaussian_1"]["tau"]

    a_elipse = d0 / 2 / sqrt(1 - tau)
    b_elipse = d0 / 2 * sqrt(1 - tau)

    elipse_samples = np.linspace(0, 2 * np.pi, 1000)
    x_elipse = a_elipse * cos(elipse_samples) * cos(rot_angle) - b_elipse * sin(elipse_samples) * sin(rot_angle) + elipse_x_offset 
    y_elipse = b_elipse * sin(elipse_samples) * cos(rot_angle) + a_elipse * cos(elipse_samples) * sin(rot_angle) + elipse_y_offset

    ring_mask      = np.zeros((N_pixels, N_pixels))
    dark_spot_mask = np.zeros((N_pixels, N_pixels))

    for px in range(N_pixels):

        for py in range(N_pixels):

            x_pixel = px / (N_pixels - 1) * FOV
            y_pixel = py / (N_pixels - 1) * FOV

            d_squared = (x_pixel - x_elipse)**2 + (y_pixel - y_elipse)**2 

            ring_mask[N_pixels - 1 - py, N_pixels - 1 - px] = np.min(d_squared) < template_std**2

            x_term = (x_pixel - elipse_x_offset) * cos(rot_angle) + (y_pixel - elipse_y_offset) * sin(rot_angle)
            y_term = (x_pixel - elipse_x_offset) * sin(rot_angle) - (y_pixel - elipse_y_offset) * cos(rot_angle)

            point_inide_ellipse = (x_term / a_elipse)**2 + (y_term / b_elipse)**2 < 1
            dark_spot_mask[N_pixels - 1 - py, N_pixels - 1 - px] = (not ring_mask[N_pixels - 1 - py, N_pixels - 1 - px]) and point_inide_ellipse
            
    return ring_mask, dark_spot_mask

def get_brigness_depression_ratio(ring_mask, dark_spot_mask, Ehtim_intensity):

    return np.min(Ehtim_intensity[dark_spot_mask != 0]) / np.mean(Ehtim_intensity[ring_mask != 0])

def get_template_slices(N_pixels: int, template: np.array, template_params: dict, FOV: float) -> tuple:

    if (template_params["Gaussian_2"] != None):

        slice_x_offset = (template_params["Gaussian_1"]["x0"] + template_params["Gaussian_2"]["x0"]) / 2 + FOV / 2
        slice_y_offset = (template_params["Gaussian_1"]["y0"] + template_params["Gaussian_2"]["y0"]) / 2 + FOV / 2
    else:
        slice_x_offset = template_params["Gaussian_1"]["x0"] + FOV / 2
        slice_y_offset = template_params["Gaussian_1"]["y0"] + FOV / 2

    x_slice_y_index = N_pixels - int( slice_y_offset / FOV * (N_pixels - 1) )
    y_slice_x_index = N_pixels - int( slice_x_offset / FOV * (N_pixels - 1) )

    x_slice = template[x_slice_y_index, 0 : N_pixels] / np.max(template)
    y_slice = np.flip(template[0 : N_pixels, y_slice_x_index]) / np.max(template)

    return x_slice, slice_x_offset, y_slice, slice_y_offset

def plot_VIDA_templte(Ehtim_Parsers: list, VIDA_parser: VIDA_params_Parser, CROP: str, crop_rel_rage: list, Plot_Brihtness_T: bool) -> None:

    Units = Units_class()
    Ehtim_Parser = Ehtim_Parsers[0]

    axes_limits     = np.array([(limit) for limit in Ehtim_Parser.WINDOW_LIMITS ]) * Units.MEGA
    Ehtim_image_FOV = np.abs(axes_limits[0] - axes_limits[1])  # Units of [uas]

    Ehtim_image_res = Ehtim_Parser.X_PIXEL_COUNT
    pixel_size      = Ehtim_image_FOV / Ehtim_image_res / Units.MEGA * Units.ARCSEC_TO_RAD  

    Intensity_ehtim_jy = np.zeros((Ehtim_image_res, Ehtim_image_res))
    Intensity_ehtim_T = np.zeros((Ehtim_image_res, Ehtim_image_res))

    for Parser in Ehtim_Parsers:

        temp_Intensity_ehtim_jy, _ = Parser.get_plottable_ehtim_data()
        # Convert the Flux from [Jy] to brightness temperature in Giga [K]
        temp_Intensity_ehtim_T = Units.Spectral_density_to_T(temp_Intensity_ehtim_jy / pixel_size**2 / Units.W_M2_TO_JY, Parser.OBS_FREQUENCY * Units.GIGA) / Units.GIGA
        
        Intensity_ehtim_jy += temp_Intensity_ehtim_jy
        Intensity_ehtim_T  += temp_Intensity_ehtim_T

    template = generate_general_gaussian_template(Ehtim_image_res, VIDA_parser.template_params["Gaussian_1"], Ehtim_image_FOV)

    if VIDA_parser.template_params["Gaussian_2"] != None:
        
        template_2 = generate_general_gaussian_template(Ehtim_image_res, VIDA_parser.template_params["Gaussian_2"], Ehtim_image_FOV)
        template = template + template_2

    if Plot_Brihtness_T:
        Intensity_ehtim = Intensity_ehtim_T
        colorbar_legend = r"Brightness Temperature [$10^9$K]"

    else:
        Intensity_ehtim = Intensity_ehtim_jy * Units.KILO
        colorbar_legend = r"Flux Per Pixel [mJy]"

    template_x_slice, slice_x_offset, template_y_slice, slice_y_offset = get_template_slices(Ehtim_image_res, template, VIDA_parser.template_params, Ehtim_image_FOV)
    Ehtim_x_slice, _, Ehtim_y_slice, _ = get_template_slices(Ehtim_image_res, Intensity_ehtim, VIDA_parser.template_params, Ehtim_image_FOV)

    template_fig = plt.figure(figsize=(21,5))

    #========================= Ehtim Image Plot =========================#

    if CROP:

        x_crop_range = np.array([(-(slice_x_offset - Ehtim_image_FOV / 2) - crop_rel_rage[0]), (-(slice_x_offset - Ehtim_image_FOV / 2) + crop_rel_rage[1])])
        y_crop_range = np.array([(-(slice_y_offset - Ehtim_image_FOV / 2) - crop_rel_rage[2]), (-(slice_y_offset - Ehtim_image_FOV / 2) + crop_rel_rage[3])])
        
        # The desired crop window could "cut" outside the simulated window
        FOV_overshoot_x = max(np.absolute(x_crop_range)) - Ehtim_image_FOV / 2
        FOV_overshoot_y = max(np.absolute(y_crop_range)) - Ehtim_image_FOV / 2

        FOV_overshoot = max(FOV_overshoot_x, FOV_overshoot_y)
        
        axes_limits = [-crop_rel_rage[0], 
                        crop_rel_rage[1], 
                       -crop_rel_rage[2], 
                        crop_rel_rage[3]]

        # If it does, crop to the end of the simulated window, while keeping the aspec ratio
        if FOV_overshoot > 0:

            x_crop_range = x_crop_range - np.array([-FOV_overshoot, FOV_overshoot])
            y_crop_range = y_crop_range - np.array([-FOV_overshoot, FOV_overshoot])

            axes_limits = [-crop_rel_rage[0] + FOV_overshoot, 
                            crop_rel_rage[1] - FOV_overshoot, 
                           -crop_rel_rage[2] + FOV_overshoot, 
                            crop_rel_rage[3] - FOV_overshoot]

        x_crop_idx = (x_crop_range / (Ehtim_image_FOV / 2) + 1) / 2 * Ehtim_image_res
        y_crop_idx = (y_crop_range / (Ehtim_image_FOV / 2) + 1) / 2 * Ehtim_image_res

        x_crop_idx = x_crop_idx.astype(int) - 1
        y_crop_idx = y_crop_idx.astype(int) - 1

    else:
    
        x_crop_idx = [0, (Ehtim_image_res - 1)]
        y_crop_idx = [0, (Ehtim_image_res - 1)]

    crop_res = x_crop_idx[1] - x_crop_idx[0]

    # The literature (for some reason) has the X axis going positive to negative, 
    # so I invert the X axis limits

    axes_limits[0] = -axes_limits[0]
    axes_limits[1] = -axes_limits[1]

    x_coords = np.linspace(axes_limits[0], axes_limits[1], crop_res)
    y_coords = np.linspace(axes_limits[2], axes_limits[3], crop_res)

    Subplot = template_fig.add_subplot(141)
    Ehtim_crop        = Intensity_ehtim[y_crop_idx[0] : y_crop_idx[1], x_crop_idx[0] : x_crop_idx[1]]
    Ehtim_crop_figure = Subplot.imshow(Ehtim_crop, cmap = "hot", extent = axes_limits)

    if CROP:

        Subplot.plot(np.zeros(crop_res), y_coords, "r")
        Subplot.plot(x_coords, np.zeros(crop_res), "b")

    else:

        Subplot.plot((slice_x_offset - Ehtim_image_FOV / 2) * np.ones(crop_res), y_coords, "r")
        Subplot.plot(x_coords, (slice_y_offset - Ehtim_image_FOV / 2) * np.ones(crop_res), "b")

    Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]')
    Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]')

    if len(Ehtim_Parsers) > 1:
        Frequency_str = "{"

        for Freq_num, Parser in enumerate(Ehtim_Parsers):

            Frequency_str += str(int(Parser.OBS_FREQUENCY))

            if Freq_num < len(Ehtim_Parsers) - 1:

                Frequency_str += ", "

        Frequency_str += "}"

    else:
        Frequency_str = str(Ehtim_Parser.OBS_FREQUENCY)

    Subplot.set_title("EHTIM Output at {}GHz".format(Frequency_str))

    colorbar = template_fig.colorbar(Ehtim_crop_figure, ax = Subplot, fraction=0.046, pad=0.04)
    colorbar.set_label(colorbar_legend)

    #========================= Template Plot =========================#

    Subplot = template_fig.add_subplot(142)
    Subplot.imshow(template[y_crop_idx[0] : y_crop_idx[1], x_crop_idx[0] : x_crop_idx[1]], cmap = "hot", extent = axes_limits)

    if CROP:

        Subplot.plot(np.zeros(crop_res), y_coords, "r--")
        Subplot.plot(x_coords, np.zeros(crop_res), "b--")

    else:

        Subplot.plot((slice_x_offset - Ehtim_image_FOV / 2) * np.ones(crop_res), y_coords, "r--")
        Subplot.plot(x_coords, (slice_y_offset - Ehtim_image_FOV / 2) * np.ones(crop_res), "b--")

    Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]')
    Subplot.set_ylabel(r'$\delta_{rel}\,\,[\mu$as]')

    ring_mask, dark_spot_mask = get_template_pixel_mask(VIDA_parser.template_params, Ehtim_image_FOV, Ehtim_image_res)
    Subplot.set_title("VIDA Template")

    colorbar = template_fig.colorbar(Ehtim_crop_figure, ax = Subplot, fraction=0.046, pad=0.04)
    colorbar.set_label(colorbar_legend)

    #------------------------ Y Slice Plot ------------------------#

    Subplot = template_fig.add_subplot(143)
    Subplot.plot(y_coords, Ehtim_y_slice[Ehtim_image_res - y_crop_idx[1] : Ehtim_image_res - y_crop_idx[0]], "r")
    Subplot.plot(y_coords, template_y_slice[Ehtim_image_res - y_crop_idx[1] : Ehtim_image_res - y_crop_idx[0]], "r--")

    Subplot.set_ylim([0,1])
    Subplot.set_xlim([axes_limits[0], axes_limits[1]])

    Subplot.set_aspect(np.absolute(axes_limits[1] - axes_limits[0]))
    Subplot.tick_params(left = False, right = False, labelleft = False,
                        labelbottom = True, bottom = True)

    Subplot.set_xlabel(r'$\delta_{rel}\,\,[\mu$as]')
    Subplot.set_title("Y Intensity Slice")

    # Subplot.imshow(ring_mask, cmap = "hot", extent = axes_limits)

    #------------------------ X Slice Plot ------------------------#

    Subplot = template_fig.add_subplot(144)
    Subplot.plot(x_coords, Ehtim_x_slice[x_crop_idx[0] : x_crop_idx[1]], "b")
    Subplot.plot(x_coords, template_x_slice[x_crop_idx[0] : x_crop_idx[1]], "b--")

    Subplot.set_ylim([0,1])
    Subplot.set_xlim([axes_limits[0], axes_limits[1]])

    Subplot.set_aspect(np.absolute(axes_limits[1] - axes_limits[0]))
    Subplot.tick_params(left = False, right = False, labelleft = False,
                            labelbottom = True, bottom = True)
        
    Subplot.set_xlabel(r'$\alpha_{rel}\,\,[\mu$as]')
    Subplot.set_title("X Intensity Slice")

    template_fig.subplots_adjust(bottom=0.22)

    # Subplot.imshow(dark_spot_mask, cmap = "hot", extent = axes_limits)
    
    Brightness_ratio_str = "f at {}GHz = {}".format(Frequency_str, 
                                                    get_brigness_depression_ratio(ring_mask, dark_spot_mask, Intensity_ehtim_jy))
    
    return template_fig, Brightness_ratio_str

