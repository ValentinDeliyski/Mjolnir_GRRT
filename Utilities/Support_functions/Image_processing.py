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

def get_template_pixel_mask(template_params, FOV, N_pixels, std_scale: float = 1):

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

            ring_mask[N_pixels - 1 - py, N_pixels - 1 - px] = np.min(d_squared) < (std_scale * template_std)**2

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

