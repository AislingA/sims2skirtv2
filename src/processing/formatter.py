# formatter.py

# To format the processed particle data into plain-text input files 
# compatible with SKIRT. This involves calculating secondary physical 
# properties (like stellar effective temperature), unit conversions (like 
# Solar Radii to km), and generating the coordinate bounds required for 
# the SKIRT spatial grid configuration.

# imports
import numpy as np
import os

def compute_stellar_temperature(luminosity, star_radius):
    """
    Computes the surface effective temperature of a star using the 
    Stefan-Boltzmann Law.

    Parameters
    ----------
    luminosity : np.ndarray
        Stellar luminosity in units of Lsun.
    star_radius : np.ndarray
        Stellar radius in units of Rsun.

    Returns
    -------
    np.ndarray
        Effective temperature in Kelvin (K).
    """
    T_sun = 5777 # K
    return T_sun * (luminosity / star_radius**2)**0.25

def print_stats(name, data_dict, verbose=False):
    """
    Utility function to print the dynamic range of datasets for verification.

    Parameters
    ----------
    name : str
        Identifier for the data group ("Source" or "Medium").
    data_dict : dict
        The dictionary of arrays to summarize.
    verbose : bool, optional
        If False, the function returns without printing. Default is False.
    """
    if not verbose:
        return
    print(f'--- {name} Data min and max values ---')
    for key, val in data_dict.items():
        if isinstance(val, np.ndarray) and val.ndim > 0:
            print(f'{key}: Min {np.min(val):.2e}, Max {np.max(val):.2e}')

def format_source_file(pt5_data, output_path, verbose=False):
    """
    Generates a SKIRT-compatible text file for stellar sources (PartType5).

    Parameters
    ----------
    pt5_data : dict
        Cleaned and transformed sink/star particle data.
    output_path : str
        The file path where the .txt source file will be saved.
    verbose : bool, optional
        Whether to print statistical summaries of the exported data.

    Returns
    -------
    str
        The path to the created source input file.
    """
    # Calculating T_eff
    temp = compute_stellar_temperature(pt5_data['StarLuminosity_Solar'], pt5_data['ProtoStellarRadius_inSolar'])

    # Structuring for SKIRT: x, y, z, h, R, T
    radius_km = pt5_data['ProtoStellarRadius_inSolar'] * 695700

    export_data = np.column_stack([
        pt5_data['Coordinates'][:, 0],
        pt5_data['Coordinates'][:, 1],
        pt5_data['Coordinates'][:, 2],
        pt5_data['BH_AccretionLength'],
        radius_km,
        temp
    ])

    print_stats("Source", pt5_data, verbose)

    header = "# x(pc) y(pc) z(pc) h(pc) R(km) T(K)"
    np.savetxt(output_path, export_data, fmt='%.6e', header=header, comments='')
    return output_path

def format_gas_file(pt0_data, output_path, verbose=False):
    """
    Generates a SKIRT-compatible text file for the gas medium (PartType0).

    Parameters
    ----------
    pt0_data : dict
        Cleaned and transformed gas particle data.
    output_path : str
        The file path where the .txt gas file will be saved.
    verbose : bool, optional
        Whether to print statistical summaries of the exported data.

    Returns
    -------
    tuple
        (str, float, float, float, float, float, float)
        The output path followed by the (xmin, xmax, ymin, ymax, zmin, zmax) 
        bounds of the gas distribution in parsecs.
    """
    # Structuring for SKIRT: x, y, z, h, M, T
    export_data = np.column_stack([
        pt0_data['Coordinates'][:, 0],
        pt0_data['Coordinates'][:, 1],
        pt0_data['Coordinates'][:, 2],
        pt0_data['SmoothingLength'],
        pt0_data['Masses'],
        pt0_data['Temperature']
    ])

    print_stats("Medium", pt0_data, verbose)

    # Calculating bounds for the SKIRT geometry setup
    coords = pt0_data['Coordinates']
    bounds = (coords[:, 0]. min(), coords[:, 0].max(),
              coords[:, 1]. min(), coords[:, 1].max(),
              coords[:, 2]. min(), coords[:, 2].max())

    header = "# x(pc) y(pc) z(pc) h(pc) M(Msun) T(K)"
    np.savetxt(output_path, export_data, fmt='%.6e', header=header, comments='')
    return (output_path, *bounds)