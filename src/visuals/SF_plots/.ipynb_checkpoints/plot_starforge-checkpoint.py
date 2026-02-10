# plot_starforge.py

# This module generates diagnostic visualizations for raw STARFORGE simulation data. It leverages
# the project's internal I/O and transformation pipeline to ensure that the data visualized here 
# matches exactly what is passed to the SKIRT runner.

# Usage:
# Run from the project root (src/) as a module:
# $ python -m visuals.SF_plots.plot_starforge

# imports
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from file_io.loader import load_snapshot, get_header_data, get_particle_data, filter_by_id
from file_io.transform import finalize_dataset
from utils.constants import MSUN_PC2_TO_G_CM2

# -------------------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------------------
def coord_to_pixel(coords, bounds, resolution):
    """
    """
    xmin, xmax, ymin, ymax = bounds

    # normalizing to 0-1
    x_norm = (coords[:, 0] - xmin) / (xmax - xmin)
    y_norm = (coords[:, 1] - ymin) / (ymax - ymin)

    # scaling to reolution
    x_pix = x_norm * resolution
    y_pix = y_norm * resolution

    return np.column_stack([x_pix, y_pix])

def get_projection_hist(coords, weights, bounds, resolution):
    """
    """
    xmin, xmax, ymin, ymax = bounds
    H, _, _ = np.histogram2d(
        coords[:, 0], coords[:, 1],
        bins=resolution,
        range=[[xmin, xmax], [ymin, ymax]],
        weights=weights
    )
    return H.T
# -------------------------------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------------------------------
# gas surface density map
def plot_gas_surface_density(gas_data, bounds, resolution=800, save_path=None):
    """
    Creates a logarithmic column density map of the gas.

    Parameters
    ----------
    gas_data: dict
        Dictionary containing 'Coordinates' (N,3) and 'Masses' (N,).
    bounds: tuple
        (xmin, xmax, ymin, ymax, zmin, zmax) limits of the box in parsecs.
    resolution: int, optional
        Pixel resolution of the 2D grid, Default is 800.
    save_path: str, optional
        If provided, saves the figure to this path.
    """
    print("Generating Gas Surface Density Map")

    # extracting the data
    coords = gas_data['Coordinates']
    mass = gas_data['Masses']
    xmin, xmax, ymin, ymax, _, _ = bounds

    # calculating surface density
    H = get_projection_hist(coords, mass, (xmin, xmax, ymin, ymax), resolution)

    # converting the mass histogram to column density (sigma)
    pixel_area = ((xmax - xmin) / resolution) * ((ymax - ymin) / resolution)
    sigma_gas = H / pixel_area
    sigma_gas_cgs = sigma_gas * MSUN_PC2_TO_G_CM2

    extent = [0, resolution, 0, resolution]

    # plotting
    fig, ax = plt.subplots(figsize=(8,8))
    im = ax.imshow(sigma_gas_cgs, origin='lower', extent=extent, cmap='inferno', norm=LogNorm(vmin=1e-4, vmax=np.max(sigma_gas_cgs)))
    cbar = plt.colorbar(im, ax=ax, label=r'$\Sigma _{gas}$ [$g \ cm^{-2}$]')
    ax.set_xlabel('X [pixels]')
    ax.set_ylabel('Y [pixels]')
    ax.set_title(f'Gas Surface Density ({resolution} x {resolution})')

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()
    plt.close()

# overplotting stars
def plot_stars_on_gas(gas_data, sink_data, bounds, resolution=800, save_path=None):
    """
    """
    print("Generating Stars-on-Gas Composite")

    coords_gas = gas_data['Coordinates']
    mass_gas = gas_data['Masses']

    coords_stars = sink_data['Coordinates']
    
    xmin, xmax, ymin, ymax, _, _ = bounds
    
    # creating background gas map
    H = get_projection_hist(coords_gas, mass_gas, (xmin, xmax, ymin, ymax), resolution)
    pixel_area_pc2 = ((xmax - xmin) / resolution) * ((ymax - ymin) / resolution)
    sigma_gas_cgs = (H / pixel_area_pc2) * MSUN_PC2_TO_G_CM2
    
    # transofrming star coords to pixels
    pixel_stars = coord_to_pixel(coords_stars, (xmin, xmax, ymin, ymax), resolution)
    
    # plotting
    fig, ax = plt.subplots(figsize=(8, 8))
    extent = [0, resolution, 0, resolution]
    
    # background
    im = ax.imshow(sigma_gas_cgs, origin='lower', extent=extent, cmap='gray', 
                   norm=LogNorm(vmin=1e-3, vmax=np.max(sigma_gas_cgs)), alpha=0.9)
    
    # stars
    ax.scatter(pixel_stars[:, 0], pixel_stars[:, 1], s=15, c='cyan', marker='*', 
               edgecolors='none', alpha=0.8, label='Stars')
    
    ax.legend()
    ax.set_xlabel('X [pixels]', fontsize=12)
    ax.set_ylabel('Y [pixels]', fontsize=12)
    ax.set_title('Star Particles over Gas Density', fontsize=14)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()
    plt.close()

# temperature phase plot
def plot_mass_weighted_temperature_map(gas_data, bounds, resolution=800, save_path=None):
    """
    """
    print("Generating Mass-Weighted Temperature Map")
    
    coords = gas_data['Coordinates']
    mass = gas_data['Masses']
    temp = gas_data['Temperature']
    xmin, xmax, ymin, ymax, _, _ = bounds
    bound_2d = (xmin, xmax, ymin, ymax)
    
    # calculating numerator
    weighted_temp_H = get_projection_hist(coords, temp * mass, bound_2d, resolution)
    
    # calculating denominator
    mass_H = get_projection_hist(coords, mass, bound_2d, resolution)
    
    # dividing
    with np.errstate(invalid='ignore'):
        avg_temp_map = weighted_temp_H / mass_H
        
    # masking empty pixels so they show up as background color
    avg_temp_map = np.ma.masked_where(mass_H == 0, avg_temp_map)
    
    # plotting
    fig, ax = plt.subplots(figsize=(8, 8))
    extent = [0, resolution, 0, resolution]
    
    im = ax.imshow(avg_temp_map, origin='lower', extent=extent, cmap='RdYlBu_r', 
                   norm=LogNorm(vmin=10, vmax=2e2))
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Average Temperature [K]', fontsize=12)
    
    ax.set_xlabel('X [pixels]', fontsize=12)
    ax.set_ylabel('Y [pixels]', fontsize=12)
    ax.set_title('Projected Temperature (Mass Weighted)', fontsize=14)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()
    plt.close()

# density vs temperature diagram
def plot_phase_diagram(gas_data, save_path=None):
    """
    """
    print("Generating Phase Diagram (Temperature-Density)")

    temp = gas_data['Temperature']
    rho = gas_data['Density'] 
    mass = gas_data['Masses']

    fig, ax = plt.subplots(figsize=(8, 8))
    
    # 2d hist
    hist = ax.hist2d(rho, temp, bins=150, norm=LogNorm(), weights=mass, cmap='magma')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.set_xlabel(r'Density', fontsize=12) 
    ax.set_ylabel(r'Temperature [K]', fontsize=12)
    ax.set_title('Phase Diagram (Temperature vs Density)', fontsize=14)
    
    cbar = plt.colorbar(hist[3], ax=ax)
    cbar.set_label('Mass in Bin [Msun]', fontsize=12)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()
    plt.close()

# -------------------------------------------------------------------------------------------
# running
# -------------------------------------------------------------------------------------------
def run_diagnostics(snapshot_path, percentage=1, resolution=800):
    """
    Driver function to load data and run all plotting routines.
    """
    # output directory handling
    data_dir = os.path.dirname(snapshot_path)
    plot_dir = os.path.join(data_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    base_name = os.path.basename(snapshot_path).replace('.hdf5', '')

    # loads and transforms data
    print(f"Loading snapshot: {snapshot_path}")

    if not os.path.exists(snapshot_path):
        print(f"Error: File {snapshot_path} not found.")
        return

    with load_snapshot(snapshot_path) as f:
        header = get_header_data(f)
        raw_gas = get_particle_data(f, 0)
        raw_sinks = get_particle_data(f, 5)

    # filtering and centering
    gas_clean = filter_by_id(raw_gas)
    pt0, pt5, center, basis = finalize_dataset(header, gas_clean, raw_sinks, percentage)

    coords = pt0['Coordinates']
    bounds = (coords[:, 0]. min(), coords[:, 0].max(),
              coords[:, 1]. min(), coords[:, 1].max(),
              coords[:, 2]. min(), coords[:, 2].max())

    # making plots
    # # Gas Surface Density
    # plot_gas_surface_density(
    #     pt0, bounds, resolution=resolution,
    #     save_path=os.path.join(plot_dir, f"{base_name}_gas_density.png")
    # )
    
    # # Stars over Gas
    # plot_stars_on_gas(
    #     pt0, pt5, bounds, resolution=resolution,
    #     save_path=os.path.join(plot_dir, f"{base_name}_stars_on_gas.png")
    # )

    # Mass-Weighted Temperature Map (Spatial)
    plot_mass_weighted_temperature_map(
        pt0, bounds, resolution=resolution, 
        save_path=os.path.join(plot_dir, f"{base_name}_temp_map.png")
    )

    # # Phase Diagram (Data Space)
    # plot_phase_diagram(
    #     pt0, 
    #     save_path=os.path.join(plot_dir, f"{base_name}_phase_diagram.png")
    # )

if __name__ == "__main__":
    snapshot_file = '../data/snapshot_150.hdf5'
    run_diagnostics(snapshot_file, percentage=1)