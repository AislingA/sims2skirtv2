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

def plot_gas_surface_density(gas_data, bounds, resolution=800, save_path=None, use_pixels=False):
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
    use_pixels: bool
        If True, plots X and Y axes in pixels instead of physical units
    """
    print("Generating Gas Surface Density Map")

    # extracting the data
    coords = gas_data['Coordinates']
    mass = gas_data['Masses']

    xmin, xmax, ymin, ymax, _, _ = bounds

    # creating a 2D histogram weighted by mass (surface density)
    # projecting all particles along the z-axis
    H, x_edges, y_edges = np.histogram2d(
        coords[:, 0], coords[:, 1],
        bins=resolution,
        range=[[xmin, xmax], [ymin, ymax]],
        weights=mass
    )

    # converting the mass histogram to column density (sigma)
    pixel_area = ((xmax - xmin) / resolution) * ((ymax - ymin) / resolution)
    sigma_gas = H.T / pixel_area
    sigma_gas_cgs = sigma_gas * MSUN_PC2_TO_G_CM2

    # determining plot extent
    if use_pixels:
        extent = [0, resolution, 0, resolution]
        xlabel = 'X [pixels]'
        ylabel = 'Y [pixels]'
    else:
        extent=[xmin, xmax, ymin, ymax]
        xlabel = 'X [pc]'
        ylabel = 'Y [pc]'

    # plotting
    fig, ax = plt.subplots(figsize=(8,8))
    im = ax.imshow(sigma_gas_cgs, origin='lower', extent=extent, cmap='inferno', norm=LogNorm(vmin=1e-4, vmax=np.max(sigma_gas_cgs)))
    cbar = plt.colorbar(im, ax=ax, label=r'$\Sigma _{gas}$ [$G \ cm^{-2}$]')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Gas Surface Density')

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()
    plt.close()

def plot_temperature_slice(gas_data, bounds, slice_thickness=0.5, resolution=800, save_path=None):
    """
    Creates a slice of the gas temperature field (projected average of a thin slab).
    Useful for visualizing hot HII regions vs cold molecular gas.

    Parameters
    ----------
    gas_data : dict
        Dictionary containing 'Coordinates'
    bounds : tuple
        (xmin, xmax, ymin, ymax, zmin, zmax).
    slice_thickness : float
        Thickness of the slice in Z (pc). Particles outside |z| < thickness/2 are ignored.
    resolution : int
        Grid resolution.
    save_path : str, optional
        File path to save the plot.
    """
    print("Generating Temperature Slice")

    coords = gas_data['Coordinates']
    temp = gas_data['Temperature']
    mass = gas_data['Masses']

    z_min_slice = -slice_thickness/2
    z_max_slice = slice_thickness/2

    mask = (coords[:, 2] > z_min_slice) & (coords[:, 2] < z_max_slice)

    sliced_x = coords[mask, 0]
    sliced_y = coords[mask, 1]
    sliced_T = temp[mask]
    sliced_m = mass[mask]

    if len(sliced_T) == 0:
        print("Warning: No particles found in the slice.")
        return

    xmin, xmax, ymin, ymax, _, _ = bounds

    # weighted histogram for temperature
    H_num, _, _ = np.histogram2d(
        sliced_x, sliced_y, bins=resolution,
        range=[[xmin, xmax], [ymin, ymax]], weights=sliced_T * sliced_m
    )
    H_den, _, _ = np.histogram2d(
        sliced_x, sliced_y, bins=resolution,
        range=[[xmin, xmax], [ymin, ymax]], weights=sliced_m
    )
    
    # avoiding division by zero
    with np.errstate(invalid='ignore'):
        avg_temp = H_num.T / H_den.T
    
        # plotting
    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(avg_temp, origin='lower', extent=[xmin, xmax, ymin, ymax], cmap='coolwarm', norm=LogNorm(vmin=10, vmax=1e4)) # 10K to 10,000K range
    cbar = plt.colorbar(im, ax=ax, label='Temperature [K]')
    ax.set_xlabel('X [pc]')
    ax.set_ylabel('Y [pc]')
    ax.set_title(f'Temperature Slice (z +/- {slice_thickness/2} pc)')
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()
    plt.close()

def plot_phase_diagram(gas_data, save_path=None):
    """
    Creates a 2D Histogram Phase Diagram (Density vs Temperature).
    
    Parameters
    ----------
    gas_data : dict
        Must contain 'Density' and 'Temperature'.
    save_path : str, optional
        File path to save the plot.
    """
    print("Generating Phase Diagram...")

    temp = gas_data['Temperature']
    rho = gas_data['Density']     # Check units! (Usually code units, convert to g/cm^3 or n_H if desired)
    mass = gas_data['Masses']

    fig, ax = plt.subplots(figsize=(8, 6))
    
    # 2D Histogram: X=Density, Y=Temperature, Color=Mass
    hist = ax.hist2d(rho, temp, bins=100, norm=LogNorm(), weights=mass, cmap='magma')
    
    # Formatting
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Density [Code Units]') 
    ax.set_ylabel(r'Temperature [K]')
    ax.set_title('Phase Diagram (Mass Weighted)')
    
    cbar = plt.colorbar(hist[3], ax=ax, label='Mass in Bin')
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()
    plt.close()

def run_diagnostics(snapshot_path, percentage=1):
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
    plot_gas_surface_density(
        pt0, bounds,
        save_path=os.path.join(plot_dir, f"{base_name}_surface_density.png"),
        use_pixels=True
    )

    plot_temperature_slice(
            pt0, bounds, slice_thickness=1.0, 
            save_path=os.path.join(plot_dir, f"{base_name}_temp_slice.png")
        )

    plot_phase_diagram(
            pt0, 
            save_path=os.path.join(plot_dir, f"{base_name}_phase_diagram.png")
        )

if __name__ == "__main__":
    snapshot_file = '../data/snapshot_150.hdf5'
    run_diagnostics(snapshot_file, percentage=1)