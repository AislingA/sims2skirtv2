# projections.py

# python -m visuals.projections

# imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
import os
from utils.constants import MSUN_AU2_TO_G_CM2, AU_TO_PC

def plot_projection(fits_path, save_path=None, cmap='inferno', vmin=None, vmax=None, unit_conv=1, spatial_conv=1, label=None, use_pixels=False):
    """
    Plots a 2D projection from a SKIRT FITS file
    
    Parameters
    ----------
    fits_path : str
        Path to the SKIRT output FITS file.
    save_path : str, optional
        If provided, saves the figure to this path.
    cmap : str
        Colormap to match SF_plots (default: 'inferno').
    vmin, vmax : float, optional
        Min and max values for normalization.
    unit_conversion : float
        Factor to multiply the FITS data (e.g. Msun/AU2) to match desired units (g/cm2).
    spatial_scale_factor : float
        Factor to convert FITS spatial units (e.g. AU) to plot units (e.g. pc).
        For AU -> pc, use (AU_TO_CM / PC_TO_CM).
    label : str
        Label for the colorbar.
    """
    with fits.open(fits_path) as hdul:
        data = hdul[0].data
        header = hdul[0].header

        nx = header['NAXIS1']
        ny = header['NAXIS2']

        bunit = header['BUNIT']
        cunit1 = header['CUNIT1']
        cunit2 = header['CUNIT2']
        
        print(f"Physical unit of the array values: {bunit}")
        print(f"Physical units of the X-axis: {cunit1}")
        print(f"Physical units of the Y-axis: {cunit2}")

        if use_pixels:
            extent = [0, nx, 0, ny]
            xlabel = 'X [pixels]'
            ylabel = 'Y [pixels]'
        else:
            dx = header['CDELT1']
            dy = header['CDELT2']
            crval1 = header['CRVAL1']
            crval2 = header['CRVAL2']
            crpix1 = header['CRPIX1']
            crpix2 = header['CRPIX2']
            
            xmin_raw = crval1 - (crpix1 - 0.5) * dx
            xmax_raw = crval1 + (nx - crpix1 + 0.5) * dx
            ymin_raw = crval2 - (crpix2 - 0.5) * dy
            ymax_raw = crval2 + (ny - crpix2 + 0.5) * dy
            
            extent = [
                xmin_raw * spatial_scale_factor,
                xmax_raw * spatial_scale_factor,
                ymin_raw * spatial_scale_factor,
                ymax_raw * spatial_scale_factor
            ]
            xlabel = 'X [pc]'
            ylabel = 'Y [pc]'


    data_converted = data * unit_conv

    if vmin is None:
        vmin = np.percentile(data_converted[data_converted > 0], 1) if np.any(data_converted > 0) else 1e-4
    if vmax is None:
        vmax = np.max(data_converted)

    # plotting
    fig, ax = plt.subplots(figsize=(8,8))
    im = ax.imshow(data_converted, origin='lower', extent=extent, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax))
    cbar = plt.colorbar(im, ax=ax, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f'SKIRT Projection: {os.path.basename(fits_path)}')

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    else:
        plt.show()
    plt.close()

if __name__ == "__main__":
    # # imported medium density probe
    # print("Plotting imported medium density probe:")
    # fits_path = '../src/snapshot_150_mediumdensityprobe_0_Sigma.fits'
    # plot_projection(
    #     fits_path,
    #     unit_conv=MSUN_AU2_TO_G_CM2,
    #     spatial_conv=AU_TO_PC,
    #     vmin=1e-4,
    #     label=r'$\Sigma _{dust}$ [$g \ cm^{-2}$]',
    #     save_path="immeddensity.png",
    #     use_pixels=True
    # )

    # # dust density probe
    # print("Plotting dust density probe:")
    # fits_path = '../src/snapshot_150_densityprobe_dust_Sigma.fits'
    # plot_projection(
    #     fits_path,
    #     unit_conv=MSUN_AU2_TO_G_CM2,
    #     spatial_conv=AU_TO_PC,
    #     vmin=1e-4,
    #     label=r'$\Sigma _{dust}$ [$g \ cm^{-2}$]',
    #     save_path="dustdensity.png",
    #     use_pixels=True
    # )

    # # imported medium temperature probe
    # print("Plotting imported medium temperature probe:")
    # fits_path = '../src/snapshot_150_mediumtemperatureprobe_0_T.fits'
    # plot_projection(
    #     fits_path,
    #     unit_conv=1,
    #     spatial_conv=AU_TO_PC,
    #     vmin=None,
    #     label=r'$Temperature _{dust}$ [$K$]',
    #     save_path="immedtemp.png",
    #     use_pixels=True
    # )

    # radiation field probe
    print("Plotting radiation field probe:")
    fits_path = '../src/snapshot_150_radiationfieldprobe_J.fits'
    plot_projection(
        fits_path,
        unit_conv=1,
        spatial_conv=1,
        vmin=None,
        label=r'$Radiation field$ [$K$]',
        save_path="radfield.png",
        use_pixels=True
    )