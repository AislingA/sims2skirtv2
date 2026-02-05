# seds.py

# 

# imports
import numpy as np
import matplotlib.pyplot as plt
from utils.constants import JY_TO_CGS

def plot_sed(file_path, instrument):
    print(f"Plotting SED from: {file_path}")
    
    try:
        data = np.loadtxt(file_path)
    except OSError:
        print(f"Error: File {file_path} not found.")
        return

    wavelength = data[:, 0]  # lambda (microns)
    total_flux = data[:, 1] * JY_TO_CGS 

    flux_components = {
        'Transparent Flux': data[:, 2] * JY_TO_CGS,
        'Direct Primary Flux': data[:, 3] * JY_TO_CGS,
        'Scattered Primary Flux': data[:, 4] * JY_TO_CGS,
        'Direct Secondary Flux': data[:, 5] * JY_TO_CGS,
        'Scattered Secondary Flux': data[:, 6] * JY_TO_CGS,
        'Transparent Secondary Flux': data[:, 7] * JY_TO_CGS
    }


    for name, flux in flux_components.items():
        plt.figure(figsize=(10, 7))
        
        # plotting Total Flux as reference
        plt.plot(
            wavelength,
            total_flux,
            label='Total Flux',
            color='black',
            linewidth=2,
            zorder=10
        )
        
        # plotting each component
        plt.plot(
            wavelength,
            flux, 
            label=name,
            linestyle='--',
            linewidth=1.5
        )

        plt.xscale('log')
        plt.yscale('log')
        
        # limits based on the data
        valid_data = total_flux[total_flux > 0]
        if len(valid_data) > 0:
            ymax = np.max(valid_data) * 5
            ymin = np.min(valid_data) / 5
            plt.ylim(max(ymin, 1e-30), ymax)
            
        plt.title(f'SED Component: {name}')
        plt.xlabel('Wavelength [microns]')
        plt.ylabel(r'Flux Density $F_{\nu}$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]')
        plt.legend()
        plt.grid(True, which="major", alpha=0.5)
        plt.grid(True, which="minor", alpha=0.2)
        
        # Save figure
        save_name = f"sed_{instrument}_{name.replace(' ', '_').lower()}.png"
        plt.savefig(save_name, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_name}")
        plt.close()

if __name__ == "__main__":
    # spec instrument
    # plot_sed('../src/snapshot_150_spec_sed.dat', 'spec')
    # # jwst instrument
    # plot_sed('../src/snapshot_150_jwst_sed.dat', 'jwst')
    # # alma instrument
    # plot_sed('../src/snapshot_150_alma_sed.dat', 'alma')
    # # rubin instrument
    # plot_sed('../src/snapshot_150_rubin_sed.dat', 'rubin')
    # # spitzer instrument
    # plot_sed('../src/snapshot_150_spitzer_sed.dat', 'spitzer')
    plot_sed('../src/snapshot_150_panchromatic_sed.dat', 'pan')