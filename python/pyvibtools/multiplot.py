import numpy as np
from typing import List
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from .calculator import vibtoolsCalculator

def process_file_multiplot(file_name, calculator_list):
    """
    Reads the file into a new vibtoolsCalculator object and appends it to the calculator_list.
    
    Parameters:
    - file_name: Name of the file to be read.
    - calculator_list: A list to which the new vibtoolsCalculator object will be appended.
    """
    # Create a new vibtoolsCalculator object
    vib_calc = vibtoolsCalculator()
    
    # Read the file into the vibtoolsCalculator object
    vib_calc.read(vibspecfile=file_name)
    
    # Append the object to the list
    calculator_list.append(vib_calc)


def multiplot(speclist: List['vibtoolsCalculator'], linewidth: int = 2, 
              figsize: tuple = (9, 6), save: str = None)
    """
    Plot the computed vibrational spectrum as a stick spectrum for each object in speclist.

    Parameters:
    - speclist: A list of vibtoolsCalculator objects.
    - linewidth: The width of the sticks (default: 2).
    - figsize: The size of the figure (default: (9, 6)).
    - save: A file name to which the spectrum is saved (needs an extension).
    """

    fig, ax = plt.subplots(figsize=figsize)

    # Get a color palette from husl based on the number of spectra
    num_spectra = len(speclist)
    colors = [mcolors.husl_to_rgb(h, 1, 0.65) for h in np.linspace(0, 360, num_spectra, endpoint=False)]

    for idx, obj in enumerate(speclist):
        color = colors[idx]  # Get the color for the current spectrum

        # Normalize the intensities using normconst
        C = obj.normalize()    
        intens_norm = obj.intens * obj.normconst

        # Temporary frequency array in case we have a scaling factor
        tmpfreq = np.copy(obj.freq) * obj.fscal

        # Check if obj.spec is not None and plot it as a continuous line
        if obj.spec is not None:
            x_values = np.linspace(obj.xmin, obj.xmax, len(obj.spec))
            y_values = obj.spec * obj.normconst
            ax.plot(x_values, y_values, color=color, linewidth=linewidth/2, label=f"Spectrum {idx + 1}")
            ax.fill_between(x_values, 0, y_values, color=color, alpha=0.25)

    # Set plot limits and labels (using the last object's attributes)
    ax.set_ylim(bottom=0)
    ax.set_xlim(speclist[-1].xmin, speclist[-1].xmax)
    ax.set_xlabel('Frequency (cm$^{-1}$)')
    ax.set_ylabel('Intensity' if speclist[-1].normconst == 1.0 else 'Relative Intensity')

    # Remove the title
    # ax.set_title('Vibrational Spectra')  # Commenting out to remove the title

    # Add the legend in the top-right corner
    ax.legend(loc='upper right')

    # Apply a more subtle grid style
    ax.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.7)

    plt.tight_layout()
    if save is not None:
        plt.savefig(save)
        print(f"Spectrum saved to {save}")
    plt.show()

