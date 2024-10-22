import numpy as np
from typing import List, Optional
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from .calculator import vibtoolsCalculator

def process_file_multiplot(file_name, calculator_list, label_list = None):
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
    if label_list:
       label_list.append(vib_calc.filename)


def multiplot(speclist: List['vibtoolsCalculator'], linewidth: int = 2, 
              figsize: tuple = (9, 6), save: str = None,
              colormap: str = None, labels: Optional[List[Optional[str]]] = None,
              xlim: tuple = None, ylim: tuple = None) -> None:
    """
    Plot the computed vibrational spectrum as a stick spectrum for each object in speclist.

    Parameters:
    - speclist: A list of vibtoolsCalculator objects.
    - linewidth: The width of the sticks (default: 2).
    - figsize: The size of the figure (default: (9, 6)).
    - save: A file name to which the spectrum is saved (needs an extension).
    - colormap: choose a different coloring scheme
    - labels: An optional list of plot labels
    - xlim/ylim: Optional x/y axis limits
    """

    fig, ax = plt.subplots(figsize=figsize)

    # Get a color palette based on the number of spectra
    num_spectra = len(speclist)
    if colormap:
        # Use the provided colormap to generate colors
        cmap = plt.get_cmap(colormap)
        colors = [cmap(i / num_spectra) for i in range(num_spectra)]
    else:
        # Use the default matplotlib color cycle
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for idx, obj in enumerate(speclist):
        color = colors[idx]  # Get the color for the current spectrum

        # Normalize the intensities using normconst
        C = obj.normalize()    
        intens_norm = obj.intens * obj.normconst

        # Temporary frequency array in case we have a scaling factor
        tmpfreq = np.copy(obj.freq) * obj.fscal

        # Assign a label: Use the provided label or default to "Spectrum idx+1"
        if labels and len(labels) > idx and labels[idx] is not None:
            label = labels[idx]
        else:
            label = f"Spectrum {idx + 1}"

        # Check if obj.spec is not None and plot it as a continuous line
        if obj.spec is not None:
            x_values = np.linspace(obj.xmin, obj.xmax, len(obj.spec))
            y_values = obj.spec * obj.normconst
            ax.plot(x_values, y_values, color=color, linewidth=linewidth/2, label=label)
            ax.fill_between(x_values, 0, y_values, color=color, alpha=0.25)

    # Set plot limits and labels (using the last object's attributes)
    if ylim: 
       ax.set_ylim(ylim)
    ax.set_ylim(bottom=0)
    if xlim:
       ax.set_xlim(xlim) 
    else:
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

