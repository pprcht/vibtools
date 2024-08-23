import numpy as np
from ase import Atoms
from ase.data import atomic_masses
from ase.units import Bohr,Hartree
from ._irtools import py_computespec_core, py_print_vib_spectrum_stdout

class IRtoolsCalculator:
    def __init__(self, atoms: Atoms, hessian: np.ndarray, dipole_gradient: np.ndarray, fscal: float = 1.0):
        """
        Initialize the IRtoolsCalculator with an ASE Atoms object, Hessian matrix, and dipole gradient matrix.

        Parameters:
        atoms (ASE Atoms object): The atomic structure.
        hessian (numpy.ndarray): The Hessian matrix (3*nat, 3*nat). Expected in Hartree/Bohr units
        dipole_gradient (numpy.ndarray): The dipole gradient matrix (3, 3*nat). expected in a.u. units
        fscal (float): The frequency scaling factor.
        """
        self.atoms = atoms
        self.hessian = hessian.astype(np.float64)
        self.dipole_gradient = dipole_gradient.astype(np.float64)
        self.hessian = np.ascontiguousarray(self.hessian)
        self.dipole_gradient = np.ascontiguousarray(self.dipole_gradient)
        self.fscal = fscal
        self.freq = None
        self.intens = None
        # Plotting parameters
        self.Cnorm = 1.0
        self.xmin = 100.0
        self.xmax = 4500.0
        self.dx = 1.0
        self.fwhm = 30.0
        self.spec = None

        # Initialize the amass array with atomic masses for elements 1-118
        self.amass = np.zeros(118, dtype=np.float64)
        for i in range(1, 119):
            self.amass[i-1] = atomic_masses[i]


    def compute(self):
        """
        Compute the vibrational spectrum using the Fortran routine.

        Returns:
        freq (numpy.ndarray): The computed frequencies.
        intens (numpy.ndarray): The computed intensities.
        """
        nat = len(self.atoms)
        at = self.atoms.get_atomic_numbers().astype(np.int32)
  
        # The Fortran/C++ code expects Bohr
        xyz = self.atoms.get_positions().astype(np.float64) / Bohr

        # Prepare output arrays
        self.freq = np.zeros(3 * nat, dtype=np.float64)
        self.intens = np.zeros(3 * nat, dtype=np.float64)

        # dipole gradient matrix needs transposing for the C++/Fortran passing
        dipole_gradient = np.ascontiguousarray(self.dipole_gradient.T)

        # Call the Fortran routine via the C++ wrapper
        py_computespec_core(nat, at, xyz, self.hessian, dipole_gradient, self.amass, self.fscal, self.freq, self.intens)

        return self.freq, self.intens


    def print(self):
        """
        Print the computed vibrational frequencies and intensities.
        """
        if self.freq is None or self.intens is None:
            self.compute()
        py_print_vib_spectrum_stdout( freq=self.freq, intens=self.intens)


    def plot(self, color='b-', linewidth=2, figsize=(9, 6), save=None):
        """
        Plot the computed vibrational spectrum as a stick spectrum.
        
        Parameters:
        - color: The color and line style of the sticks (default: 'b-').
        - linewidth: The width of the sticks (default: 2).
        - figsize: The size of the figure (default: (9, 6)).
        - save: A file name if given to which the spectrum is saved. (needs an extension)
        """
        if self.freq is None or self.intens is None:
            self.compute()
    
        import matplotlib.pyplot as plt
    
        # Normalize the intensities using Cnorm
        intens_norm = self.intens * self.Cnorm
    
        fig, ax = plt.subplots(figsize=figsize)
    
        # Create a stick plot using stem
        markerline, stemlines, baseline = ax.stem(self.freq, intens_norm, linefmt=color, basefmt=" ")
    
        # Adjust the markerline and stemlines
        plt.setp(markerline, 'marker', 'o')  # Remove markers at the top of the sticks
        plt.setp(stemlines, 'linewidth', linewidth)  # Set the linewidth of the sticks
    
        # Set plot limits and labels
        ax.set_ylim(bottom=0)  # Set the lower y-limit to exactly zero
        ax.set_xlim(self.xmin,self.xmax)
        ax.set_xlabel('Frequency (cm$^{-1}$)')
        ax.set_ylabel('Intensity')
        ax.set_title('Vibrational Spectrum')
    
        # Apply a more subtle grid style
        ax.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
    
        plt.tight_layout()
        if save is not None:
           plt.savefig(save)
           print(f"Spectrum saved to {save}")
        plt.show()
    

