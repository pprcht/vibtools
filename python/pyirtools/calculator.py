import numpy as np
from ase import Atoms
from ase.data import atomic_masses
from ase.units import Bohr,Hartree
from ._irtools import py_computespec_core, py_print_vib_spectrum_stdout, py_lorentzian_broadening

class IRtoolsCalculator:
    def __init__(self, atoms: Atoms=None, hessian: np.ndarray=None, dipole_gradient: np.ndarray=None, fscal: float = 1.0):
        """
        Initialize the IRtoolsCalculator with an ASE Atoms object, 
        Hessian matrix, and dipole gradient matrix.

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
        py_computespec_core(nat, at, xyz, self.hessian, dipole_gradient, 
                            self.amass, self.fscal, self.freq, self.intens)
        return self.freq, self.intens


    def print(self):
        """
        Print the computed vibrational frequencies and intensities.
        """
        if self.freq is None or self.intens is None:
            self.compute()
        py_print_vib_spectrum_stdout( freq=self.freq, intens=self.intens)



    def broaden(self, xmin=None, xmax=None, dx=None, fwhm=None):
        """
        Apply Lorentzian functions to make a stick spectrum into a continuous one.
        
        Parameters:
        - xmin: Minimum x value for broadening. Overwrites self.xmin if provided.
        - xmax: Maximum x value for broadening. Overwrites self.xmax if provided.
        - dx: Step size for x values. Overwrites self.dx if provided.
        - fwhm: Full width at half maximum for the Lorentzian broadening. 
                Overwrites self.fwhm if provided.
        """
        if self.freq is None or self.intens is None:
            self.compute()
    
        # Overwrite self attributes if arguments are provided
        if xmin is not None:
            self.xmin = xmin
        if xmax is not None:
            self.xmax = xmax
        if dx is not None:
            self.dx = dx
        if fwhm is not None:
            self.fwhm = fwhm
    
        nmodes = len(self.freq)
        npoints = int(np.round(np.abs(self.xmin - self.xmax) / self.dx))
        self.spec = np.ascontiguousarray(np.zeros(npoints, dtype=np.float64))
        py_lorentzian_broadening(nmodes=nmodes, freq=self.freq, intens=self.intens,
                                 xmin=self.xmin, xmax=self.xmax, dx=self.dx,
                                 fwhm=self.fwhm, npoints=npoints, plt=self.spec)
    

    def plot(self, color='b', linewidth=2, figsize=(9, 6), save=None, 
             sticks=False, stickmarker='None'):
        """
        Plot the computed vibrational spectrum as a stick spectrum.
        
        Parameters:
        - color: The color and line style of the sticks (default: 'b').
        - linewidth: The width of the sticks (default: 2).
        - figsize: The size of the figure (default: (9, 6)).
        - save: A file name if given to which the spectrum is saved. (needs an extension)
        """
        if self.freq is None or self.intens is None:
            self.compute()
        if self.spec is None:
            self.broaden() 
    
        import matplotlib.pyplot as plt
    
        # Normalize the intensities using Cnorm
        intens_norm = self.intens * self.Cnorm
    
        fig, ax = plt.subplots(figsize=figsize)
    
        # Check if self.spec is not None and plot it as a continuous line
        if self.spec is not None:
            x_values = np.linspace(self.xmin, self.xmax, len(self.spec))
            y_values = self.spec*self.Cnorm 
            ax.plot(x_values, y_values, color, linewidth=linewidth/2) 
            ax.fill_between(x_values, 0, y_values, color=color, alpha=0.25)

        if sticks:
          # Create a stick plot using stem
          markerline, stemlines, baseline = ax.stem(self.freq, intens_norm, 
          linefmt=color, basefmt=" ")
          # Adjust the markerline and stemlines
          plt.setp(markerline, 'marker', stickmarker)  # Remove markers at the top of the sticks
          plt.setp(stemlines, 'linewidth', linewidth)  # Set the linewidth of the sticks

        # Set plot limits and labels
        ax.set_ylim(bottom=0)  # Set the lower y-limit to exactly zero
        ax.set_xlim(self.xmin,self.xmax)
        ax.set_xlabel('Frequency (cm$^{-1}$)')
        if self.Cnorm == 1.0:
          ax.set_ylabel('Intensity')
        else:
          ax.set_ylabel('Relative Intensity')
        ax.set_title('Vibrational Spectrum')
    
        # Apply a more subtle grid style
        ax.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
    
        plt.tight_layout()
        if save is not None:
           plt.savefig(save)
           print(f"Spectrum saved to {save}")
        plt.show()
 
        # Return pyplot objects in case we want to add more later
        return fig,ax     

    def normalize(self):
        """
        Normalize the spectrum to the root of its numerical integral 
        """
        if self.freq is None or self.intens is None:
            self.compute()
        if self.spec is None:
            self.broaden()
        # Initialize nfac and tmpa
        self.Cnorm = 1.0
        tmpa = 0.0
        # Sum the elements in self.spec
        for value in self.spec:
            tmpa += value*self.dx
        # Compute the normalization factor
        tmpa = np.sqrt(tmpa)
        self.Cnorm = 1.0 / tmpa
        return self.Cnorm


#########################################################################################


def matchscore(calc1, calc2):
    """
    Calculate the match score (r_msc), Euclidean norm (r_euc), 
    and Pearson correlation coefficient (r_pcc) between 
    two IRtoolsCalculator objects using their self.spec arrays.
    
    Parameters:
    - calc1: The first IRtoolsCalculator object.
    - calc2: The second IRtoolsCalculator object.
    
    Returns:
    - A dictionary containing r_msc, r_euc, r_pcc
    """
    # Ensure the first calculator's spectrum is computed and broadened
    if calc1.freq is None or calc1.intens is None:
        calc1.compute()
    if calc1.spec is None:
        calc1.broaden()
    C1 = calc1.normalize()

    # Ensure the second calculator's spectrum is computed and broadened
    if calc2.freq is None or calc2.intens is None:
        calc2.compute()
    if calc2.spec is None:
        calc2.broaden()
    C2 = calc2.normalize()

    # Get the spectra arrays (both normalized)
    u = calc1.spec * C1
    v = calc2.spec * C2

    # Match Score (r_msc)
    numerator_msc = np.sum(u * v) ** 2
    denominator_msc = np.sum(u ** 2) * np.sum(v ** 2)
    r_msc = numerator_msc / denominator_msc

    # Euclidean Norm (r_euc)
    numerator_euc = np.sum((u - v) ** 2)
    denominator_euc = np.sum(v ** 2)
    r_euc = (1.0 + numerator_euc / denominator_euc) ** -1

    # Pearson Correlation Coefficient (r_pcc)
    u_mean = np.mean(u)
    v_mean = np.mean(v)
    numerator_pcc = np.sum((u - u_mean) * (v - v_mean))
    denominator_pcc = np.sqrt(np.sum((u - u_mean) ** 2) * np.sum((v - v_mean) ** 2))
    r_pcc = numerator_pcc / denominator_pcc

    return r_msc, r_euc, r_pcc


