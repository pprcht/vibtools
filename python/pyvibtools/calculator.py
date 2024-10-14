import numpy as np
from ase import Atoms
from ase.data import atomic_masses
from ase.units import Bohr,Hartree
from ._vibtools import py_computespec_core, py_print_vib_spectrum_stdout, py_lorentzian_broadening

from .readers import read_freqint, read_hessian, read_dipgrad, read_ASE, read_jdx
from .filetypes import check_ASE_readable
from .filetypes import FileFormatChecker,register_default_formats
from .utils import SIS,process_exp_spectrum

class vibtoolsCalculator:
    def __init__(self, atoms: Atoms=None, hessian: np.ndarray=None, 
                 dipole_gradient: np.ndarray=None, fscal: float = 1.0):
        """
        Initialize the vibtoolsCalculator with an ASE Atoms object, 
        Hessian matrix, and dipole gradient matrix.

        Parameters:
        atoms (ASE Atoms object): The atomic structure.
        hessian (numpy.ndarray): The Hessian matrix (3*nat, 3*nat). Expected in Hartree/Bohr units
        dipole_gradient (numpy.ndarray): The dipole gradient matrix (3, 3*nat). expected in a.u. units
        fscal (float): The frequency scaling factor.
        """
        self.atoms = atoms
        if hessian is not None:
           self.hessian = hessian.astype(np.float64)
           self.hessian = np.ascontiguousarray(self.hessian)
        else:
           self.hessian = None
        if dipole_gradient is not None:
           self.dipole_gradient = dipole_gradient.astype(np.float64)
           self.dipole_gradient = np.ascontiguousarray(self.dipole_gradient)
        else:
           self.dipole_gradient = None
        self.fscal = fscal
        self.freq = None
        self.intens = None
        self.filename = None
        self.expspec = False
        # Plotting parameters
        self.normconst = 1.0
        self.xmin = 100.0   # in cm⁻¹
        self.xmax = 5000.0  # in cm⁻¹
        self.dx = 1.0       # in cm⁻¹
        self.fwhm = 30.0    # for Lorentzian broadening
        self.IR_CUTOFF = 10.0  # in km/mol (just for determining IR activity)
        self.spec = None

        # Initialize the amass array with atomic masses for elements 1-118
        self.amass = np.zeros(118, dtype=np.float64)
        for i in range(1, 119):
            self.amass[i-1] = atomic_masses[i]


    def read(self, xyzfile=None, hessfile=None, dipfile=None, vibspecfile=None):
        """
        Read and overwrite data of a given vibtoolsCalculator
        """
        self.hessian = None
        self.dipole_gradient = None
        self.freq = None
        self.intens = None
        self.spec = None
        self.filename = None  
         
        if xyzfile is not None:
           if check_ASE_readable(xyzfile):
              self.atoms = read_ASE(xyzfile)

        if hessfile is not None:
           self.hessian = read_hessian(hessfile)
 
        if dipfile is not None:
           self.dipole_gradient = read_dipgrad(dipfile)

        if vibspecfile is not None:
          checker = FileFormatChecker()
          register_default_formats(checker)
          file_type = checker.check_format(vibspecfile)
      
          if file_type is None:
              print(f"File '{vibspecfile}' does not exist or has an unknown format.")
          elif file_type == "TM_vibspectrum":
              self.freq, self.intens = read_freqint(vibspecfile)  
              if self.freq is not None and self.intens is not None:
                 self.filename=vibspecfile 
          elif file_type == "JDX_experimental":
              print('JDX file')
              _, spectral_data = read_jdx(vibspecfile)
              self.freq, self.intens = zip(*spectral_data)
              if self.freq is not None and self.intens is not None:
                 self.filename=vibspecfile
              self.expspec = True 
              self.freq = np.array(self.freq)
              self.intens = np.array(self.intens)

    def compute(self):
        """
        Compute the vibrational spectrum using the Fortran routine.
        Does nothing for read-in experimental spectra

        Returns:
        freq (numpy.ndarray): The computed frequencies.
        intens (numpy.ndarray): The computed intensities.
        """
        if self.expspec:
          # Do nothing for experimental spectra
          return self.freq, self.intens
        if self.atoms is None:
          raise ValueError("The ASE atoms object is required for frequency calculation")
        if self.hessian is None:
          raise ValueError("Non-massweighted Hessian matrix is required for frequency calculation")
        if self.dipole_gradient is None:  
          raise ValueError("Dipole gradient is required for intensity calculation") 

        nat = len(self.atoms)
        at = self.atoms.get_atomic_numbers().astype(np.int32)
  
        # The Fortran/C++ code expects Bohr
        xyz = self.atoms.get_positions().astype(np.float64) / Bohr

        # Prepare output arrays
        self.freq = np.zeros(3 * nat, dtype=np.float64)
        self.intens = np.zeros(3 * nat, dtype=np.float64)

        # dipole gradient matrix needs transposing for the C++/Fortran passing
        dipole_gradient = np.ascontiguousarray(self.dipole_gradient.T)

        # Copy the Hessian as to NOT OVERWRITE IT
        hessian = np.copy(self.hessian)
       
        # Call the Fortran routine via the C++ wrapper
        py_computespec_core(nat, at, xyz, hessian, dipole_gradient, 
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
        # Temporary frequency array in case we have a scaling factor
        tmpfreq = np.copy(self.freq) * self.fscal
        npoints = int(np.round(np.abs(self.xmin - self.xmax) / self.dx)) + 1
        self.spec = np.ascontiguousarray(np.zeros(npoints, dtype=np.float64))
        if self.expspec:
           spectral_data = list(zip(self.freq,self.intens))
           spectral_data = process_exp_spectrum(spectral_data, 
                                                self.dx, xmin=self.xmin,
                                                xmax=self.xmax)
           _,self.spec = zip(*spectral_data)
           self.spec = np.array(self.spec)
        else:
           py_lorentzian_broadening(nmodes=nmodes, freq=tmpfreq, intens=self.intens,
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
    
        # Normalize the intensities using normconst
        intens_norm = self.intens * self.normconst

        # Temporary frequency array in case we have a scaling factor    
        tmpfreq = np.copy(self.freq) * self.fscal

        fig, ax = plt.subplots(figsize=figsize)
    
        # Check if self.spec is not None and plot it as a continuous line
        if self.spec is not None:
            x_values = np.linspace(self.xmin, self.xmax, len(self.spec))
            y_values = self.spec*self.normconst 
            ax.plot(x_values, y_values, color, linewidth=linewidth/2) 
            ax.fill_between(x_values, 0, y_values, color=color, alpha=0.25)

        if sticks:
          # Create a stick plot using stem
          markerline, stemlines, baseline = ax.stem(tmpfreq, intens_norm, 
          linefmt=color, basefmt=" ")
          # Adjust the markerline and stemlines
          plt.setp(markerline, 'marker', stickmarker)  # Remove markers at the top of the sticks
          plt.setp(stemlines, 'linewidth', linewidth)  # Set the linewidth of the sticks

        # Set plot limits and labels
        ax.set_ylim(bottom=0)  # Set the lower y-limit to exactly zero
        ax.set_xlim(self.xmin,self.xmax)
        ax.set_xlabel('Frequency (cm$^{-1}$)')
        if self.normconst == 1.0:
          ax.set_ylabel('Intensity')
        else:
          ax.set_ylabel('Relative Intensity')

        if self.filename is not None:
          ax.set_title(f'Vibrational Spectrum: {self.filename}')
        else:
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

    def normalize(self, scheme='integral_root'):
        """
        Normalize the spectrum according to the selected scheme.
    
        Parameters:
        - scheme: A string indicating the normalization scheme to use. 
                  The default is 'integral_root', which normalizes the spectrum 
                  to the root of its numerical integral. Other options can be 
                  added as needed.
        """
        if self.freq is None or self.intens is None:
            self.compute()
        if self.spec is None:
            self.broaden()
        
        # Initialize nfac and tmpa
        self.normconst = 1.0
    
        if scheme == 'integral_root':
            # Sum the elements in self.spec and compute the normalization factor
            tmpa = 0.0
            for value in self.spec:
                tmpa += (value**2) * self.dx
            tmpa = np.sqrt(tmpa)
            self.normconst = 1.0 / tmpa
        
        elif scheme == 'max_value':
            # Normalize to the maximum value in the spectrum
            max_val = max(self.spec)
            self.normconst = 1.0 / max_val
    
        elif scheme == 'sum':
            # Normalize by the sum of the intensities
            sum_val = sum(self.spec) * self.dx
            self.normconst = 1.0 / sum_val
      
        elif scheme == 'msc':
            # A normalization scheme that was in the orignal specmatch code     
            # I think it should be equivalent to the sum scheme if dx=1, no idae why it exists
            summe = 0.0
            sqrt_spec = np.sqrt(self.spec)
            for value in sqrt_spec:
                summe += (value**2)
            self.normconst = 1.0 / np.sqrt(summe)

        else:
            raise ValueError(f"Unknown normalization scheme: {scheme}")
    
        return self.normconst


#########################################################################################


def matchscore(calc1, calc2):
    """
    Calculate the match score (r_msc), Euclidean norm (r_euc), 
    and Pearson correlation coefficient (r_pcc) between 
    two vibtoolsCalculator objects using their self.spec arrays.
    
    Parameters:
    - calc1: The first vibtoolsCalculator object.
    - calc2: The second vibtoolsCalculator object.
    
    Returns:
    - A dictionary containing r_msc, r_euc, r_pcc
    """
    # Ensure the first calculator's spectrum is computed and broadened
    if calc1.freq is None or calc1.intens is None:
        calc1.compute()
    if calc1.spec is None:
        calc1.broaden()
    # The two normalizations should be consistent with the newspecmatch code
    C1 = calc1.normalize()
    u = (np.array(np.copy(calc1.spec), dtype=np.float64) * C1)
    y_pred = u/np.sum(u) # different normalization for SIS

    # Ensure the second calculator's spectrum is computed and broadened
    if calc2.freq is None or calc2.intens is None:
        calc2.compute()
    if calc2.spec is None:
        calc2.broaden()
    # Ensure the second calculator's spectrum is computed and broadened
    C2 = calc2.normalize()
    v = (np.array(np.copy(calc2.spec), dtype=np.float64) * C2)
    y_target = v/np.sum(v) # different normalization for SIS

    # Get the sqrt of the spectra arrays (consistency with the newspecmatch code)
    u = (np.sqrt(u))
    v = (np.sqrt(v))

    # Initialize sums
    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    sum4 = 0.0
    
    # Iterate through u and v simultaneously
    for u_i, v_i in zip(u, v):
        sum1 += u_i * v_i
        sum2 += u_i**2
        sum3 += v_i**2
        sum4 += (v_i - u_i)**2 

    # Calculate the Matchscore (msc) using the Cauchy-Schwarz inequality
    r_msc = (sum1**2) / (sum2 * sum3) if sum2 * sum3 != 0 else 0.0
    
    # Calculate the Euclidean norm (euc)
    #r_euc = np.sqrt(sum4 / sum2) if sum2 != 0 else 0.0
    r_euc = 1.0/(1.0 + (sum4/sum2)) if sum2 != 0 else 1.0 #as in original newspecmatch code

    # Pearson Correlation Coefficient (r_pcc)
    u_mean = np.mean(u)
    v_mean = np.mean(v)
    numerator_pcc = np.sum((u - u_mean) * (v - v_mean))
    denominator_pcc = np.sqrt(np.sum((u - u_mean) ** 2) * np.sum((v - v_mean) ** 2))
    r_pcc = numerator_pcc / denominator_pcc

    # Spectral Information Similarity Metric from https://doi.org/10.1021/acs.jcim.1c00055
    # Note, the spectra broadening is different compared to this publication!
    r_sis = SIS(y_pred, y_target)


    return r_msc, r_euc, r_pcc, r_sis


