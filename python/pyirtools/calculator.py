import numpy as np
from ase import Atoms
from ase.data import atomic_masses
from ase.units import Bohr,Hartree
from ._pyirtools import computespec_core

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
        self.fscal = fscal
        self.freq = None
        self.intens = None

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

        # Call the Fortran routine via the C++ wrapper
        computespec_core(nat, at, xyz, self.hessian, self.dipole_gradient, self.amass, self.fscal, self.freq, self.intens)

        return self.freq, self.intens

    def plot(self):
        """
        Plot the computed vibrational spectrum.
        """
        if self.freq is None or self.intens is None:
            self.compute()

        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 6))
        plt.plot(self.freq, self.intens, 'b-', lw=2)
        plt.xlabel('Frequency (cm^-1)')
        plt.ylabel('Intensity')
        plt.title('Vibrational Spectrum')
        plt.grid(True)
        plt.show()

    def print(self):
        """
        Print the computed vibrational frequencies and intensities.
        """
        if self.freq is None or self.intens is None:
            self.compute()

        print("Frequencies (cm^-1) and Intensities:")
        for f, i in zip(self.freq, self.intens):
            print(f"Frequency: {f:.2f} cm^-1, Intensity: {i:.2f}")

