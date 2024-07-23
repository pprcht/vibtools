#!/usr/bin/env python
import argparse
from mace.calculators import MACECalculator
from ase import build
from ase.io import read, write
from ase.units import Hartree, Bohr, Ang
from ase.optimize import LBFGS
import numpy as np
import time

au_to_debye = 2.54174692380177980056


##################################################################################################

def get_dipole(atoms, calc):
    # implelemtation goes here
    calc.get_property('dipole', atoms=atoms)
    dipole = calc.results['dipole']
    return dipole

##################################################################################################

def setup_singlepoint(input_file, model_type='DipoleMACE', model_path=None):
    atoms = read(input_file)
    calc = MACECalculator(model_path, device='cpu', default_dtype='float64', model_type=model_type)
    start=time.time()
    dipole = get_dipole(atoms, calc)
    finish=time.time()
    totdipole = np.linalg.norm(dipole)
    dipole = dipole 
    print(f"Molecular dipole moment (x,y,z): {dipole[0]} {dipole[1]} {dipole[2]} \n")
    print(f"Total dipole moment: {totdipole} \n")
    print(f"time {finish-start:.3f} s")
    return atoms,calc

###################################################################################################

def calculate_numerical_hessian(atoms, calc, displace=1e-3):
    """
    Calculate Cartesian dipole derivative matrix for a given set of atoms using central differences.
    
    Parameters:
        atoms (ase.Atoms): The Atoms object with an attached calculator.
        calc (MACECalculator): The MACE calculator
        dx (float): The displacement in Angstroms for numerical derivatives.
    
    Returns:
        np.ndarray: The numerical dipole derivatives
    """
    n_atoms = len(atoms)
    dipgrad = np.zeros((3 , 3*n_atoms))
    dx = displace * Bohr

    # Save the original positions
    original_positions = atoms.get_positions().copy()

    print("\nNumerical Dipole Derivatives")
    print("-----------------")
    print(f"Displacement : {displace:.4e} a.u. ({dx:.4e} AA)")
    start=time.time()
    for i in range(n_atoms):
        print(f"Displacing atom {i+1} / {n_atoms}")
        for j in range(3):
            # Positive displacement
            atoms.positions[i][j] += dx
            #forces_plus = -atoms.get_forces()
            dipole_plus = get_dipole(atoms, calc) / au_to_debye            

            # Negative displacement
            atoms.positions[i][j] -= 2 * dx
            #forces_minus = -atoms.get_forces()
            dipole_minus = get_dipole(atoms, calc) / au_to_debye
            
            # Reset the position
            atoms.positions[i][j] = original_positions[i][j]
            
            # Fill the numerical dipole derivatives
            ii = i*3+j
            for l in range(3): 
               dipgrad[l,ii] = (dipole_plus[l] - dipole_minus[l]) / (2*dx)

    finish=time.time()    
    print(f"Done. {finish-start:.3f} s")

    return dipgrad

##################################################################################################

def write_dipgrad(dipgrad, filename):
    """
    Write plaintext file with dipole gradients
    
    Args:
    dipgrad (np.ndarray): A 3 by  3*n_atoms matrix.
    filename (str): The name of the file to write to.
    """
    n_atoms3 = dipgrad.shape[1] 
    with open(filename, 'w') as f:
        for i in range(n_atoms3):
            for j in range(3):
                f.write(f"{dipgrad[j, i]:25.15f}")
            f.write("\n")

##################################################################################################

def main():
    parser = argparse.ArgumentParser(description='Run different MACE tasks via ASE.')
    parser.add_argument('input_file', help='The path to the input file.')
    parser.add_argument('--hess','--dmu', action='store_true', help='Calculate the numerical dipole derivatives')
    parser.add_argument('-dx', default=1.0e-3 , type=float, help='Set numerical Hessian displacement (in Bohr)')
    args = parser.parse_args()

    print("\n MACE script")
    print(" ===========\n")

    #model_path = '/home/philipp/code/MACE-mu/models/SPICE_medium_dipole.model'
    model_path = '/home/philipp/code/MACE-mu/models/SPICE_small_dipole.model'
    model_type = 'DipoleMACE'

    # Load atoms
    atoms,calc = setup_singlepoint(args.input_file, model_type, model_path)

    wrgeo=False
    wrhess=False
    # Run calculations based on the specified run type
    if args.hess:
       dipole_derivatives = calculate_numerical_hessian(atoms, calc, displace=args.dx)
       print("\nDipole derivatives")
       print(dipole_derivatives)
       write_dipgrad(dipole_derivatives, f"dipgrad.{model_type}")
       wrhess=True

    
    if wrgeo or wrhess:
        print("\nFiles written:")
        if wrhess:
           print(f" - dipgrad.{model_type}")

#################################################################################################
if __name__ == '__main__':
    main()

