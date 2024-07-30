#!/usr/bin/env python
import argparse
from mace.calculators import mace_off
from ase import build
from ase.io import read, write
from ase.units import Hartree, Bohr, Ang
from ase.optimize import LBFGS
import numpy as np
import time

au_to_debye = 2.54174692380177980056

##################################################################################################

def setup_singlepoint(input_file, model_type='medium'):
    atoms = read(input_file)
    calc = mace_off(model=model_type, default_type="float64", device='cpu')
    atoms.calc = calc
    print(f"\nInitial single point energy (Hartree): {atoms.get_potential_energy() / Hartree}") 
    dipole = get_dipole()
    totdipole = np.linalg.norm(dipole, ord=2) * au_to_debye
    #print(f"Molecular dipole moment (x,y,z)/au, tot/Debye : {dipole}, {totdipole} \n")
    return atoms,calc

##################################################################################################

def run_optimization(atoms, model_type, rmsforce_au=1e-3):
    # Use rmsforce_au as an input argument for convergence condition
    rmsforce_eVAA = rmsforce_au * Hartree / Bohr
    
    print("\nGeometry optimization")
    print("---------------------")
    print(f"RMS force (Ha/Bohr, eV/AA): {rmsforce_au}, {rmsforce_eVAA}\n")
    optimizer = LBFGS(atoms)
    optimizer.run(fmax=rmsforce_eVAA)
    atoms.write(f'{model_type}.xyz')
    energy = atoms.get_potential_energy() / Hartree
    print(f"Final single point energy (Hartree): {energy}\n")
    return energy

##################################################################################################

def get_dipole( ):
    # implelemtation goes here
    dipole = (0.0, 0.0, 0.0)
    return dipole

###################################################################################################

def calculate_numerical_hessian(atoms, displace=1e-3):
    """
    Calculate and symmetrize the numerical Hessian matrix for a given set of atoms using central differences.
    
    Parameters:
        atoms (ase.Atoms): The Atoms object with an attached calculator.
        dx (float): The displacement in Angstroms for numerical derivatives.
    
    Returns:
        np.ndarray: The symmetrized Hessian matrix.
    """
    n_atoms = len(atoms)
    hessian = np.zeros((3 * n_atoms, 3 * n_atoms))
    dipgrad = np.zeros((3 , 3*n_atoms))
    dx = displace * Bohr


    # Save the original positions
    original_positions = atoms.get_positions().copy()

    print("\nNumerical Hessian")
    print("-----------------")
    print(f"Displacement : {displace:.4e} a.u. ({dx:.4e} AA)")
    start=time.time()
    for i in range(n_atoms):
        print(f"Displacing atom {i+1} / {n_atoms}")
        for j in range(3):
            # Positive displacement
            atoms.positions[i][j] += dx
            forces_plus = -atoms.get_forces()
            #dipole_plus = get_dipole()            

            # Negative displacement
            atoms.positions[i][j] -= 2 * dx
            forces_minus = -atoms.get_forces()
            #dipole_minus = get_dipole()
            
            # Reset the position
            atoms.positions[i][j] = original_positions[i][j]
            
            # Compute second derivatives
            derivative = (forces_plus - forces_minus) / (2 * dx)
            
            # convert eV/AA to au
            derivative = derivative / Hartree * Bohr *Bohr           

            # Fill the Hessian matrix: each element is a second derivative
            for k in range(n_atoms):
                for l in range(3):
                    hessian[3 * i + j, 3 * k + l] = derivative[k][l]

            # Fill the numerical dipole derivatives
            #ii = i*3+j
            #for l in range(3): 
            #   dipgrad[l,ii] = (dipole_plus[l] - dipole_minus[l]) / (2*dx)

    finish=time.time()    
    print(f"Done. {finish-start:.3f} s")
    # Symmetrize the Hessian matrix
    hessian = 0.5 * (hessian + hessian.T)

    return hessian#,dipgrad

###################################################################################################

def calculate_autograd_hessian(atoms,calc):
    """
    Calculate and symmetrize the numerical Hessian matrix for a given set of atoms using central differenc
    
    Parameters:
        atoms (ase.Atoms): The Atoms object with an attached calculator.
        calc  (MACECalculator): MACE calculator    

    Returns:
        np.ndarray: The Hessian matrix.
    """
    n_atoms = len(atoms)
    hessian = np.zeros((3 * n_atoms, 3 * n_atoms))
    #dipgrad = np.zeros((3 , 3*n_atoms))

    print("\nAutograd Hessian")
    print("-----------------")
    print("Calculating Hessian via autograd...")
    start=time.time()
    hessian_autograd = calc.get_hessian(atoms=atoms)
    finish=time.time()
    # reshape
    hessian = hessian_autograd.reshape((3 * n_atoms, 3 * n_atoms))

    # transform from eV/AA^2 to a.u.
    hessian = hessian / Hartree * Bohr * Bohr
 
    print(f"Done. {finish-start:.3f} s")
    # Symmetrize the Hessian matrix
    hessian = 0.5 * (hessian + hessian.T)

    return hessian

##################################################################################################

def write_tm_hessian(hessian, filename):
    """
    Write the hessian matrix to a file in the Turbomole format.
    
    Args:
    hessian (np.ndarray): A 3*n_atoms by 3*n_atoms matrix.
    filename (str): The name of the file to write to.
    """
    n_atoms3 = hessian.shape[0]  # Assuming hessian is a square matrix
    with open(filename, 'w') as f:
        f.write(' $hessian\n')
        for i in range(n_atoms3):
            k = 0
            for j in range(n_atoms3):
                k += 1
                # Print up to four elements on the same line without a newline character
                if k <= 4:
                    f.write(f"{hessian[i, j]:20.10f}")
                else:
                    # After four elements, print with a newline and reset counter
                    f.write(f"{hessian[i, j]:20.10f}\n")
                    k = 0
            # If the last line of the loop didn't end with a newline (not multiple of 4)
            if k != 0:
                f.write('\n')
        f.write(' $end\n')

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
    parser.add_argument('-s', '--small', action='store_true', help='Use the small model')
    parser.add_argument('-m', '--medium', action='store_true', help='Use the medium model')
    parser.add_argument('-l', '--large', action='store_true', help='Use the large model')
    parser.add_argument('--opt', action='store_true', help='Run geometry optimization')
    parser.add_argument('--hess', action='store_true', help='Calculate the numerical Hessian')
    parser.add_argument('--ohess', action='store_true', help='Optimize and then calculate the Hessian')
    parser.add_argument('-dx', default=1.0e-3 , type=float, help='Set numerical Hessian displacement (in Bohr)')
    parser.add_argument('--num', action='store_true', help='Select NUMERICAL Hessian calculation')
    args = parser.parse_args()

    print("\n MACE script")
    print(" ===========\n")

    # Determine model size
    if args.small:
        model_type = 'small'
    elif args.medium:
        model_type = 'medium'
    elif args.large:
        model_type = 'large'
    else:
        # medium model as default
        model_type = 'medium' 

    # Load atoms
    atoms,calc = setup_singlepoint(args.input_file, model_type)

    wrgeo=False
    wrhess=False
    # Run calculations based on the specified run type
    if args.opt or args.ohess:
        run_optimization(atoms, model_type)
        atoms.write('test.json')
        atoms.write('test.json')
        wrgeo=True 
    if args.hess or args.ohess:
        if args.num:
          #hessian_matrix,dipole_derivatives = calculate_numerical_hessian(atoms, displace=args.dx)
          hessian_matrix = calculate_numerical_hessian(atoms, displace=args.dx)
          print("\nHessian Matrix (numerical):")
          print(hessian_matrix)
          write_tm_hessian(hessian_matrix, f"numhess.{model_type}")
          #print("\nDipole derivatives")
          #print(dipole_derivatives)
          #write_dipgrad(dipole_derivatives, f"dipgrad.{model_type}")
          wrhess=True
        else:
          hessian_matrix = calculate_autograd_hessian(atoms,calc)
          print("\nHessian matrix (autograd):")
          print(hessian_matrix)
          write_tm_hessian(hessian_matrix, f"autogradhess.{model_type}")
          wrhess=True

    
    if wrgeo or wrhess:
        print("\nFiles written:")
        if wrgeo:
           print(f" - {model_type}.xyz")
        if wrhess:
           if args.num:
             print(f" - numhess.{model_type}")
             #print(f" - dipgrad.{model_type}")
           else: 
             print(f" - autogradhess.{model_type}")

#################################################################################################
if __name__ == '__main__':
    main()

