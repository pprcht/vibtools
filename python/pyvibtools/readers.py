# pyvibtools/readers.py
#
# Attempts to implements a variety of reader functions for common file types
# from atomic simulation software. Will be expanded over time.
#
import numpy as np
from ase.io import read
from .filetypes import FileFormatChecker,register_default_formats

##########################################################################################
####################### auto format reader interfaces ####################################
##########################################################################################
def read_freqint(filename):
    """
    Read a file and return frequency and intensity pairs
    """
    checker = FileFormatChecker()
    register_default_formats(checker)
    file_type = checker.check_format(filename)

    if file_type is None:
        print(f"File '{filename}' does not exist or has an unknown format.")
    elif file_type == "TM_vibspectrum":
        wave_numbers, intensities = read_vibspectrum(filename)

    return np.array(wave_numbers), np.array(intensities)



def read_hessian(filename):
    """
    Read a file and return a Hessian matrix
    """
    checker = FileFormatChecker()
    register_default_formats(checker) 
    file_type = checker.check_format(filename)

    if file_type is None:
        print(f"File '{filename}' does not exist or has an unknown format.")
    elif file_type == "TM_hessian":
        hessian = read_tm_hessian(filename)
    elif file_type == 'plain':
        hessian = read_plain_hessian(filename)

    return hessian



def read_dipgrad(filename):
    """
    Read a file and return a dipole derivative matrix (3 x 3*nat elements)
    """
    checker = FileFormatChecker()
    register_default_formats(checker)
    file_type = checker.check_format(filename)

    if file_type is None:
        print(f"File '{filename}' does not exist or has an unknown format.")
    elif file_type == "plain":
        dipgrad = read_plain_dipgrad(filename)

    return dipgrad


##########################################################################################
##########################################################################################
##########################################################################################

def read_vibspectrum(filename):
    """
    Reads a Turbomole vibspectrum file and extracts the wave number and IR intensity columns.

    Parameters:
    - filename: Path to the vibspectrum file.

    Returns:
    - A tuple of two lists: (wave_numbers, intensities)
    """
    wave_numbers = []
    intensities = []
    in_block = False  # Flag to indicate whether we are inside the $vibrational spectrum block
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("$vibrational"):
                in_block = True
                continue
            if line.startswith("$end"):
                in_block = False
                continue
            if in_block:
                if line.startswith("#") or not line:
                    # Skip comment lines and empty lines
                    continue
                # Split the line into columns
                columns = line.split()
                # Determine the format based on the number of columns
                if len(columns) >= 4: 
                    try:
                        float(columns[1])
                        _, wave_number, intensity, *_ = columns
                    except ValueError:  # mode, symmetry, wave number, IR intensity
                        _, _, wave_number, intensity, *_ = columns
                elif len(columns) == 3:  # mode, wave number, IR intensity
                    _, wave_number, intensity = columns
                # Convert to float and append to the lists
                wave_numbers.append(float(wave_number))
                intensities.append(float(intensity))
    return wave_numbers, intensities

#########################################################################################
#########################################################################################

def read_plain_hessian(filename):
    """
    Reads a plain text Hessian matrix file containing numbers and converts it into a 2D matrix 
    with dimensions 3*nat x 3*nat, where nat is derived from the number of elements.

    Parameters:
    - filename: Path to the Hessian file.

    Returns:
    - A 2D NumPy array with dimensions 3*nat x 3*nat.
    """
    # Read the file and extract all numbers into a list
    hessian_data = []
    
    with open(filename, 'r') as file:
        for line in file:
            # Split the line into components and convert them to floats
            hessian_data.extend(map(float, line.split()))
    # Calculate the total number of elements
    total_elements = len(hessian_data)
    # Determine the number of atoms (nat)
    if total_elements % 9 != 0:
        raise ValueError("Total number of elements is not divisible by 9. File may be malformed.")
    nat = int(np.sqrt(total_elements // 9))  # Since we need (3*nat) x (3*nat) elements
    if total_elements != 9 * nat * nat:
        raise ValueError(f"Total number of elements ({total_elements}) does not match the expected (3*nat) x (3*nat) structure.")
    # Reshape the hessian_data list into a 2D NumPy array with dimensions 3*nat x 3*nat
    hessian_matrix = np.array(hessian_data).reshape(3 * nat, 3 * nat)
    
    return hessian_matrix



def read_tm_hessian(filename):
    """
    Reads a Turbomole Hessian file and extracts the Hessian matrix.

    Parameters:
    - filename: Path to the Hessian file.

    Returns:
    - A 2D NumPy array representing the Hessian matrix.
    """
    hessian_data = []
    in_block = False  # Flag to indicate whether we are inside the $hessian block
    
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("$hessian"):
                in_block = True
                continue
            if line.startswith("$end"):
                in_block = False
                continue
            if in_block:
                if not line:
                    # Skip empty lines
                    continue
                # Split the line into floating-point numbers and extend the hessian_data list
                hessian_data.extend(map(float, line.split()))

    # Convert the list of hessian_data to a 2D NumPy array
    # Since the matrix is 3*nat x 3*nat, the number of elements must be a perfect square
    size = int(np.sqrt(len(hessian_data)))
    if size * size != len(hessian_data):
        raise ValueError("The number of elements in the Hessian data does not form a square matrix.")
    
    hessian_matrix = np.array(hessian_data).reshape(size, size)
    
    return hessian_matrix

##########################################################################################
##########################################################################################
##########################################################################################

def read_plain_dipgrad(filename):
    """
    Reads a plain text dipgrad file containing numbers and converts it into a 2D matrix 
    with dimensions 3 x 3*nat, where nat is derived from the number of elements.

    Parameters:
    - filename: Path to the dipgrad file.

    Returns:
    - A 2D NumPy array with dimensions 3 x 3*nat.
    """
    # Read the file and extract all numbers into a list
    dipgrad_data = []
    
    with open(filename, 'r') as file:
        for line in file:
            # Split the line into components and convert them to floats
            dipgrad_data.extend(map(float, line.split()))
    
    # Calculate the total number of elements
    total_elements = len(dipgrad_data)
    
    # Determine the number of atoms (nat)
    if total_elements % 3 != 0:
        raise ValueError("Total number of elements is not divisible by 3. File may be malformed.")
    
    nat = total_elements // 9  # Since we need 3 x 3*nat elements
    
    if total_elements != 3 * 3 * nat:
        raise ValueError(f"Total number of elements ({total_elements}) does not match the expected 3 x 3*nat structure.")

    # Reshape the dipgrad_data list into a 2D NumPy array with dimensions 3 x 3*nat
    dipgrad_matrix = np.array(dipgrad_data).reshape(3 * nat , 3)
    dipgrad_matrix = np.ascontiguousarray(dipgrad_matrix.T)
    
    return dipgrad_matrix

##########################################################################################
##########################################################################################
##########################################################################################

def read_ASE(filename):
    structure = read(filename)
    return structure   

##########################################################################################
##########################################################################################
##########################################################################################

def read_jdx(filename):
    with open(filename, 'r') as f:
        metadata = {}
        spectral_data = []
        deltax = None
        in_spectral_data = False

        for line in f:
            line = line.strip()

            # Handle metadata (key-value pairs)
            if '=' in line and not in_spectral_data:
                key, value = line.split('=', 1)
                metadata[key.strip('# ').upper()] = value.strip()

                # Capture DELTAX value for interpolation
                if 'DELTAX' in metadata:
                    deltax = float(metadata['DELTAX'])

                # Start collecting spectral data
                if '##XYDATA' in line or '##PEAKTABLE' in line:
                    in_spectral_data = True
                    continue

            # Handle spectral data (X followed by multiple Y values)
            elif in_spectral_data:
                if line.startswith('##END'):
                    in_spectral_data = False
                    continue

                # Parse the X value and subsequent Y values
                values = line.split()
                x_value = float(values[0])  # First value is the X value
                y_values = [float(y) for y in values[1:]]  # Remaining values are Y values

                # Calculate corresponding X values based on DELTAX
                for i, y in enumerate(y_values):
                    x = x_value + i * deltax
                    spectral_data.append((x, y))

        return metadata, spectral_data

