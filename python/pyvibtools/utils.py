import numpy as np


# Implement metrics SID and SIS from https://doi.org/10.1021/acs.jcim.1c00055
def SID(y_pred, y_target):
    # Ensure inputs are numpy arrays
    y_pred = np.asarray(y_pred)
    y_target = np.asarray(y_target)

    # Avoid division by zero and log of zero by adding a small epsilon
    epsilon = 1e-12
    y_pred = np.clip(y_pred, epsilon, None)
    y_target = np.clip(y_target, epsilon, None)

    # Calculate SID
    term1 = y_pred * np.log(y_pred / y_target)
    term2 = y_target * np.log(y_target / y_pred)
    
    return np.sum(term1 + term2)


def SIS(y_pred, y_target):
    # Ensure inputs are numpy arrays
    y_pred = np.asarray(y_pred)
    y_target = np.asarray(y_target)

    sid_calc = SID(y_pred, y_target)

    return 1.0/(1.0 + sid_calc)
 



# Process the jdx spectrum to fit xmin, xmax and dx
def process_exp_spectrum(spectral_data, dx, xmin=None, xmax=None):
    # Step 1: Round X values to nearest integer but keep them as floats in the data
    rounded_spectrum = [(round(x), y) for x, y in spectral_data]

    # Step 2: Recalculate x_values and y_values
    x_values = [x for x, y in rounded_spectrum]
    y_values = [y for x, y in rounded_spectrum]


    # Step 3: Check if interpolation is needed
    current_deltax = np.mean(np.diff(x_values))
    if dx >= current_deltax:
        print(f"No interpolation needed, current deltax is {current_deltax}, which is smaller than or equal to desired dx {dx}.")
        interpolated_spectrum = rounded_spectrum
    else:
        # Interpolate data to match the desired dx
        # Generate new X values with the desired dx, using np.linspace to include the last value
        new_x_values = np.linspace(x_values[0], x_values[-1], int((x_values[-1] - x_values[0]) / dx) + 1)

        # Interpolate Y values for the new X values
        new_y_values = np.interp(new_x_values, x_values, y_values)

        # Combine new X and Y values into the result
        interpolated_spectrum = list(zip(new_x_values, new_y_values))

    # Step 4: Handle xmin interpolation if needed
    x_values_interpolated = [x for x, y in interpolated_spectrum]
    if xmin is not None and xmin < x_values_interpolated[0]:
        # Add a point at (xmin, 0) and interpolate between xmin and current min x using an exponential function
        x_min = x_values_interpolated[0]-1
        y_min = interpolated_spectrum[0][1]

        # Generate X values between xmin and the current minimum x
        extended_x_values = np.linspace(xmin, x_min, int((x_min - xmin) / dx) + 1)
  
        # Apply an exponential decay function for new points
        decay_constant = 1.0 / 100.0
        extended_y_values = y_min * np.exp(-decay_constant * (x_min - extended_x_values)**2)

        # Combine the extended values with the interpolated spectrum
        extended_spectrum = list(zip(extended_x_values, extended_y_values))
        interpolated_spectrum = extended_spectrum + interpolated_spectrum

    # Step 5: Handle xmax interpolation if needed
    x_values_interpolated = [x for x, y in interpolated_spectrum]
    if xmax is not None and xmax > x_values_interpolated[-1]:
        # Add a point at (xmax, 0) and interpolate between current max x and xmax
        x_max = x_values_interpolated[-1]+1
        y_max = interpolated_spectrum[-1][1]

        # Generate X values between current max x and xmax
        extended_x_values = np.linspace(x_max, xmax, int((xmax - x_max) / dx) + 1)
        xtmp = extended_x_values[0]

        # Linearly interpolate between the current max x and xmax
        decay_constant = 1.0 / 100.0
        extended_y_values = y_max * np.exp(-decay_constant * (extended_x_values - xtmp)**2)

        # Combine the extended values with the interpolated spectrum
        extended_spectrum = list(zip(extended_x_values, extended_y_values))
        interpolated_spectrum = interpolated_spectrum + extended_spectrum

    # Step 6: Recalculate x_values and y_values after xmin/xmax handling
    x_values_interpolated = [x for x, y in interpolated_spectrum]
    y_values_interpolated = [y for x, y in interpolated_spectrum]

    # Step 7: Filter out points outside of the range [xmin, xmax]
    if xmin is not None:
        interpolated_spectrum = [(x, y) for x, y in interpolated_spectrum if x >= xmin]

    if xmax is not None:
        interpolated_spectrum = [(x, y) for x, y in interpolated_spectrum if x <= xmax]

    return interpolated_spectrum
 

def get_symnum(sym):
    sym = sym.lower()  # Convert to lowercase
    symnum = 1 # Assume C1 as standard
    if 'c2' in sym or 's4' in sym:
        symnum = 2
    elif 'c3' in sym or 's6' in sym:
        symnum = 3
    elif 'c4' in sym or 's8' in sym:
        symnum = 4
    elif 'c5' in sym:
        symnum = 5
    elif 'c6' in sym:
        symnum = 6
    elif 'c7' in sym:
        symnum = 7
    elif 'c8' in sym:
        symnum = 8
    elif 'c9' in sym:
        symnum = 9
    elif 'd2' in sym:
        symnum = 4
    elif 'd3' in sym:
        symnum = 6
    elif 'd4' in sym:
        symnum = 8
    elif 'd5' in sym:
        symnum = 10
    elif 'd6' in sym:
        symnum = 12
    elif 'd7' in sym:
        symnum = 14
    elif 'd8' in sym:
        symnum = 16
    elif 'd9' in sym:
        symnum = 18
    elif 't' in sym or 'td' in sym or 'th' in sym:
        symnum = 12
    elif 'o' in sym or 'oh' in sym:
        symnum = 24
    elif 'ih' in sym:
        symnum = 60
    return symnum
  

def print_thermo_setup(num_frequencies: int, num_imag: int, symmetry: str,
    rotational_number: int, scaling_factor: float, rotor_cutoff: float,
    imaginary_cutoff: float, temperature: float):
    # Header and footer for the box
    top_border = "   ..................................................."
    title = "   :                  THERMO SETUP                   :"
    separator = "   :.................................................:"

    # Body with formatted entries
    body = [
        f"   :  # frequencies                          {num_frequencies:<6}  :",
        f"   :  # imaginary                            {num_imag:<6}  :",
        f"   :  temperature                        {temperature:.2f} K    :",
        f"   :  symmetry                               {symmetry:<6}  :",
        f"   :  rotational number                       {rotational_number:<6} :",
        f"   :  scaling factor                  {scaling_factor:<10.7f}     :",
        f"   :  rotor cutoff                   {rotor_cutoff:<10.7f} cm⁻¹ :",
        f"   :  imag. cutoff                  {imaginary_cutoff:<10.7f} cm⁻¹ :",
    ]
    # Assemble the printout
    printout = "\n".join([top_border, title, separator] + body + [separator])
    print(printout+"\n")



