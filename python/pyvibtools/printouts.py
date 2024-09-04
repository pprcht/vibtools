# pyvibtools/printouts

from .calculator import vibtoolsCalculator

def print_vibspectrum(calc: vibtoolsCalculator, outfile=None):
    """"
    Process the vibspectrum from a given vibtoolsCalculator

    Parameters:
    - calc: The vibtoolsCalculator
    - outfile: the vibspectrum file to which to output is written
               If it is None, the spectrum si instead printed to stdout
    """
    # Ensure the calculation has been performed, otherwise don't execute the routine
    if calc.freq is None or calc.intens is None:
        return

    # Prepare the spectrum data as a list of formatted strings
    spectrum_lines = []
    if outfile:
       spectrum_lines.append("$vibrational spectrum")
    if calc.fscal != 1.0:
       spectrum_lines.append(f"# WARNING: frequencies are scaled by factor {calc.fscal:.5f}")
    spectrum_lines.append("#  mode     symmetry     wave number   IR intensity    selection rules")
    spectrum_lines.append("#                         cm**(-1)        km/mol         IR     RAMAN")
    for i, (freq, intens) in enumerate(zip(calc.freq, calc.intens), start=1):
        if abs(freq) < 1e-3: # 6 (bzw. 5) modes should be zero 
            symmetry = ""
            ir = " - "  # Both columns show '-' if both values are small 
        else:
            symmetry = "a"  # Example symmetry label, maybe implement at later date
            if abs(intens) > calc.IR_CUTOFF: # Determine IR activity with cutoff from calculator
               ir = "YES"
            else:
               ir = "no "  
        raman = " - " # Raman activity not implemented at the moment

        line = f"{i:6} {symmetry:>8} {freq:18.2f} {intens:15.5f} {ir:>9} {raman:>7}"
        spectrum_lines.append(line)
    if outfile: 
       spectrum_lines.append("$end")

    # Output to a file if outfile is provided, otherwise print to stdout
    if outfile:
        with open(outfile, 'w') as f:
            for line in spectrum_lines:
                f.write(line + "\n")
        print(f"Vibrational spectrum written to {outfile}")
    else:
        for line in spectrum_lines:
            print(line)
