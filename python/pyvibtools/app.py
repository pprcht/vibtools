# pyvibtools/app.py

import argparse
from pyvibtools.calculator import vibtoolsCalculator, matchscore
from pyvibtools.printouts import print_vibspectrum

def main():
    parser = argparse.ArgumentParser(
        description="A program to process vibrational spectra."
    )

    parser.add_argument('-i', '--file', type=str, metavar='FILE', default=None, required=True,
                        help=('Specify input file. Depending on the application, '
                             'this file can be variety of file types, e.g. a molecular '
                             'structure or a vibspectrum file.' ))

    general_group = parser.add_argument_group('General Options')
    general_group.description = '\nGeneral options for input and output operations'
    general_group.add_argument('-hess','--hessian', type=str, metavar='FILE', default=None,
                        help=('Specify non-massweighted Hessian input file. '
                              'The Hessian must be in atomic units (Hartree, Bohr)'))
    general_group.add_argument('-dip','--dipole-gradient', type=str, metavar='FILE', default=None,
                        help=('Specify Cartesian dipole derivative input file. '
                              'Dipole gradient must be in atomic units (charge, Bohr)'))
#    general_group.add_argument('-m', '--mass', type=str, metavar='<str>',
#                        help=('Specify input file containing user-defined atomic masses, '
#                              'which can be used for atomic mass scaling as discussed in '
#                              'https://dx.doi.org/10.1021/acs.jctc.0c00877. '
#                              'Needs two columns: the atomic number and associated mass.'))
    general_group.add_argument('-s', '--scal', type=float, metavar='<float>', default=None,
                        help='Specify a linear frequency scaling factor')
    general_group.add_argument('-o', '--output', type=str, metavar='FILE', default=None,
                        help=("Specify the output file (written in TM vibspectrum format)"))


    function_group = parser.add_argument_group('Processing/Plotting Options')
    function_group.description = ('\nThese options will be applied to the spectrum'
                                  ' and are useful for additional processing, like plotting'
                                  ' or comparing spectra via match scores')
    function_group.add_argument('--plot', action='store_true',
                        help=("Apply Lorentzian line shapes to a the computed/read spectrum"
                              "and plot via matplotlib."))
    function_group.add_argument('-msc','--matchscore', type=str, metavar='FILE2', default=None,
                        help="Specify vibrational spectrum file for matchscore calculation "
                             "comparing the computed/read-in (-i) spectrum with FILE2")

    # Parse arguments
    args = parser.parse_args()

    # Initialize empty
    ircalc1 = vibtoolsCalculator()

    # Read files
    ircalc1.read(xyzfile=args.file, hessfile=args.hessian, dipfile=args.dipole_gradient)

    # More arguments from argparser that may affect further processing    
    if args.scal:
       ircalc1.fscal = args.scal

    # Calculate and print frequencies:
    if ircalc1.atoms is not None and ircalc1.hessian is not None:
       ircalc1.compute()
       print_vibspectrum(ircalc1)
    else:
       # Try to read in the file as vibspectrum
       ircalc1.read(vibspecfile=args.file)

    # Save the vibspectrum to file (TM format)
    if args.output:
       print_vibspectrum(ircalc1, args.output)

    # If requested, plot
    if args.plot:
       ircalc1.plot()

    # matchscore calculation
    if args.matchscore:
       file2 = args.matchscore
       ircalc2 = vibtoolsCalculator()
       ircalc2.read(vibspecfile=file2)
       mscs = matchscore(ircalc1, ircalc2)
       print(f"\nMatchscores:")
       print('%10s %10s %10s' % ('MSC','EUC','PCC')) 
       print('%10.4f %10.4f %10.4f' % (mscs[0], mscs[1], mscs[2])) 


if __name__ == "__main__":
    main()
