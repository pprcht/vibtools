# pyvibtools/app.py

import argparse
from pyvibtools.calculator import vibtoolsCalculator, matchscore, calc_autorange
from pyvibtools.printouts import print_vibspectrum, export_vibspectrum_to_csv
from pyvibtools.multiplot import *

########################################################################################
########################################################################################
def main():
########################################################################################
########################################################################################
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
    general_group.add_argument('-o', '--output', type=str, metavar='FILE', default=None,
                        help=("Specify the output file (written in TM vibspectrum format)"))

#######################

    function_group = parser.add_argument_group('Processing/Plotting Options')
    function_group.description = ('\nThese options will be applied to the spectrum'
                                  ' and are useful for additional processing, like plotting'
                                  ' or comparing spectra via match scores.')
    function_group.add_argument('--plot', action='store_true',
                        help=("Apply Lorentzian line shapes to a the computed/read spectrum"
                              "and plot via matplotlib."))
    function_group.add_argument('--multiplot', '-mp', metavar='FILE', type=str, nargs='+',
                        help='One or more files to plot together with -i', required=False)
    function_group.add_argument('-ocsv', action='store_true',
                        help=("Write the spectrum (additionally) to vibspectrum.csv"))
    function_group.add_argument('-msc','--matchscore', type=str, metavar='FILE2', default=None,
                        help="Specify vibrational spectrum file for matchscore calculation "
                             "comparing the computed/read-in (-i) spectrum with FILE2")

#######################

    data_group = parser.add_argument_group('Spectral Data Processing Parameters')
    data_group.description = ('\nThese options are used to define parameters'
                              ' for the spectra processing via the CLI.'
                              ' These options will be primarily applied to the spectrum'
                              ' defined via the -i argument.' )
    data_group.add_argument('-s', '--scal', type=float, metavar='<float>', default=None,
                        help='Specify a linear frequency scaling factor') 
    data_group.add_argument('-xmin', type=float, metavar='<float>', default=None,
                        help=("Specify minimum frequency value considered for spectra processing"))
    data_group.add_argument('-xmax', type=float, metavar='<float>', default=None,
                        help=("Specify maximum frequency value considered for spectra processing"))
    data_group.add_argument('--autox', action='store_true',
                        help='Allow automatic determination of frequency range if at least one spectrum is experimental')
    data_group.add_argument('-dx', type=int, metavar='<int>', default=None,
                        help=("Specify frequency interval considered for spectra processing."
                              " The minimum value is 1 cm⁻¹"))
#    data_group.add_argument('-m', '--mass', type=str, metavar='<str>',
#                        help=('Specify input file containing user-defined atomic masses, '
#                              'which can be used for atomic mass scaling as discussed in '
#                              'https://dx.doi.org/10.1021/acs.jctc.0c00877. '
#                              'Needs two columns: the atomic number and associated mass.'))

########################################################################################
########################################################################################

    # Parse arguments
    args = parser.parse_args()

    # Initialize empty
    vibspec1 = vibtoolsCalculator()

    # Read files
    vibspec1.read(xyzfile=args.file, hessfile=args.hessian, dipfile=args.dipole_gradient)

    # More arguments from argparser that may affect further processing    
    if args.scal:
       vibspec1.fscal = args.scal

    # Calculate and print frequencies:
    if vibspec1.atoms is not None and vibspec1.hessian is not None:
       vibspec1.compute()
       print_vibspectrum(vibspec1)
    else:
       # Try to read in the file as vibspectrum
       vibspec1.read(vibspecfile=args.file)

    # Save the vibspectrum to file (TM format)
    if args.output:
       print_vibspectrum(vibspec1, args.output)
    # Additional output to vibpectrum.csv 
    if args.ocsv:
       export_vibspectrum_to_csv(vibspec1)    

########################################################################################
########################################################################################

    # Bounds, for further applications in the following
    if args.xmin:
       vibspec1.xmin = args.xmin
    if args.xmax:
       vibspec1.xmax = args.xmax
    if args.dx:
       vibspec1.dx = float(args.dx)

    # If requested, plot
    if args.plot and args.multiplot is None:
       vibspec1.plot()
    elif args.multiplot and not args.plot:
       spectra_list = []
       spectra_list.append(vibspec1)
       spectra_labels = []
       spectra_labels.append(vibspec1.filename)
       for file_name in args.multiplot:
          process_file_multiplot(file_name, spectra_list, spectra_labels)
       multiplot(spectra_list, labels=spectra_labels)
    elif args.plot and args.multiplot:
       raise ValueError("Please use --plot or --multiplot, but not both")

    # matchscore calculation
    if args.matchscore:
       file2 = args.matchscore
       vibspec2 = vibtoolsCalculator()
       vibspec2.read(vibspecfile=file2)
       # Check and print bounds
       calc_autorange(vibspec1,vibspec2,args.autox)
       print('%12s : %8.2f' % ('xmin / cm⁻¹', vibspec1.xmin))
       print('%12s : %8.2f' % ('xmax / cm⁻¹', vibspec1.xmax))
       print('%12s : %8.2f' % ('dx / cm⁻¹', vibspec1.dx))
       # Using the two spectra, calculate:
       mscs = matchscore(vibspec1, vibspec2)
       print(f"\nMatchscores:")
       print('%10s %10s %10s %10s' % ('MSC','EUC','PCC','SIS')) 
       print('%10.4f %10.4f %10.4f %10.4f' % (mscs[0], mscs[1], mscs[2], mscs[3])) 


########################################################################################
########################################################################################

if __name__ == "__main__":
    main()

