# pyirtools/app.py

import argparse
import pyirtools.calculator

def main():
    parser = argparse.ArgumentParser(
        description="A program to process vibrational spectra."
    )
    parser.add_argument('-i', '--file', type=str, metavar='FILE', default=None,
                        help=('Specify input file. Depending on the application, '
                             'this file can be variety of file types, e.g. a molecular '
                             'structure or a vibspectrum file.' ))
    parser.add_argument('-hess','--hessian', type=str, metavar='FILE', default=None,
                        help=('Specify non-massweighted Hessian input file. '
                              'The Hessian must be in atomic units (Hartree, Bohr)'))
    parser.add_argument('-dip','--dipole-gradient', type=str, metavar='FILE', default=None,
                        help=('Specify Cartesian dipole derivative input file. '
                              'Dipole gradient must be in atomic units (charge, Bohr)'))
#    parser.add_argument('-m', '--mass', type=str, metavar='<str>',
#                        help=('Specify input file containing user-defined atomic masses, '
#                              'which can be used for atomic mass scaling as discussed in '
#                              'https://dx.doi.org/10.1021/acs.jctc.0c00877. '
#                              'Needs two columns: the atomic number and associated mass.'))
    parser.add_argument('-s', '--scal', type=float, metavar='<float>', default=None,
                        help='Specify a linear frequency scaling factor')
    parser.add_argument('-o', '--output', type=str, metavar='FILE',
                        help=("Specify the output file (if not provided, the output is simply "
                              "written to 'vibspectrum')"))
    parser.add_argument('--plot', action='store_true',
                        help=("Apply Lorentzian line shapes to a the computed/read spectrum"
                              "and plot via matplotlib."))
    parser.add_argument('-msc','--matchscore', nargs=2, metavar=('FILE1', 'FILE2'),
                        help="Specify two vibrational spectrum files for matchscore calculation")

    # Parse arguments
    args = parser.parse_args()

    #
    # TODO
    #

if __name__ == "__main__":
    main()
