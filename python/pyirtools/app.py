# pyirtools/app.py

import argparse
import pyirtools.calculator

def main():
    parser = argparse.ArgumentParser(description="Python CLI tool for pyirtools")
    parser.add_argument("input", help="Input file")
    parser.add_argument("output", help="Output file")
    parser.add_argument("--verbose", action="store_true", help="Increase output verbosity")
    
    args = parser.parse_args()
    
    #
    # TODO
    #
   
    if args.verbose:
        print(f"Processing input file: {args.input}")
    
    if args.verbose:
        print(f"Output written to: {args.output}")

if __name__ == "__main__":
    main()
