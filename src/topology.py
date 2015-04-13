#!/usr/bin/env python

import sys, os
import argparse
import numpy as np
import xyz

def read_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename",
                        help="molecular geometry in .XYZ format")
    parser.add_argument('-out', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="optional output filename,\
                        if not, default is [filename].top")
    parser.add_argument("-crf", "--covRadFactor", type=float,
                        help="optional covalent radius factor,\
                        equal to 1 by default")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    args = parser.parse_args()
    return args



def main():
    # read inputs
    args = read_arguments()
    path_to_file = os.path.abspath(args.filename)
    covRadFactor = 1.
    if (args.covRadFactor == None):
        print "no factor for bond distance specified \
        > default covalent radius factor = {}\
        ".format(covRadFactor)
    else:
        covRadFactor = args.covRadFactor
        print "Covalent radius factor set to ",covRadFactor
        if args.verbose:
            print "Approximate the molecular topology stored in {} \n  \
            with connections detected as covalent bonds if pair-atomic \
            distance goes below {} times the sum of the covalent radii.\
            ".format(args.filename, args.covRadFactor)
    # parse_molecule_XYZ()
    molecule = xyz.parse_XYZ(path_to_file)
    print molecule

    # get_atoms()
    # build_matrix_distance()
    # detect_covalent_bonds() # build unique connected pairs and add connectedAtoms to each atom 
    # get_angles()
    # get_dihedral_angles()
    # print_topology()

if __name__ == "__main__":
    main()
