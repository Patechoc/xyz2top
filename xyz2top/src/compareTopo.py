#!/usr/bin/env python

import sys, os
import argparse
#import numpy as np
#import xyz
#import atomsinmolecule as mol
import topology as topo
import math


class compare_topologies(object):
    def __init__(self, molecule1, molecule2, covRadFactor=1.3):
        comparison_requirements(molecule1, molecule2)
        self.topology1 = topo.topology(molecule1, covRadFactor)
        self.topology2 = topo.topology(molecule2, covRadFactor)

        # self.atomEntities = [atomEntity(ai,i) for i,ai in enumerate(self.molecule.listAtoms)]
        # self.atomicPairs = [] # contains all atomPairs
        # self.covalentBonds = [] # contains only atomPairs detected as connected
        # self.covalentBondAngles = []
        # self.covalentDihedralAngles = []
        # self.covBonds_built = False
        # self.covBondAngles_built = False
        # self.covBondDihedrals_built = False
        # self.build_topology()

    def get_object(self):
        obj = {}
        obj["molecule1"] = self.molecule1.get_object()
        obj["molecule2"] = self.molecule2.get_object()
        # obj["atomEntities"] = [e.get_object() for e in self.atomEntities]
        # obj["atomicPairs"] = [p.get_object() for p in self.atomicPairs]
        # obj["covalentBonds"] = [b.get_object() for b in self.covalentBonds]
        # obj["covalentBondAngles"] = [b.get_object() for b in self.covalentBondAngles]
        # obj["covalentDihedralAngles"] = [b.get_object() for b in self.covalentDihedralAngles]
        return obj

    def __str__(self):
        return "COMPARISON OF TOPOLOGIES (summary):\
        \n\tmolecule: {} ({} atoms)\
        \n\tCovalent radius factor: {}\
        \n\tTotal nb. of possible atomic pairs : {}\
        \n\tTotal nb. of pairs detected as bonds: {}\
        \n\tTotal nb. of angles between bonds: {}\
        \n\tTotal nb. of dihedral angles: {}\
        ".format(self.molecule.shortname,
                 self.molecule.nbAtomsInMolecule,
                 self.covRadFactor)

    def get_as_JSON(self):
        topoComparison = self.get_object()
        import json
        return json.dumps(topo, sort_keys=True, indent=4)

def comparison_requirements(molecule1, molecule2):
    ## the molecules should have the same atoms, provided in the same order
    if self.molecule1.nbAtomsInMolecule != self.molecule2.nbAtomsInMolecule:
        msg =  "Not the same number of atoms comparing {} and {}".format(self.molecule1.shortname, self.molecule2.shortname)
        sys.exit(msg)
    if self.molecule1.charge != self.molecule2.charge:
        msg = "Not the same molecular charge comparing {} and {}".format(self.molecule1.shortname, self.molecule2.shortname)
    sys.exit(msg)
    for atom1, atom2 in zip(self.molecule1.listAtoms, self.molecule2.listAtoms):
        if self.atom1.atomSymbol != self.atom2.atomSymbol:
            msg = "Not the same atom symbols: comparing {} and {}".format(str(atom1), str(atom2))
            sys.exit(msg)
        if self.atom1.atomCharge != self.atom2.atomCharge:
            msg = "Not the same atom charge: comparing {} and {}".format(str(atom1), str(atom2))
            sys.exit(msg)
        if self.atom1.unitDistance != self.atom2.unitDistance:
            msg = "Not the same atom unitDistance: comparing {} and {}".format(str(atom1), str(atom2))
            sys.exit(msg)


def read_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("mol1.xyz",
                        help="First molecular geometry in .XYZ format.")
    parser.add_argument("mol2.xyz",
                        help="Second molecular geometry in .XYZ format.")
    parser.add_argument('-out', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="optional output filename,\
                        if not, default is mol1_vs_mol2.top")
    parser.add_argument("-crf", "--covRadFactor", type=float,
                        help="optional covalent radius factor,\
                        equal to 1 by default")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    args = parser.parse_args()
    return args


# def compare_topologies(filepath1, prefix1, filepath2, prefix2):
#         covRadFactor = -1.
#         list_pairs = []
#         list_triples = []
#         list_quads = []
#         # check that the topologies used the same configuration (= same covRadFactor)
#         config_topo = {"covRadFactor":covRadFactor}
#         # compare the number of covalent bonds
#         # for identical covalent bonds, provides stats (ErrorMaxAbs, ErrorMean, ErrorStd, ErrorRMS)
#         error_bonds ={}

#         # compare the number of angles between covalent bonds
#         # for identical angles btw bonds, provides stats (ErrorMaxAbs, ErrorMean, ErrorStd, ErrorRMS)
#         error_angles ={}

#         # compare the number of dihedral angles between 3 covalent bonds
#         # for identical dihedrals btw 3 bonds, provides stats (ErrorMaxAbs, ErrorMean, ErrorStd, ErrorRMS)
#         error_dihedrals ={}

#         errors = {"config_topo":config_topo,
#                   "error_bonds":error_bonds,
#                   "error_angles":error_angles,
#                   "error_dihedrals":error_dihedrals}
#         return errors

def main():
    # read inputs
    args = read_arguments()
#     path_to_file = os.path.abspath(args.filename)
#     if (args.covRadFactor == None):
#         print "no factor for bond distance specified\n>> default covalent radius factor will apply.\n(Run './main.py --help' for more options.)"
#     else:
#         print "Covalent radius factor set to ", args.covRadFactor
#         if args.verbose:
#             print "Approximate the molecular topology stored in {} \n  \
#             with connections detected as covalent bonds if pair-atomic \
#             distance goes below {} times the sum of the covalent radii.\
#             ".format(args.filename, args.covRadFactor)
#     ### parse_molecule_XYZ()
#     molecule = xyz.parse_XYZ(path_to_file)
#     #print molecule.get_object()
#     #print molecule

#     ### compute the topology
#     if (args.covRadFactor != None):
#         molecular_topology = topology(molecule, args.covRadFactor)
#     else:
#         molecular_topology = topology(molecule)
#     molecular_topology.build_topology()
# #    print molecular_topology.get_as_JSON()
#     print molecular_topology

#     ### print topology to file
#     jsonString = molecular_topology.get_as_JSON()
#     with open('./topology.json', 'w') as outfile:
#         outfile.write(jsonString)

#     print "\nZmatrix format: (not done yet)"
#     print molecular_topology.get_as_Zmatrix()

#     print "\nZmatrix format with variables: (not done yet)"
#     print molecular_topology.get_as_Zmatrix(useVariables=True)

#     print  "\nCreate 3 topology files for bonds, angles and dihedrals + config.txt"
#     [filename_config, filename_bonds, filename_angles, filename_dihedralAngles] = molecular_topology.write_topology_files()
#     print "files generated:" \
#         + "\n\t- " + filename_config \
#         + "\n\t- " + filename_bonds \
#         + "\n\t- " + filename_angles \
#         + "\n\t- " + filename_dihedralAngles;


if __name__ == "__main__":
    main()
