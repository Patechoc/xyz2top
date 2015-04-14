#!/usr/bin/env python

import sys, os
import argparse
import numpy as np
import xyz
import math
import atomInMolecule as mol
import elements

class atomEntity(object):
    def __init__(self, atomObject, atomIndex):
        self.atomIndex     = atomIndex
        self.atomObject    = atomObject
        self.indicesNeighbours = []
    def add_neighbourAtom(self, atomEntityObject):
        self.indicesNeighbours.append(atomEntityObject.atomIndex)
    def __str__(self):
        str = "ATOM ENTITY: index={}, infos= {}".format(self.atomIndex,self.atomObject)
        if len(self.indicesNeighbours) > 0:
            str += ", #neighbours= {} ({})".format(len(self.indicesNeighbours),self.indicesNeighbours)
        else:
            str += ", #neighbours= 0"
        return str

class atomPair(object):
    def __init__(self, atomEntity_i, atomEntity_j):
        self.atomEntity_i = atomEntity_i
        self.atomEntity_j = atomEntity_j
        # INTER-ATOMIC DISTANCE
        self.distance = get_atomic_separation(self.atomEntity_i.atomObject,
                                            self.atomEntity_j.atomObject)
        # COVALENT BOND DISTANCE
        # http://chemwiki.ucdavis.edu/Theoretical_Chemistry/Chemical_Bonding/General_Principles/Covalent_Bond_Distance,_Radius_and_van_der_Waals_Radius
        self.covDist = self.sum_covalent_radii()
    def sum_covalent_radii(self):
        ele_ai = elements.ELEMENTS[self.atomEntity_i.atomObject.atomSymbol]
        ele_aj = elements.ELEMENTS[self.atomEntity_j.atomObject.atomSymbol]
        covRad_ai = ele_ai.covrad
        covRad_aj = ele_aj.covrad
        return covRad_ai + covRad_aj
    def __str__(self):
        return "ATOM PAIR between:\n\
        \t {}\n\
        \t {}\n\
        \t Covalent bond distance: {}\n\
        \t Interatomic distance : {}\
        ".format(self.atomEntity_i,  self.atomEntity_j,
                 self.covDist, self.distance)


class topology(object):
    def __init__(self, molecule, covRadFactor=1.3):
        self.molecule = molecule
        # inter-atomic distance as an upper triangular matrix of atom-pairs
        self.distanceMatrix = None
        # connections detected as covalent bonds
        self.covRadFactor = covRadFactor
        self.atomicPairs = [] # contains all atomPairs
        self.connections = [] # contains only atomPairs detected as connected
        self.atomNeighbours = [atomEntity(ai,i) for i,ai in enumerate(self.molecule.listAtoms)]
    def is_connected(self, pair):
        isConnected = False
        if ( get_atomic_separation(pair.atomEntity_i.atomObject,
                                   pair.atomEntity_j.atomObject)
             < self.covRadFactor * pair.covDist):
            isConnected = True
        return isConnected

    def check_connection(self, ai, aj, i, j):
        entity_i = atomEntity(ai, i)
        entity_j = atomEntity(aj, j)
        pair = atomPair(entity_i, entity_j)
        # add the pair to the list of atomic pairs
        self.atomicPairs.append(pair)
        # if the atoms are 'close' enough, add the pair to the list of connections too
        if self.is_connected(pair):
            self.connections.append(pair)
            self.atomNeighbours[i].add_neighbourAtom(entity_j)
            self.atomNeighbours[j].add_neighbourAtom(entity_i)
            print atomPair(self.atomNeighbours[i], self.atomNeighbours[j])
        #print pair

    def get_atomicConnections(self):
        print "covRadFactor: ",self.covRadFactor
        # go through all unique pairs of atoms
        # compare distance to covalent bond distance (scaled)
        [[self.check_connection(ai, aj, i, j) for j, aj in enumerate(self.molecule.listAtoms) if j>i] for i, ai in enumerate(self.molecule.listAtoms[:-1])]
        print "Nb. of total unique pair of atoms: ",len(self.atomicPairs)
        print "Nb. of covalent bond detected: ",len(self.connections)

    def get_distanceMatrix(self):
        mat = [[get_atomic_separation(ai, aj) for j, aj in enumerate(self.molecule.listAtoms) if j>i] for i, ai in enumerate(self.molecule.listAtoms[:-1])]
        self.distanceMatrix = mat
        return mat

    def read_distance(self, index_i, index_j):
        return self.distanceMatrix[index_i][index_j]

    def build_topology(self):
        self.get_distanceMatrix() # DEPRECATED
        self.get_atomicConnections()

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

def get_atomic_separation(ai,aj):
    """
    aj and aj are atomInfos objects and this function
    returns the interatomic distance which separates them
    """
    return  math.sqrt((ai.xCoord-aj.xCoord)**2
                         +(ai.yCoord-aj.yCoord)**2
                         +(ai.zCoord-aj.zCoord)**2)

def main():
    # read inputs
    args = read_arguments()
    path_to_file = os.path.abspath(args.filename)
    if (args.covRadFactor == None):
        print "no factor for bond distance specified\n\t\
        >> default covalent radius factor will apply.\nRun './main.py --help' for more options."
    else:
        print "Covalent radius factor set to ", args.covRadFactor
        if args.verbose:
            print "Approximate the molecular topology stored in {} \n  \
            with connections detected as covalent bonds if pair-atomic \
            distance goes below {} times the sum of the covalent radii.\
            ".format(args.filename, args.covRadFactor)
    # parse_molecule_XYZ()
    molecule = xyz.parse_XYZ(path_to_file)
    #print molecule

    # build_matrix_distance()
    if (args.covRadFactor != None):
        molecular_topology = topology(molecule, args.covRadFactor)
    else:
        molecular_topology = topology(molecule)
    molecular_topology.build_topology()

    # detect_covalent_bonds() # build unique connected pairs and add connectedAtoms to each atom 
    # get_angles()
    # get_dihedral_angles()
    # print_topology()

if __name__ == "__main__":
    main()
