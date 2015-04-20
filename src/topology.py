#!/usr/bin/env python

import sys, os
import argparse
import numpy as np
import xyz
import math
import atomInMolecule as mol
import elements

class atomEntity(object):
    def __init__(self, atomInfos, atomIndex):
        self.atomIndex     = atomIndex
        self.atomInfos    = atomInfos
        self.neighbourIndices = []
    def add_neighbourAtom(self, atomEntityObject):
        self.neighbourIndices.append(atomEntityObject.atomIndex)
    def __str__(self):
        str = "ATOM ENTITY: index={}, infos= {}".format(self.atomIndex,self.atomInfos)
        if len(self.neighbourIndices) > 0:
            str += ", #neighbours= {} ({})".format(len(self.neighbourIndices),self.neighbourIndices)
        else:
            str += ", #neighbours= 0"
        return str

class atomPair(object):
    def __init__(self, atomEntity_i, atomEntity_j):
        self.atomEntity_i = atomEntity_i
        self.atomEntity_j = atomEntity_j
        # INTER-ATOMIC DISTANCE
        self.distance = get_interatomic_distance(self.atomEntity_i.atomInfos,
                                            self.atomEntity_j.atomInfos)
        # COVALENT BOND DISTANCE
        # http://chemwiki.ucdavis.edu/Theoretical_Chemistry/Chemical_Bonding/General_Principles/Covalent_Bond_Distance,_Radius_and_van_der_Waals_Radius
        self.covDist = self.sum_covalent_radii()
    def __str__(self):
        return "ATOM PAIR between:\n\
        \t {}\n\
        \t {}\n\
        \t Covalent bond distance: {}\n\
        \t Interatomic distance : {}\
        ".format(self.atomEntity_i,  self.atomEntity_j,
                 self.covDist, self.distance)
    def sum_covalent_radii(self):
        ele_ai = elements.ELEMENTS[self.atomEntity_i.atomInfos.atomSymbol]
        ele_aj = elements.ELEMENTS[self.atomEntity_j.atomInfos.atomSymbol]
        covRad_ai = ele_ai.covrad
        covRad_aj = ele_aj.covrad
        return covRad_ai + covRad_aj

class atomTriple(object):
    def __init__(self, atomEntity_i, atomEntity_j, atomEntity_k):
        self.atomEntity_i = atomEntity_i
        self.atomEntity_j = atomEntity_j
        self.atomEntity_k = atomEntity_k
        self.vector_ji = np.array(atomEntity_i.atomInfos.coordinates()) - np.array(atomEntity_j.atomInfos.coordinates())
        self.vector_jk = np.array(atomEntity_k.atomInfos.coordinates()) - np.array(atomEntity_j.atomInfos.coordinates())
        self.distance_ji = math.sqrt(np.dot(self.vector_ji, self.vector_ji))
        self.distance_jk = math.sqrt(np.dot(self.vector_jk, self.vector_jk))
        self.cos_angle_ijk = np.dot(self.vector_ji,self.vector_jk)/self.distance_ji/self.distance_jk
        self.angle_ijk = np.arccos(self.cos_angle_ijk)
    def __str__(self):
        return "ATOM TRIPLE between:\n\
        \tI: {}\n\
        \tJ: {}\n\
        \tK: {}\n\
        \tdistance JI: {}\n\
        \tdistance JK: {}\n\
        \tangle IJK=acos(JI.JK): {} radians, {} degres\
        ".format(self.atomEntity_i,  self.atomEntity_j, self.atomEntity_k,
                 self.distance_ji, self.distance_jk,
                 self.angle_ijk, self.get_angle_ijk())
    def get_angle_ijk(self, inDegree=True):
        if inDegree:
            return self.angle_ijk*180./math.pi
        else:
            return self.angle_ijk

class atomQuadruple(object):
    def __init__(self, atomEntity_i, atomEntity_j, atomEntity_k, atomEntity_l):
        self.atomEntity_i = atomEntity_i
        self.atomEntity_j = atomEntity_j
        self.atomEntity_k = atomEntity_k
        self.atomEntity_l = atomEntity_l
        self.vector_ji = np.array(atomEntity_i.atomInfos.coordinates()) - np.array(atomEntity_j.atomInfos.coordinates())
        #print "vecors"
        #print self.vector_ji
        self.vector_jk = np.array(atomEntity_k.atomInfos.coordinates()) - np.array(atomEntity_j.atomInfos.coordinates())
        #print self.vector_jk
        self.vector_kl = np.array(atomEntity_l.atomInfos.coordinates()) - np.array(atomEntity_k.atomInfos.coordinates())
        #print self.vector_jk
        self.distance_ji = self.get_distance(self.vector_ji)
        #print "distances"
        #print self.distance_ji
        self.distance_jk = self.get_distance(self.vector_jk)
        #print self.distance_jk
        self.distance_kl = self.get_distance(self.vector_kl)
        #print self.distance_kl
        self.vector_orthoJK_JI = np.cross(self.vector_jk/self.distance_jk, self.vector_ji/self.distance_ji)
        self.dist_orthoJK_JI = self.get_distance(self.vector_orthoJK_JI)

        self.vector_orthoKJ_IJ = np.cross(-self.vector_jk/self.distance_jk, -self.vector_ji/self.distance_ji)
        self.dist_orthoKJ_IJ = self.get_distance(self.vector_orthoKJ_IJ)
        #print "n1 n2"
        #print self.vector_orthoJK_JI
        #print "norm orthoJK_JI: ",np.dot(self.vector_orthoJK_JI,self.vector_orthoJK_JI)
        self.vector_orthoKJ_KL = np.cross(-self.vector_jk/self.distance_jk, self.vector_kl/self.distance_kl)
        self.dist_orthoKJ_KL = self.get_distance(self.vector_orthoKJ_KL)
        #print self.vector_orthoKJ_KL
        #print "norm orthoKJ_KL: ",np.dot(self.vector_orthoKJ_KL,self.vector_orthoKJ_KL)
        self.cos_angle_ijkl = np.dot(self.vector_orthoJK_JI/self.dist_orthoJK_JI, self.vector_orthoKJ_KL/self.dist_orthoKJ_KL)
        self.sin_angle_ijkl = np.cross(self.vector_orthoJK_JI/self.dist_orthoJK_JI, self.vector_orthoKJ_KL/self.dist_orthoKJ_KL)
        #print "cos="
        #print self.cos_angle_ijkl
        if abs(self.cos_angle_ijkl) > 1.:
            print "|cos(X)|>1 ????"
            sys.exit(1)
        self.angle_ijkl = np.arccos(self.cos_angle_ijkl)
        self.angleM_ijkl = math.asin(math.sqrt(np.dot(self.sin_angle_ijkl,self.sin_angle_ijkl)))
        #print "arccos:"
        #print self.angle_ijkl
        #print self.get_dihedral_angle_ijkl 

        ##  The vectors \mathbf{b}_1 and \mathbf{b}_2 define the first plane,
        ## whereas \mathbf{b}_2 and \mathbf{b}_3 define the second plane.
        self.vector_B1 = -self.vector_ji / self.get_distance(self.vector_ji)
        self.vector_B2 =  self.vector_jk / self.get_distance(self.vector_jk)
        self.vector_B3 =  self.vector_kl / self.get_distance(self.vector_kl)
        # n1= b1 x b2
        self.vector_N1 = np.cross(self.vector_B1, self.vector_B2)
        # n2= b2 x b3
        self.vector_N2 = np.cross(self.vector_B2, self.vector_B3)
        # Dihedral= atan2( dot(n1xn2,b2/|b2|) , dot(n1,n2))
        self.dihedral = np.arctan2( np.dot(np.cross(self.vector_N1, self.vector_N2),self.vector_B2), np.dot(self.vector_N1, self.vector_N2) )
    def get_distance(self, vector):
        return math.sqrt(np.dot(vector, vector))
    def __str__(self):
        return "ATOM QUADRUPLE between:\n\
        \tI: {}\n\
        \tJ: {}\n\
        \tK: {}\n\
        \tL: {}\n\
        \tdistance JI: {}\n\
        \tdistance JK: {}\n\
        \tdistance KL: {}\n\
        \tsin(...) {}\n\
        \tangle IJKL=acos(n1,n2): {} radians, {} degres\n\
        \tangle IJKL=acos(n1,n2): {} radians, {} degres (math module)\n\
        \tdihedral =atan2(y,x): {} radians, {} degres\n\
        \t  where n1 and n2 are vector orthogonal to planes (IJK) and (JKL) respectively\
        ".format(self.atomEntity_i,  self.atomEntity_j, self.atomEntity_k, self.atomEntity_l,
                 self.distance_ji, self.distance_jk, self.distance_kl,
                 self.sin_angle_ijkl,
                 self.angle_ijkl, self.get_dihedral_angle_ijkl(), self.angleM_ijkl, self.angleM_ijkl*180./math.pi,
                 self.dihedral, self.dihedral*180./math.pi)
    def get_dihedral_angle_ijkl(self, inDegree=True):
        if inDegree:
            return self.angle_ijkl*180./math.pi
        else:
            return self.angle_ijkl

class topology(object):
    def __init__(self, molecule, covRadFactor=1.3):
        self.molecule = molecule
        # inter-atomic distance as an upper triangular matrix of atom-pairs
        self.distanceMatrix = None
        # connections detected as covalent bonds
        self.covRadFactor = covRadFactor
        self.atomEntities = [atomEntity(ai,i) for i,ai in enumerate(self.molecule.listAtoms)]
        self.atomicPairs = [] # contains all atomPairs
        self.covalentBonds = [] # contains only atomPairs detected as connected
        self.covalentBondAngles = []
        self.covalentDihedralAngles = []
    def get_indices_neighbouringAtoms(self, indexAtomEntity):
        entity = self.get_atomEntity_by_index(indexAtomEntity)
        return entity.neighbourIndices
    def get_atomEntity_by_index(self, indexAtomEntity):
        return [ai for i,ai in enumerate(self.atomEntities) if i == indexAtomEntity][0]
    def __str__(self):
        return "TOPOLOGY:\
        \n\tmolecule: {} ({} atoms)\
        \n\tCovalent radius factor: {}\
        \n\tTotal nb. of possible atomic pairs : {}\
        \n\tTotal nb. of pairs detected as bonds: {}\
        \n\tTotal nb. of angles between bonds: {}\
        \n\tTotal nb. of dihedral angles: {}\
        ".format(self.molecule.shortname, self.molecule.nbAtomsInMolecule,
                 self.covRadFactor, len(self.atomicPairs),
                 len(self.covalentBonds),
                 len(self.covalentBondAngles),
                 len(self.covalentDihedralAngles))

    def is_connected(self, pair):
        isConnected = False
        if ( get_interatomic_distance(pair.atomEntity_i.atomInfos,
                                   pair.atomEntity_j.atomInfos)
             < self.covRadFactor * pair.covDist):
            isConnected = True
        return isConnected

    def check_covalentBond(self, ai, aj, i, j):
        entity_i = atomEntity(ai, i)
        entity_j = atomEntity(aj, j)
        pair = atomPair(entity_i, entity_j)
        # add the pair to the list of atomic pairs
        self.atomicPairs.append(pair)
        # if the atoms are 'close' enough, add the pair to the list of covalentBonds too
        if self.is_connected(pair):
            self.covalentBonds.append(pair)
            self.atomEntities[i].add_neighbourAtom(entity_j)
            self.atomEntities[j].add_neighbourAtom(entity_i)
            #print atomPair(self.atomEntities[i], self.atomEntities[j])
        #print pair

    def get_covalentBonds(self):
        # go through all unique pairs of atoms
        # compare distance to covalent bond distance (scaled)
        [[self.check_covalentBond(ai, aj, i, j) for j, aj in enumerate(self.molecule.listAtoms) if j>i] for i, ai in enumerate(self.molecule.listAtoms[:-1])]
        #print "Nb. of total unique pair of atoms: ",len(self.atomicPairs)
        #print "Nb. of covalent bond detected: ",len(self.covalentBonds)

    def get_covalentBondAngles(self):
        # reduce the search to the atoms that have 'at least' 2 neighbours
        indicesAtomsWithEnoughBonds =[j for j,aj in enumerate(self.atomEntities)
                                      if len(self.get_atomEntity_by_index(j).neighbourIndices) > 1]
        for j in indicesAtomsWithEnoughBonds:
            for i in self.get_indices_neighbouringAtoms(j):
                for k in self.get_indices_neighbouringAtoms(j):
                    if k>i:
                        ai = self.get_atomEntity_by_index(i)
                        aj = self.get_atomEntity_by_index(j)
                        ak = self.get_atomEntity_by_index(k)
                        self.covalentBondAngles.append(atomTriple(ai,aj,ak))
                        #print  atomTriple(ai,aj,ak)
                        #if aj.atomInfos.atomSymbol == "C" and j==9:
                        #    print atomTriple(ai,aj,ak)

    def get_covalentDihedralAngles(self):
        ### reduce the search to the bonding atoms that have 'at least' 2 neighbours each
        indicesPairsWithEnoughNeighbours = []
        for indPair,pair in enumerate(self.covalentBonds):
            #print "index of the pair: ", indPair
            #print atomPair(self.get_atomEntity_by_index(pair.atomEntity_i.atomIndex),\
            #self.get_atomEntity_by_index(pair.atomEntity_j.atomIndex))
            indicesNeighboursAtom2 = self.get_indices_neighbouringAtoms(pair.atomEntity_i.atomIndex)
            indicesNeighboursAtom3 = self.get_indices_neighbouringAtoms(pair.atomEntity_j.atomIndex)
            if len(indicesNeighboursAtom2)>1 and len(indicesNeighboursAtom3)>1:
                indicesPairsWithEnoughNeighbours.append(indPair)
        ### for each such bond, find all possible plans to compare and get dihedral angles
        for indPair in indicesPairsWithEnoughNeighbours:
            pairJK = self.covalentBonds[indPair]
            j = pairJK.atomEntity_i.atomIndex
            k = pairJK.atomEntity_j.atomIndex
            atomJ  = self.get_atomEntity_by_index(j)
            atomK  = self.get_atomEntity_by_index(k)
            indicesAtom1 = [i for i in self.get_indices_neighbouringAtoms(j) if i!=k]
            indicesAtom4 = [l for l in self.get_indices_neighbouringAtoms(k) if l!=j]
            for i in indicesAtom1:
                atomI = self.get_atomEntity_by_index(i)
                for l in indicesAtom4:
                    atomL = self.get_atomEntity_by_index(l)
                    self.covalentDihedralAngles.append(atomQuadruple(atomI,atomJ,atomK,atomL))
                    #print atomQuadruple(atomI,atomJ,atomK,atomL)
                    if atomJ.atomInfos.atomSymbol == "C" and j==9 and k==19:
                        print atomQuadruple(atomI,atomJ,atomK,atomL)

    def build_topology(self):
        self.get_covalentBonds()
        self.get_covalentBondAngles()
        self.get_covalentDihedralAngles()

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

def get_interatomic_distance(atomInfos_i,atomInfos_j):
    """
    atomInfos_i and atomInfos_j are atomInfos objects and this function
    returns the interatomic distance which separates them
    """
    return  math.sqrt((atomInfos_i.xCoord-atomInfos_j.xCoord)**2
                      +(atomInfos_i.yCoord-atomInfos_j.yCoord)**2
                      +(atomInfos_i.zCoord-atomInfos_j.zCoord)**2)

def main():
    # read inputs
    args = read_arguments()
    path_to_file = os.path.abspath(args.filename)
    if (args.covRadFactor == None):
        print "no factor for bond distance specified\n>> default covalent radius factor will apply.\n(Run './main.py --help' for more options.)"
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
    print molecular_topology
    # detect_covalent_bonds() # build unique connected pairs and add connectedAtoms to each atom 
    # get_angles()
    # get_dihedral_angles()
    # print_topology()

if __name__ == "__main__":
    main()
