#!/usr/bin/env python

import sys, os
import argparse
import numpy as np
#import atomsinmolecule as mol
import topology as topo
import math
import pandas as pd

class topologyDiff(object):
    def __init__(self, molecule1, molecule2, covRadFactor=1.3):
        errors = {}
        requirements_for_comparison(molecule1, molecule2)
        self.molecule1 = molecule1
        self.molecule2 = molecule2
        self.topology1 = topo.topology(molecule1, covRadFactor)
        self.topology2 = topo.topology(molecule2, covRadFactor)
        self.orderedBonds1  = self.topology1.order_convalentBondDistances()
        self.orderedBonds2  = self.topology2.order_convalentBondDistances()
        #print "\n".join([str(elem) for elem in self.orderedBonds2])
        self.orderedAngles1 = self.topology1.order_angles()
        self.orderedAngles2 = self.topology2.order_angles()
        self.orderedDihedral1 = self.topology1.order_dihedralAngles()
        self.orderedDihedral2 = self.topology2.order_dihedralAngles()
        self.error_bonds  = self.compare_bonds(percentLargest = 1.5)
        self.error_angles = self.compare_angles()
        self.error_dihedrals = self.compare_dihedralAngles()
        # print "error_bonds", self.error_bonds
        # print "error_angles", self.error_angles
        # print "error_dihedrals", self.error_dihedrals

    def compare_bonds(self, percentLargest = -1):
        ## Keep all data toghether and filter/sort on it
        nameCol_i      = "index_i"
        nameCol_j      = "index_j"
        nameCol_IDs    = "uniquePairID"
        nameCol_dist   = "distance [A]"
        nameCol_dist1  = "Mol1 dist. [A]"
        nameCol_dist2  = "Mol2 dist. [A]"
        nameCol_errors = "Dist.error [A]"
        nameCol_absError  = "absError [A]"

        # same nb. of bonds?
        if len(self.orderedBonds1) != len(self.orderedBonds2):
            msg =  "Not as many covalents bonds detected in both structures:\n - {}".format(molecule1.shortname, molecule2.shortname)
            sys.exit(msg)
        ## error in distance (Angstrom)  for each bond
        ## checking that the unique ID is the same, if not as many bonds, exit with an error
        id1 = np.array(self.orderedBonds1[1:])[:,0]
        id2 = np.array(self.orderedBonds2[1:])[:,0]
        diffIDs = np.sum(np.absolute(np.subtract(id1, id2)))
        if diffIDs > 0:
            msg =  "As many covalents bonds detected, but not between the same atoms comparing structures:\n - {}".format(molecule1.shortname, molecule2.shortname)
            sys.exit(msg)
        ## Pandas Dataframe
        df1 = pd.DataFrame(self.orderedBonds1[1:], columns=self.orderedBonds1[0])
        df2 = pd.DataFrame(self.orderedBonds2[1:], columns=self.orderedBonds2[0])
        ## convert string to float/int
        for header in [nameCol_dist]:
            df1[header] = df1[header].astype('float64')
            df2[header] = df2[header].astype('float64')
        for header in [nameCol_IDs, nameCol_i, nameCol_j]:
            df1[header] = df1[header].astype('int')
            df2[header] = df2[header].astype('int')
        df1 = df1.rename(columns={nameCol_dist:nameCol_dist1})
        df2 = df2.rename(columns={nameCol_dist:nameCol_dist2})
        df = df1
        df[nameCol_dist2] = df2[nameCol_dist2]
        df[nameCol_errors] = df[nameCol_dist1] - df[nameCol_dist2]
        ###df = df.sort([nameCol_errors, nameCol_IDs], ascending=[False,True])
        df[nameCol_absError] = df[nameCol_errors].abs()
        df = df.sort([nameCol_absError], ascending=[False])
        # print df
        ## STATISTICS
        return get_statistics(df, nameCol_errors, unit="angstrom")

    def compare_angles(self):
        ## Keep all data toghether and filter/sort on it
        nameCol_IDs    = "uniqueID"
        nameCol_i      = "index_i"
        nameCol_j      = "index_j"
        nameCol_k      = "index_k"
        nameCol_anglDeg = 'Angle IJK [deg]'
        nameCol_anglDeg1 = 'Angle1 IJK [deg]'
        nameCol_anglDeg2 = 'Angle2 IJK [deg]'
        nameCol_errors   = "Angle error [deg]"
        nameCol_absError = "absError [deg]"
        nameCol_relError = "relError [deg]"
        # same nb. of angles?
        if len(self.orderedAngles1) != len(self.orderedAngles2):
            msg =  "Not as many covalents angles detected in both structures:\n - {}".format(molecule1.shortname, molecule2.shortname)
            sys.exit(msg)
        ## Pandas Dataframe
        df1 = pd.DataFrame(self.orderedAngles1[1:], columns=self.orderedAngles1[0])
        df2 = pd.DataFrame(self.orderedAngles2[1:], columns=self.orderedAngles2[0])
        ## convert string to float/int
        for header in [nameCol_IDs, nameCol_i, nameCol_j, nameCol_k]:
            df1[header] = df1[header].astype('int')
            df2[header] = df2[header].astype('int')
        for header in [nameCol_anglDeg]:
            df1[header] = df1[header].astype('float64')
            df2[header] = df2[header].astype('float64')
        df1 = df1.rename(columns={nameCol_anglDeg:nameCol_anglDeg1})
        df2 = df2.rename(columns={nameCol_anglDeg:nameCol_anglDeg2})
        df = df1
        df[nameCol_anglDeg2] = df2[nameCol_anglDeg2]
        df[nameCol_errors] = df[nameCol_anglDeg1] - df[nameCol_anglDeg2]
        ## checking that the unique ID is the same, if not as many angles, exit with an error
        diffIDs = pd.DataFrame(df1[nameCol_IDs].values - df2[nameCol_IDs].values).abs().sum()
        if diffIDs.values[0] > 0:
            msg =  "As many covalents angles detected, but not between the same atoms comparing structures:\n - {}".format(molecule1.shortname, molecule2.shortname)
            sys.exit(msg)
        ###df = df.sort([nameCol_errors, nameCol_IDs], ascending=[False,True])
        df[nameCol_absError] = df[nameCol_errors].abs()
        df[nameCol_relError] = df[nameCol_errors].map( lambda x: x if abs(x) < 180. else np.sign(x)*abs(360.-x))
        df = df.sort([nameCol_relError, nameCol_IDs], ascending=[False,True])
        #print pd.DataFrame(d8.values-d7.values)
        ## STATISTICS
        return get_statistics(df, nameCol_relError, unit="degrees")

    def compare_dihedralAngles(self):
        ## Keep all data toghether and filter/sort on it
        nameCol_IDs    = "uniqueID"
        nameCol_i      = "index_i"
        nameCol_j      = "index_j"
        nameCol_k      = "index_k"
        nameCol_l      = "index_l"
        nameCol_dihedDeg  = "Dihedral IJ-KL [deg]"
        nameCol_dihedDeg1 = "Dihedral1 IJ-KL [deg]"
        nameCol_dihedDeg2 = "Dihedral2 IJ-KL [deg]"
        nameCol_errors   = "Dihedral angle error [deg]"
        nameCol_absError = "absError [deg]"
        nameCol_relError = "relError [deg]"
        # same nb. of dihedral angles?
        if len(self.orderedDihedral1) != len(self.orderedDihedral2):
            msg =  "Not as many covalents dihedral angles detected in both structures:\n - {}".format(molecule1.shortname, molecule2.shortname)
            sys.exit(msg)
        ## Pandas Dataframe
        df1 = pd.DataFrame(self.orderedDihedral1[1:], columns=self.orderedDihedral1[0])
        df2 = pd.DataFrame(self.orderedDihedral2[1:], columns=self.orderedDihedral2[0])
        ## convert string to float/int
        for header in [nameCol_IDs, nameCol_i, nameCol_j, nameCol_k, nameCol_l]:
            df1[header] = df1[header].astype('int')
            df2[header] = df2[header].astype('int')
        for header in [nameCol_dihedDeg]:
            df1[header] = df1[header].astype('float64')
            df2[header] = df2[header].astype('float64')
        df1 = df1.rename(columns={nameCol_dihedDeg:nameCol_dihedDeg1})
        df2 = df2.rename(columns={nameCol_dihedDeg:nameCol_dihedDeg2})
        df = df1
        df[nameCol_dihedDeg2] = df2[nameCol_dihedDeg2]
        df[nameCol_errors] = df[nameCol_dihedDeg1] - df[nameCol_dihedDeg2]
        ## checking that the unique ID is the same, if not as many angles, exit with an error
        diffIDs = pd.DataFrame(df1[nameCol_IDs].values - df2[nameCol_IDs].values).abs().sum()
        if diffIDs.values[0] > 0:
            msg =  "As many covalents dihedral angles detected, but not between the same atoms comparin structures:\n - {}".format(molecule1.shortname, molecule2.shortname)
            sys.exit(msg)
        df[nameCol_absError] = df[nameCol_errors].abs()
        df[nameCol_relError] = df[nameCol_errors].map( lambda x: x if abs(x) < 180. else np.sign(x)*(360.-abs(x)))
        df = df.sort([nameCol_relError, nameCol_IDs], ascending=[False,True])
        #print pd.DataFrame(d8.values-d7.values)
        ## STATISTICS
        return get_statistics(df, nameCol_relError, unit="degrees")

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
        \n\tmolecules compared:\
        \n\t\t- {} ({} atoms)\
        \n\t\t- {} ({} atoms)\
        \n\tCovalent radius factor: {}\
        \n\tCovalents bonds errors:\
        \n\t\t- mean: {:-.1e} {}\
        \n\t\t- std:  {:-.1e} {}\
        \n\tCovalents angles errors:\
        \n\t\t- mean: {:-.1e} {}\
        \n\t\t- std:  {:-.1e} {}\
        \n\tDihedral angles errors:\
        \n\t\t- mean: {:-.1e} {}\
        \n\t\t- std:  {:-.1e} {}\
        ".format(self.molecule1.shortname,
                 self.molecule1.nbAtomsInMolecule,
                 self.molecule2.shortname,
                 self.molecule2.nbAtomsInMolecule,
                 self.topology1.covRadFactor,
                 self.error_bonds['mean'],self.error_bonds['unit'],
                 self.error_bonds['stdDev'],self.error_bonds['unit'],
                 self.error_angles['mean'],self.error_angles['unit'],
                 self.error_angles['stdDev'],self.error_angles['unit'],
                 self.error_dihedrals['mean'],self.error_dihedrals['unit'],
                 self.error_dihedrals['stdDev'],self.error_dihedrals['unit'])

    def get_as_JSON(self):
        topoComparison = self.get_object()
        import json
        return json.dumps(topo, sort_keys=True, indent=4)

def requirements_for_comparison(molecule1, molecule2):
    msg = ""
    ## the molecules should have the same atoms, provided in the same order
    if molecule1.nbAtomsInMolecule != molecule2.nbAtomsInMolecule:
        msg =  "Not the same number of atoms comparing:\n-{} and\n-{}".format(molecule1.shortname, molecule2.shortname)
        sys.exit(msg)
    if molecule1.charge != molecule2.charge:
        msg = "Not the same molecular charge comparing:\n-{} and\n-{}".format(molecule1.shortname, molecule2.shortname)
        sys.exit(msg)
    for atom1, atom2 in zip(molecule1.listAtoms, molecule2.listAtoms):
        if atom1.atomSymbol != atom2.atomSymbol:
            msg = "Not the same atom symbols: comparing:\n-{} and\n-{}".format(str(atom1), str(atom2))
            sys.exit(msg)
        if atom1.atomCharge != atom2.atomCharge:
            msg = "Not the same atom charge: comparing:\n-{} and\n-{}".format(str(atom1), str(atom2))
            sys.exit(msg)
        if atom1.unitDistance != atom2.unitDistance:
            msg = "Not the same atom unitDistance: comparing:\n-{} and\n-{}".format(str(atom1), str(atom2))
            sys.exit(msg)


def read_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("file_mol1",
                        help="First molecular geometry in .XYZ format.")
    parser.add_argument("file_mol2",
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


def get_statistics(dataFrame, nameData, unit=""):
    mean     = dataFrame[nameData].mean()
    variance = dataFrame[nameData].var()
    stdDev   = dataFrame[nameData].std()
    mad      = dataFrame[nameData].mad()
    maxAbs   = dataFrame[nameData].abs().max()
    return {
        "unit":unit,
        "mean":mean,
        "variance":variance,
        "stdDev":stdDev,
        "mad":mad,  ## mean/average absolute deviation
        "maxAbs":maxAbs}

def example_valinomycin_pureLinK_vs_LinKwithDF():
    # read inputs
    # args = read_arguments()
    # path_to_file1 = os.path.abspath(args.file_mol1)
    # path_to_file2 = os.path.abspath(args.file_mol2)
    path_to_file1 = "/home/ctcc2/Documents/CODE-DEV/xyz2top/xyz2top/tests/files/valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ.xyz"
    path_to_file2 = "/home/ctcc2/Documents/CODE-DEV/xyz2top/xyz2top/tests/files/valinomycin_geomOpt_DFT-b3lyp-noDF_cc-pVTZ.xyz"
    import xyz2molecule as xyz
    molecule1 = xyz.parse_XYZ(path_to_file1)
    molecule2 = xyz.parse_XYZ(path_to_file2)
    diff = topologyDiff(molecule1, molecule2, covRadFactor=1.3)


if __name__ == "__main__":
    example_valinomycin_pureLinK_vs_LinKwithDF()
