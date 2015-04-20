#!/usr/bin/env python

import sys, os
import re
import numpy.testing as npt


Atomic_NUMBERS = {"H": 1,"He": 2,
                  "Li": 3,"Be": 4,"B": 5,"C": 6,"N": 7,"O": 8,"F": 9,"Ne": 10,
                  "Na": 11,"Mg": 12,"Al": 13,"Si": 14,"P": 15,"S": 16,"Cl": 17,"Ar": 18,
                  "K": 19,"Ca": 20,"Sc": 21,"Ti": 22,"V": 23,"Cr": 24,"Mn": 25,"Fe": 26,"Co": 27,"Ni": 28,"Cu": 29,"Zn": 30,"Ga": 31,"Ge": 32,"As": 33,"Se": 34,"Br": 35,"Kr": 36,
                  "Rb": 37,"Sr": 38,"Y": 39,"Zr": 40,"Nb": 41,"Mo": 42,"Tc": 43,"Ru": 44,"Rh": 45,"Pd": 46,"Ag": 47,"Cd": 48,"In": 49,"Sn": 50,"Sb": 51,"Te": 52,"I": 53,"Xe": 54,
                  "Cs": 55,"Ba": 56,"La": 57,"Hf": 72,"Ta": 73,"W": 74,"Re": 75,"Os": 76,"Ir": 77,"Pt": 78,"Au": 79,"Hg": 80,"Tl": 81,"Pb": 82,"Bi": 83,"Po": 84,"At": 85,"Rn": 86,
                  "Fr": 87,"Ra": 88,"Ac": 89,"Ku": 104,"Ha": 105}


class atomInfos(object):
    def __init__(self, atomSymbol="", atomCharge=None):
        self.bohr_in_angstrom = 0.52917721092   #0.5291772083
        self.atomSymbol = atomSymbol
        if atomCharge is None:
            self.atomCharge = float(Atomic_NUMBERS[atomSymbol.strip()])
        else:
            self.atomCharge = float(atomCharge)
        self.unitDistance = None
        self.xCoord = None
        self.yCoord = None
        self.zCoord = None
    def get_object(self):
        json = {}
        json["atomSymbol"] = self.atomSymbol
        json["atomCharge"] = self.atomCharge
        json["bohr_in_angstrom"] = self.bohr_in_angstrom
        json["unitDistance"] = self.unitDistance
        json["coordinates"] = self.coordinates()
        return json
    def __str__(self):
        return self.get_content_atomCoord()
    def setSymbol(self, symbol):
        self.atomSymbol = symbol
    def coordinates(self):
        return [self.xCoord,self.yCoord,self.zCoord]
    def setUnitDistance(self, unitDistance):
        self.unitDistance = unitDistance
    def setCharge(self, charge):
        self.atomCharge = float(charge)
    def setAtomCoord(self, x,y,z):
        self.xCoord = float(x)
        self.yCoord = float(y)
        self.zCoord = float(z)
    def get_content_atomCoord(self, printUnit=True, newUnitDistance=None):
        assert self.unitDistance != None, "The unit distance of an atom should be known to be printed"
        s = ""
        if newUnitDistance != None:
            # atom coordinates already have the expected unit distance
            if  re.match(newUnitDistance, str(self.unitDistance), re.I):
                s += get_content_atomCoord(self, printUnit)
            # atom coordinates are in Angstrom and expected to be shown in Bohr
            elif re.match("bohr", newUnitDistance, re.I) and re.match("angstrom", self.unitDistance, re.I):
                #1 bohr=0.5291772083 Angstrom
                s += '{0:s}     {1: 7.14f}     {2: 7.14f}     {3: 7.14f}'.format(self.atomSymbol, self.xCoord/self.bohr_in_angstrom, self.yCoord/self.bohr_in_angstrom, self.zCoord/self.bohr_in_angstrom)
            # atom coordinates are in Bohr and expected to be shown in Angstrom
            elif re.match("angstrom", newUnitDistance, re.I) and re.match("bohr", str(self.unitDistance), re.I):
                s += '{0:s}     {1: 7.14f}     {2: 7.14f}     {3: 7.14f}'.format(self.atomSymbol, self.xCoord*self.bohr_in_angstrom, self.yCoord*self.bohr_in_angstrom, self.zCoord*self.bohr_in_angstrom)
        else:
            s += '{0:s}     {1: 7.14f}     {2: 7.14f}     {3: 7.14f}'.format(self.atomSymbol, self.xCoord, self.yCoord, self.zCoord)
        if printUnit:
            s += ' (in {0:s})'.format(str(self.unitDistance))
        return s

    def print_atomCoord(self):
        print self.get_content_atomCoord()

class groupSameAtoms(atomInfos):
    def __init__(atomSymbol, atomCharge=None):
        atomInfos.__init__(self, atomSymbol, atomCharge=None)
        self.nbAtomsInGroup = 0
        self.listAtomsCoord = []
    def addAtomInfo(self,atom):
        assert isinstance(atom,atomInfos), 'Trying to add something which is not an atomInfos object to a groupSameAtoms object'
        npt.assertEqual(self.atomsSymbol == atom.atomSymbol,
                        err_msg="Adding an atom to a group with not the same symbol")
        npt.assertEqual(self.atomCharge == atom.atomCharge,
                        err_msg="Adding an atom to a group with not the same charge")
        self.nbAtomsInGroup += 1
        self.listAtomsCoord.append(atom)
    def setAtomCharge(self, charge):
        self.atomCharge = float(charge)
    def getContent_DALTON_groupSameAtoms(self):
        s = ''
        s += 'Charge={0} Atoms={1:d}\n'.format(self.atomCharge, self.nbAtomsInGroup)
        for a in self.listAtomsCoord:
            s += a.get_content_atomCoord()
        return s
    def print_DALTON_groupSameAtoms(self):
        print '{0}'.format(self.getContent_DALTON_groupSameAtoms())


class molecule(object):
    def __init__(self, shortname="", name="", comments=""):
        self.name = name
        self.shortname = shortname
        if (name.strip() == ""):
            self.name = shortname
        self.nbAtomsInMolecule = 0
        self.unitDistance = None
        self.charge= None
        self.comments = comments
        self.listAtoms = []
    def get_object(self):
        json = {}
        json["name"] = self.name
        json["shortname"] = self.shortname
        json["nbAtomsInMolecule"] = self.nbAtomsInMolecule
        json["unitDistance"] = self.unitDistance
        json["charge"] = self.charge
        json["charge"] = self.charge
        json["comments"] = self.comments
        json["listAtoms"] = [atom.get_object() for atom in self.listAtoms]
        return json
    def setunitDistance(self,unitDistance):
        self.unitDistance = str(unitDistance)
    def addAtomInfo(self,atom):
        #print "atomInfos: "+str(atom)
        assert (atom.atomSymbol != ""), "Trying to add an atom without symbol to a molecule"
        assert (atom.unitDistance != None), "Trying to add an atom without unit distance to a molecule"
        if self.nbAtomsInMolecule == 0 and self.unitDistance == None:
            self.unitDistance = atom.unitDistance
        elif self.nbAtomsInMolecule > 0:
            assert (str(atom.unitDistance).strip() == str(self.unitDistance).strip()),\
                'Trying to add an atom with different unitDistance than the molecule it is added to:\natom_unit={0:s} molecule_unit={1:s}'.format(str(atom.unitDistance).strip(), str(self.unitDistance).strip())
        assert isinstance(atom,atomInfos), 'Trying to add something which is not an atomInfos object to a molecule object'
        self.nbAtomsInMolecule += 1
        self.listAtoms.append(atom)
        return self


    def __str__(self):
        return self.moleculeAsString()

    def moleculeAsString(self):
        strMol =   "shortname: " + self.shortname +"\n"\
                   + "name:" + self.name +"\n" \
                   + "(" + str(self.nbAtomsInMolecule) + " atoms)\n"\
                   + "molecular charge: " + str(self.charge)+"\n" \
                   + "distances in: " + self.unitDistance+"\n" \
                   + "comments: " + self.comments+"\n"
        strMol = strMol + "".join(['{0}\n'.format(atom.get_content_atomCoord()) for atom in self.listAtoms])
        return strMol

    def getContent_format_XYZ(self):
        s = ""
        s += str(self.nbAtomsInMolecule) + "\n"
        s += self.name  + " (in Angstrom)\n"
        s += "\n".join(['{0}'.format(atom.get_content_atomCoord(printUnit=False, newUnitDistance="angstrom")) for atom in self.listAtoms])
        return s

def main():
    print "hello molecule!"
    # create an atom
    atomSymbol = 'O'
    atom = atomInfos(atomSymbol)
    print "Charge of "+ atomSymbol + " is: ",atom.atomCharge
    # create a molecule
    myMolecule = molecule("Patrickyne", name="Patrickyne Merlotusine", comments="Highly toxic large protein", nbAtomsInMolecule=5000)
    print myMolecule


if __name__ == "__main__":
    main()
