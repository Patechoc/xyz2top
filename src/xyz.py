#!/usr/bin/env python 

import sys, os
import parser as prs
import atomInMolecule as mol



def parse_XYZ(path_to_file):
    stringXYZ = ""
    with open(path_to_file, 'r') as f:
        stringXYZ = f.read()

    ### PyParsing grammar
    coordinate = prs.Combine(prs.Optional(prs.Literal("-")) + prs.integer+prs.Optional(prs.Literal(".") + prs.integer))
    oneAtomCoordinates = (prs.element).setResultsName("atomAbrev") + coordinate.setResultsName("xCoord") + coordinate.setResultsName("yCoord") + coordinate.setResultsName("zCoord")
    get_infos = prs.LineStart() + (prs.integer).setResultsName("nbAtomsTotal") + prs.EOL + prs.StrangeName.setResultsName("comments")

    ## build the molecule object
    shortname = os.path.splitext(os.path.basename(path_to_file))[0]
    name = shortname
    comments = path_to_file

    ## header of XYZ file
    for tokens in get_infos.searchString(stringXYZ):
        if tokens.nbAtomsTotal:
            nbAtomsInMolecule = tokens.nbAtomsTotal
        if tokens.comments:
            comments = tokens.comments
    molecule = mol.molecule(shortname, name, comments)
    molecule.setunitDistance("Angstrom")
    ## adding atoms
    for tokens in oneAtomCoordinates.searchString(stringXYZ):
        #print tokens.dump()
        if tokens.atomAbrev != "":
            atom = mol.atomInfos(str(tokens.atomAbrev))
            atom.setUnitDistance("Angstrom")
            atom.setAtomCoord(float(tokens.xCoord), float(tokens.yCoord), float(tokens.zCoord))
            molecule.addAtomInfo(atom)
    return molecule


if __name__ == "__main__":
    ## TEST
    path_to_file = "../files/histidine.xyz"
    molecule = parse_XYZ(path_to_file)
    print str(molecule)
