#!/usr/bin/env python

import unittest
import topology
import xyz
import filecmp

class check_generated_topology_files(unittest.TestCase):
    def setUp(self):
        self.path_to_XYZ_format = "../files/histidine.xyz"
        self.histidine_moleculeObject = xyz.parse_XYZ(self.path_to_XYZ_format)
        self.histidine_topologyObject = topology.topology(self.histidine_moleculeObject)
        [self.filename_config, self.filename_bonds] = self.histidine_topologyObject.get_topology_files()
        self.path_to_stored_files = "./files/"

    def test_compare_topology_configFile(self):
        #with open(self.path_to_stored_files + self.filename_config, 'r') as f_stored:
        #    with open('./'+self.filename_config, 'r') as f_generated:
        #print "self.path_to_stored_files: ", self.path_to_stored_files
        #print "self.filename_config: ", self.filename_config
        res = filecmp.cmp(self.path_to_stored_files+self.filename_config, "./"+self.filename_config)
        self.assertTrue(res)

    # def test_compare_topology_covBondDistFile(self):
    #     with open(self.path_to_stored_files + self.filename_bonds, 'r') as f_stored:
    #         with open('./' + self.filename_bonds, 'r') as f_generated:
    #             filecmp.cmp(f_stored, f_generated)

class compared_topology_of_Histidine_with_known_results_from_Avogadro(unittest.TestCase):
    def setUp(self):
        self.path_to_XYZ_format = "../files/histidine.xyz"
        self.maxDiff = None
        self.histidine_moleculeObject = xyz.parse_XYZ(self.path_to_XYZ_format)
        self.histidine_topologyObject = topology.topology(self.histidine_moleculeObject)

    def test_some_bondDistances_against_AVOGADRO_results(self):
        # Histidine topology computed from XYZ input
        with open(self.path_to_XYZ_format, 'r') as f:
            self.histidine_topologyObject.get_covalentBonds()
        # check distance H2-O2 (AVOGADRO) against H1-O19 (this module)
        bondLength_Avogadro =  1.00643
        bondLength_xyz2top  = [elem.distance for elem in self.histidine_topologyObject.covalentBonds if elem.atomEntity_i.atomIndex == 1 and elem.atomEntity_j.atomIndex == 19][0]
        self.assertAlmostEqual(bondLength_Avogadro, bondLength_xyz2top, places=4)
        # check distance C-O (AVOGADRO) against C9-O19 (this module)
        bondLength_Avogadro = 1.36028
        bondLength_xyz2top  = [elem.distance for elem in self.histidine_topologyObject.covalentBonds if elem.atomEntity_i.atomIndex == 9 and elem.atomEntity_j.atomIndex == 19][0]
        self.assertAlmostEqual(bondLength_Avogadro, bondLength_xyz2top, places=4)
        # check distance C-O (AVOGADRO) against C9-O18 (this module)
        bondLength_Avogadro = 1.23144
        bondLength_xyz2top  = [elem.distance for elem in self.histidine_topologyObject.covalentBonds if elem.atomEntity_i.atomIndex == 9 and elem.atomEntity_j.atomIndex == 18][0]
        self.assertAlmostEqual(bondLength_Avogadro, bondLength_xyz2top, places=4)

    def test_some_angles_against_AVOGADRO_results(self):
        # Histidine topology computed from XYZ input
        with open(self.path_to_XYZ_format, 'r') as f:
            self.histidine_topologyObject.get_covalentBonds()
            self.histidine_topologyObject.get_covalentBondAngles()
        # check angle H2-O2-C1 (AVOGADRO) against H1-O19-C9 (this module)
        angle_Avogadro = 105.4387
        angle_xyz2top  = [elem.get_angle() for elem in self.histidine_topologyObject.covalentBondAngles if elem.atomEntity_j.atomIndex == 19][0]
        self.assertAlmostEqual(angle_Avogadro, angle_xyz2top, places=4)
        # check angle O1-C1-O2 (AVOGADRO) against O18-C9-O19 (this module)
        angle_Avogadro = 123.9856
        angle_xyz2top  = [elem.get_angle() for elem in self.histidine_topologyObject.covalentBondAngles if elem.atomEntity_j.atomIndex == 9 and elem.atomEntity_i.atomIndex == 18][0]
        self.assertAlmostEqual(angle_Avogadro, angle_xyz2top, places=4)
        # check angle O2-C1-C4 (AVOGADRO) against O19-C9-C12 (this module)
        angle_Avogadro = 113.2043
        angle_xyz2top  = [elem.get_angle() for elem in self.histidine_topologyObject.covalentBondAngles if elem.atomEntity_i.atomIndex == 12 and elem.atomEntity_j.atomIndex == 9 and elem.atomEntity_k.atomIndex == 19][0]
        self.assertAlmostEqual(angle_Avogadro, angle_xyz2top, places=4)

    def test_some_dihedral_angles_against_AVOGADRO_results(self):
        # Histidine topology computed from XYZ input
        with open(self.path_to_XYZ_format, 'r') as f:
            self.histidine_topologyObject.build_topology()
        # check Dihedral H2-O2-C1-O1 (AVOGADRO) against H1-O19-C9-O18 (this module)
        dihedral_Avogadro = -179.8814
        dihedral_xyz2top  = [elem.get_dihedral_angle() for elem in self.histidine_topologyObject.covalentDihedralAngles if elem.atomEntity_i.atomIndex == 18 and elem.atomEntity_j.atomIndex == 9 and elem.atomEntity_k.atomIndex == 19 and elem.atomEntity_l.atomIndex == 1 ][0]
        self.assertAlmostEqual(dihedral_Avogadro, dihedral_xyz2top, places=4)

        # check Dihedral H2-O2-C1-C4 (AVOGADRO) against H1-O19-C9-C12 (this module)
        dihedral_Avogadro = 2.4383
        dihedral_xyz2top  = [elem.get_dihedral_angle() for elem in self.histidine_topologyObject.covalentDihedralAngles if elem.atomEntity_i.atomIndex == 12 and elem.atomEntity_j.atomIndex == 9 and elem.atomEntity_k.atomIndex == 19 and elem.atomEntity_l.atomIndex == 1 ][0]        
        self.assertAlmostEqual(dihedral_Avogadro, dihedral_xyz2top, places=4)


if __name__ == '__main__':
    unittest.main()
