#!/usr/bin/env python

import unittest
import xyz


class read_xyz_format(unittest.TestCase):
    def setUp(self):
        self.path_to_file1 = "./files/histidine.data"

    def test_extracting_data_from_XYZ_format(self):
        str1 = str(xyz.parse_XYZ(self.path_to_file1))
        stringDATA = ""
        with open(self.path_to_file1, 'r') as f:
            stringDATA = f.read()
        self.assertEqual(str1.strip(), stringDATA.strip())

    def test_jsonOutput_Histidine_molecule(self):
        self.jsonHistidine = {'listAtoms': [{'coordinates': [-1.14902471332943, 1.00712023098862, 0.13591620816877], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'H', 'unitDistance': 'Angstrom', 'atomCharge': 1.0}, {'coordinates': [-1.99491807877866, -0.24117276176644, 2.89986635578188], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'H', 'unitDistance': 'Angstrom', 'atomCharge': 1.0}, {'coordinates': [-1.06662207264791, 1.82960564281643, 2.34031146027655], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'H', 'unitDistance': 'Angstrom', 'atomCharge': 1.0}, {'coordinates': [0.30614058355617, 0.88602453617217, 2.48051165179861], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'H', 'unitDistance': 'Angstrom', 'atomCharge': 1.0}, {'coordinates': [-0.12748455878107, -1.1863180134193, -0.42200736830855], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'H', 'unitDistance': 'Angstrom', 'atomCharge': 1.0}, {'coordinates': [0.33590918869435, -1.43117173867306, 1.25744790059564], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'H', 'unitDistance': 'Angstrom', 'atomCharge': 1.0}, {'coordinates': [2.83606289990766, -0.37106494576807, 1.99674765933953], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'H', 'unitDistance': 'Angstrom', 'atomCharge': 1.0}, {'coordinates': [3.48463424069596, 2.05765029096117, -1.37181243893832], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'H', 'unitDistance': 'Angstrom', 'atomCharge': 1.0}, {'coordinates': [1.14863091785471, 1.01574781834498, -1.60593641137974], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'H', 'unitDistance': 'Angstrom', 'atomCharge': 1.0}, {'coordinates': [-2.2481083752515, -0.58352224399409, 1.05302316149422], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'C', 'unitDistance': 'Angstrom', 'atomCharge': 6.0}, {'coordinates': [2.63102713882836, 0.13428360747356, 1.06613047675826], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'C', 'unitDistance': 'Angstrom', 'atomCharge': 6.0}, {'coordinates': [3.02693381274811, 1.3956465802414, -0.65561792915948], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'C', 'unitDistance': 'Angstrom', 'atomCharge': 6.0}, {'coordinates': [-0.96887358386731, 0.2678783158116, 0.92794236599846], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'C', 'unitDistance': 'Angstrom', 'atomCharge': 6.0}, {'coordinates': [0.19844811513605, -0.66422149912031, 0.48533561352453], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'C', 'unitDistance': 'Angstrom', 'atomCharge': 6.0}, {'coordinates': [1.49423218608971, 0.05048815004851, 0.28293527485411], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'C', 'unitDistance': 'Angstrom', 'atomCharge': 6.0}, {'coordinates': [-0.67900595386139, 0.89724997805987, 2.22898932368553], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'N', 'unitDistance': 'Angstrom', 'atomCharge': 7.0}, {'coordinates': [3.57939992810719, 0.97162866594325, 0.47571927252474], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'N', 'unitDistance': 'Angstrom', 'atomCharge': 7.0}, {'coordinates': [1.75985825580935, 0.86910491315057, -0.81740130110534], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'N', 'unitDistance': 'Angstrom', 'atomCharge': 7.0}, {'coordinates': [-2.81628989262937, -1.07972561038044, 0.07968358117985], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'O', 'unitDistance': 'Angstrom', 'atomCharge': 8.0}, {'coordinates': [-2.64668003898623, -0.76603191500358, 2.34072514768128], 'bohr_in_angstrom': 0.52917721092, 'atomSymbol': 'O', 'unitDistance': 'Angstrom', 'atomCharge': 8.0}], 'name': 'histidine', 'nbAtomsInMolecule': 20, 'unitDistance': 'Angstrom', 'comments': 'lsdalton20140924_geomOpt-b3lyp_Vanlenthe_6-31G_df-def2_Histidine_2CPU_16OMP_2014_10_28T1007', 'charge': None, 'shortname': 'histidine'}
        with open(self.path_to_file1, 'r') as f:
            strFile = f.read()
            histidine = xyz.parse_XYZ(self.path_to_file1)
            jsonHistidine = histidine.get_json()
            for key, value in self.jsonHistidine.iteritems():
                print key
                if (key not in ['comments']):
                    self.assertEqual(value,  jsonHistidine[key])

if __name__ == '__main__':
    unittest.main()
