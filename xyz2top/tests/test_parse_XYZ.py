#!/usr/bin/env python

import unittest
from xyz2top import xyz


class read_xyz_format(unittest.TestCase):
    def setUp(self):
        self.path_to_XYZ_format = "./files/histidine.data"
        self.path_to_json_topology = "./files/HistidineTopology.json"
    def test_extracting_data_from_XYZ_format(self):
        str1 = str(xyz.parse_XYZ(self.path_to_XYZ_format))
        stringDATA = ""
        with open(self.path_to_XYZ_format, 'r') as f:
            stringDATA = f.read()
        self.assertEqual(str1.strip(), stringDATA.strip())

    def test_jsonOutput_Histidine_molecule(self):
        import json
        with open(self.path_to_json_topology) as json_file:
            json_topology = json.load(json_file)
        with open(self.path_to_XYZ_format, 'r') as f:
            histidine = xyz.parse_XYZ(self.path_to_XYZ_format)
            histidineXYZ_Object = histidine.get_object()
            for key, value in histidineXYZ_Object.iteritems():
                if (key not in ['comments']):
                    self.assertEqual(value,  json_topology["molecule"][key])

if __name__ == '__main__':
    unittest.main()
