#!/usr/bin/env python

import unittest
import topology
import xyz

class compared_topology_of_Histidine_with_known_results_from_Avogadro(unittest.TestCase):
    def setUp(self):
        self.path_to_XYZ_format = "../files/histidine.xyz"

    def test_entire_topology_against_original_one(self):
        import json
        # Histidine topology read from file
        path_to_json_topology = "./files/HistidineTopology.json"
        with open(path_to_json_topology) as json_file:
            json_topology = json.load(json_file)
        # Histidine topology computed from XYZ input
        with open(self.path_to_XYZ_format, 'r') as f:
            histidine_moleculeObject = xyz.parse_XYZ(self.path_to_XYZ_format)
            histidine_topology = topology.topology(histidine_moleculeObject)
            histidine_topology.build_topology()
            obj = histidine_topology.get_object()
            with open('./files/HistidineTopology_new.json', 'w') as outfile:
                outfile.write(str(obj))
#            for key, value in histidine_topology.iteritems():
#                #if (key not in ['comments']):
#                self.assertEqual(value,  json_topology[key])




if __name__ == '__main__':
    unittest.main()
