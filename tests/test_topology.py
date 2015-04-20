#!/usr/bin/env python

import unittest
import topology


class compared_topology_of_Histidine_with_known_results_from_Avogadro(unittest.TestCase):
    def setUp(self):
        self.path_to_file1 = "../files/histidine.data"
    def compare_bond_distances(self):
        print "hello"
        #self.assertEqual(str1.strip(), stringDATA.strip())


if __name__ == '__main__':
    unittest.main()
