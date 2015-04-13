#!/usr/bin/env python

import unittest
import xyz
import numpy as np

# http://docs.python-guide.org/en/latest/writing/tests/
# http://stackoverflow.com/questions/1896918/running-unittest-with-typical-test-directory-structure

class read_xyz_format(unittest.TestCase):
    def setUp(self):
        self.path_to_file1 = "./files/histidine.data"
    def test_extracting_data_from_XYZ_format(self):
        str1 = str(xyz.parse_XYZ(self.path_to_file1))
        stringDATA = ""
        with open(self.path_to_file1, 'r') as f:
            stringDATA = f.read()
            print "-----------------------------------------"
            print "toto:\n"
            print stringDATA
            print "-----------------------------------------"
        print "-----------------------------------------"
        print "tata:\n"
        print str1
        print "-----------------------------------------"
        self.assertEqual(str1.strip(), stringDATA.strip())
if __name__ == '__main__':
    unittest.main()
