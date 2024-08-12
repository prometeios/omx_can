"""
Test the initial generation of the INCAR, KPOINTS, KPATH and POTCAR files.
"""
import unittest
import os

from src.omx_can import omx_io
from src.omx_can import mag_conv

SAMPLE_DIR = './tests/sample'

class TestSample(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        print("Setting up for test")

    def test_alter_init(self):
        mag_conv.alter_init_magmom(os.path.join(SAMPLE_DIR, "NbSeTe.dat"), 4)
        self.assertTrue(True, "No error")


    @classmethod
    def tearDownClass(cls):
        print("Tearing down after test")
        # Find the all the openmx output files except the restart files and delete them
        cls.after_all_tests(cls)()
    
    @classmethod
    def after_all_tests(cls):
        print("After all tests")
        # Find the all the openmx output files except the restart files and delete them
        os.sys("rm *.*")
        cls.after_all_tests(cls)()
        

if __name__ == '__main__':
    unittest.main()
    

    

    