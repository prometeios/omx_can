"""
@brief: Tests for omx_can package.
@usage: python -m unittest discover
"""
import os, sys

import unittest

def load_tests(loader, standard_tests, pattern):
    start_dir = './tests'
    pattern = 'test_*.py'
    print(f"Discovering tests in start_dir='{start_dir}' with pattern='{pattern}'")
    
    # Get the package's test modules
    if not hasattr(loader, '_discovered'):
        loader._discovered = True
        package_tests = loader.discover(start_dir=start_dir, pattern=pattern)
        standard_tests.addTests(package_tests)
    
    return standard_tests

def setup_package():
    """
    Set up any required resources or state for the entire test package.
    
    This function is called once before running any tests in the package.
    """
    print("Setting up the test package")
    global SAMPLE_DIR
    SAMPLE_DIR = './test/sample'

def teardown_package():
    """
    Clean up any resources or state used by the test package.
    
    This function is called once after all tests in the package have run.
    """
    files_to_remove = []
    
    for file_name in files_to_remove:
        file_path = os.path.join(SAMPLE_DIR, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)
            
class CustomTestRunner(unittest.TextTestRunner):
    def run(self, test):
        setup_package()
        result = super().run(test)
        teardown_package()
        return result

if __name__ == '__main__':
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite = load_tests(loader, suite, pattern='test_*.py')
    
    runner = CustomTestRunner()
    runner.run(suite)
    