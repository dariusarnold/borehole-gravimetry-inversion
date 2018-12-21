import filecmp
import os
import unittest
import warnings

import numpy as np

import pyGravInv


class TestInversionOutput(unittest.TestCase):
    """
    This class calls the inversion function for the three norms and compares the fileoutput
    to a "known good" one. This allows easy comparison whether something broke the implementation.
    """

    def setUp(self):
        self.norms = {"L2": 0,
                      "W12": 1,
                      "Seminorm": 2}
        self.input_file = os.path.join("data", "input_inversion.dat")
        self.partial_correct_result_file = os.path.join("data", "output_inversion_{norm}.dens")

    @staticmethod
    def call_inversion_function(filename, norm_id):
        """
        Call the inversion program and return the filename of the output file
        :param filename:
        :param norm_id:
        :return:
        """
        (depth, dens), params = pyGravInv.inversion(filename, pyGravInv.Norm(norm_id), discretization_steps=10000)
        return depth, dens

    def test_norms(self):
        for norm_name, norm_id in self.norms.items():
            with self.subTest(norm_name=norm_name, norm_id=norm_id):
                correct_file = self.partial_correct_result_file.format(norm=norm_name)
                correct_depth, correct_dens = np.loadtxt(correct_file, unpack=True)
                actual_depth, actual_dens = self.call_inversion_function(self.input_file, norm_id)
                np.testing.assert_array_almost_equal(correct_depth, actual_depth, decimal=2)
                np.testing.assert_array_almost_equal(correct_dens, actual_dens, decimal=2)


class TestInversionOutputErrors(unittest.TestCase):
    """
    Test whether the inversion with errors gives the same output as expected
    """

    def setUp(self):
        self.norms = {"L2_errors": 0,
                      "Seminorm_errors": 1
                      }
        self.input_file = os.path.join("data", "input_inversion_errors.dat")
        self.partial_correct_result_file = os.path.join("data", "output_inversion_{norm}.dens")

    @staticmethod
    def call_inversion_function(filename, norm_id):
        (depth, dens), params = pyGravInv.inversion_error(filename, pyGravInv.ErrorNorm(norm_id), discretization_steps=10000)
        return depth, dens

    def test_norms(self):
        for norm_name, norm_id in self.norms.items():
            with self.subTest(norm_name=norm_name, norm_id=norm_id):
                correct_file = self.partial_correct_result_file.format(norm=norm_name)
                correct_depth, correct_dens = np.loadtxt(correct_file, unpack=True)
                actual_depth, actual_dens = self.call_inversion_function(self.input_file, norm_id)
                np.testing.assert_array_almost_equal(correct_depth, actual_depth, decimal=2)
                np.testing.assert_array_almost_equal(correct_dens, actual_dens, decimal=2)


class TestInterpolationOutput(unittest.TestCase):
    """
    This calls the interpolation function and compares the file output to a "known good" one.
    """
    def setUp(self):
        self.input_file = os.path.join("data", "input_interpolation.dat")
        self.correct_result_file = os.path.join("data", "output_interpolation.int")

    @staticmethod
    def call_interpolation_function(filename):
        """
        Call the interpolation program and return the x, f(x) values of the interpolated function
        """
        (x, y), params = pyGravInv.interpolate(filename, pyGravInv.InterpolationNorm.LinearInterpolationNorm, 0, 10, 10000)
        return x, y

    def test_interpolation(self):
        correct_file = os.path.join("data", "output_interpolation.int")
        actual_x, actual_y = self.call_interpolation_function(self.input_file)
        correct_x, correct_y = np.loadtxt(correct_file, unpack=True)
        np.testing.assert_array_almost_equal(correct_x, actual_x, decimal=2)
        np.testing.assert_array_almost_equal(correct_y, actual_y, decimal=2)

def main():
    # Calling the script from the outside wont set cwd automatically, so relative paths only work when
    # this script is called from the folder where it is saved, not from any other path. Change cwd manually to
    # make this script callable from anywhere.
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
    unittest.main()


if __name__ == '__main__':
    main()
