import unittest
import os
import filecmp
import warnings


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
        progname = "Inversion.exe" if os.name == 'nt' else "Inversion"
        prog_path = os.path.join("..", "cmake-build-debug", progname)
        call_string = "{programm_path} {input_path} {discretization_steps} {norm_id}"
        # discretzation_steps: Number of steps that are applied to discretize the resulting density distribution.
        call_string = call_string.format(programm_path=prog_path,
                                         input_path=filename,
                                         discretization_steps=10000,
                                         norm_id=norm_id)
        os.system(call_string)
        output_filename = filename.replace(".dat", ".dens")
        return output_filename

    def test_norms(self):
        for norm_name, norm_id in self.norms.items():
            with self.subTest(norm_name=norm_name, norm_id=norm_id):
                correct_file = self.partial_correct_result_file.format(norm=norm_name)
                output_filename = self.call_inversion_function(self.input_file, norm_id)
                error_message = """
                {norm_name}: Results are not the same. Compare files:
                Expected: {expected}
                Got: {got}""".format(norm_name=norm_name, expected=correct_file, got=output_filename)
                self.assertTrue(filecmp.cmp(correct_file, output_filename, shallow=False), error_message)
                try:
                    os.remove(output_filename)
                except OSError:
                    warnings.warn("Couldn't delete file {}".format(output_filename))


class TestInversionOutputErrors(unittest.TestCase):
    """
    Test whether the inversion with errors gives the same output as expected
    """

    def setUp(self):
        self.norms = {"L2_errors": 0
                      }
        self.input_file = os.path.join("data", "input_inversion_errors.dat")
        self.partial_correct_result_file = os.path.join("data", "output_inversion_{norm}.dens")

    @staticmethod
    def call_inversion_function(filename, norm_id):
        progname = "Inversion_with_Errors.exe" if os.name == 'nt' else "Inversion_with_Errors"
        prog_path = os.path.join("..", "cmake-build-debug", progname)
        call_string = "{prog_path} {input_path} {discretization_steps} {norm_id}"
        call_string = call_string.format(prog_path=prog_path,
                                         input_path=filename,
                                         discretization_steps=10000,
                                         norm_id=norm_id)
        os.system(call_string)
        output_filename = filename.replace(".dat", ".dens")
        return output_filename

    def test_norms(self):
        for norm_name, norm_id in self.norms.items():
            with self.subTest(norm_name=norm_name, norm_id=norm_id):
                correct_file = self.partial_correct_result_file.format(norm=norm_name)
                output_filename = self.call_inversion_function(self.input_file, norm_id)
                error_message = """
                    {norm_name}: Results are not the same. Compare files:
                    Expected: {expected}
                    Got: {got}""".format(norm_name=norm_name, expected=correct_file, got=output_filename)
                self.assertTrue(filecmp.cmp(correct_file, output_filename, shallow=False), error_message)
                try:
                    pass
                    #os.remove(output_filename)
                except OSError:
                    warnings.warn("Couldn't delete file {}".format(output_filename))


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
        Call the interpolation program and return the filename of the output file
        :param filename:
        :return:
        """
        progname = "Interpolation.exe" if os.name == "nt" else "Interpolation"
        prog_path = os.path.join("..", "cmake-build-debug", progname)
        call_string = "{programm_path} {input_path} {discretization_steps} {lower_bound} {upper_bound}"
        call_string = call_string.format(programm_path=prog_path,
                                         input_path=filename,
                                         discretization_steps=10000,
                                         lower_bound=0,
                                         upper_bound=10)
        os.system(call_string)
        output_filename = filename.replace(".dat", ".int")
        return output_filename

    def test_interpolation(self):
        output_filename = self.call_interpolation_function(self.input_file)
        error_message = """
        Interpolation: Results are not the same. Compare files:
        Expected: {expected}
        Got: {got}""".format(expected=self.correct_result_file, got=output_filename)
        print(output_filename, self.correct_result_file)
        self.assertTrue(filecmp.cmp(self.correct_result_file, output_filename, shallow=False), error_message)
        try:
            os.remove(output_filename)
        except OSError:
            warnings.warn("Couldn't delete file {}".format(output_filename))


def main():
    # Calling the script from the outside wont set cwd automatically, so relative paths only work when
    # this script is called from the folder where it is saved, not from any other path. Change cwd manually to
    # make this script callable from anywhere.
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
    unittest.main()


if __name__ == '__main__':
    main()
