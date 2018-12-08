import scipy.integrate as integrate
import numpy as np
from DensityModel import DensityModel, DiscretizedDensityModel


def convert_ms2_to_mGal(val):
    """
    Convert acceleration from m/s² to mGal
    """
    return val*1E5


class SyntheticGenerator:
    """Generate synthetic borehole gravity measurement data from a DensityModel"""

    FREE_AIR_GRADIENT = 0.308 # mGal/m
    GAMMA = 4*np.pi*6.674E-11 # m³/(kg*s²)

    @staticmethod
    def calculate_gravity_measurement_continuous(density_model, measurement_depths):
        """
        Calculate synthetic gravity measurements for the continuous density model at the given measurement depths.
        :param density_model: Instance of type DensityModel
        :param measurement_depths: ndarray of depths in m where measurement data should be calculated
        :return: measurement values in mGal at the measurement depths
        """
        # calculate delta_g(z_n) = f*z_n - \gamma \integral_0^n \rho(z) dz
        # this is f * z_n, unit mGal
        first_term = SyntheticGenerator.FREE_AIR_GRADIENT * measurement_depths
        # now integrate density but mind the singularity at the spike borders
        singularity_points = (density_model.spike_position_top, density_model.spike_position_bottom)
        # integral and gamma is in units of kg/m²
        integral = np.array([integrate.quad(density_model._eval, 0, upper_limit, points=singularity_points)[0] for upper_limit in measurement_depths])
        # that means we have to convert them
        result = first_term - convert_ms2_to_mGal(SyntheticGenerator.GAMMA * integral)
        return result

    @staticmethod
    def calculate_gravity_measurement_discrete(discrete_density_model, measurement_depths):
        """
        Calculate synthetic gravity measurements for the discrete density model at the given measurement depths.
        :param discrete_density_model: Instance of type DiscreteDensityModel
        :param measurement_depths: ndarray of depths in m where measurement data should be calculated
        :return: measurement values in mGal at the measurement depths
        """
        ddm = discrete_density_model
        # calculate delta_g(z_n) = f*z_n - \gamma \sum_(i=0)^n \rho(z_i) \delta z
        # this is f * z_n, unit mGal
        first_term = SyntheticGenerator.FREE_AIR_GRADIENT * measurement_depths
        # now sum up density
        integral = np.array([integrate.quad(ddm._eval, 0, upper_limit)[0] for upper_limit in measurement_depths])
        #integral and gamma is in SI units that means we have to convert it
        result = first_term - convert_ms2_to_mGal(SyntheticGenerator.GAMMA * integral)
        return result


    @staticmethod
    def save_to_file(fname, density_model, measurement_depths, measurement_errors=None):
        """
        Save the acceleration values for the given density model to a text file.
        If errors are given, save them as well.
        Data is saved in three columns: depth (m), gravity acceleration (mGal), measurement error (mGal)
        :param fname: Filename in which synthetic data is saved
        :param density_model: Density model for which snythetic measurements will be calculated
        :param measurement_depths: depths were synthetic measurements should be calculated
        :param measurement_errors: if not None, this param describes the measurement errors. If scalar value,
        assume all data has the same error, if a list of values is given it has to have the same number of
        entries as measurement_depths and it is assumed that every entry describes the error at the associated depth.
        :return:
        """
        if isinstance(density_model, DensityModel):
            measurement_data = SyntheticGenerator.calculate_gravity_measurement_continuous(density_model, measurement_depths)
        elif isinstance(density_model, DiscretizedDensityModel):
            measurement_data = SyntheticGenerator.calculate_gravity_measurement_discrete(density_model, measurement_depths)
        if measurement_errors is not None:
            try:
                num_elements = len(measurement_errors)
                if num_elements != len(measurement_depths):
                    raise IndexError("measurement_errors must have the same length as measurement_depths")
            except TypeError:
                # construct tuple of errors from single constant value
                measurement_errors = (measurement_errors,) * len(measurement_depths)
            data = (measurement_depths, measurement_data, measurement_errors)
        else:
            data = (measurement_depths, measurement_data)
        np.savetxt(fname, np.transpose(data), delimiter=' ', fmt="%.2f")