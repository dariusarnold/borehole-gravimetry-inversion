"""
Genereate snythetic density measurements
"""

import numpy as np

from DensityModel import DensityModel, make_model_fromfile
from SyntheticGenerator import SyntheticGenerator


def save_synthetic_data(fname, density_background, density_spike, spike_top, spike_width, measurement_depths):
    """
    Generate synthetic data for a density model with constant background density and one density spike.
    Data is saved to file
    :param fname: filename
    :param density_background: background density in kg/m³
    :param density_spike: density of spike in kg/m³
    :param spike_top: depth of the top of the spike layer in m
    :param spike_width: in m
    :param measurement_depths: iterable of depths at which synthetic data should be created
    :return: None
    """
    dm = DensityModel(density_background, spike_top, spike_width, density_spike)
    SyntheticGenerator.save_to_file(fname, dm, measurement_depths)


def generate_synthetic_measurement():
    depths = np.array((25., 66., 143.))
    save_synthetic_data("test.dat", 2700, 3500, 40, 5, depths)


def recreate_measurement_data():
    # read density model that was created as an inversion result
    dm = make_model_fromfile("input_dens_model.dens")
    depths = np.array((25., 66., 143.))
    # evaluate this density model at the same points
    SyntheticGenerator.save_to_file("output_meas_data.dat", dm, depths)


if __name__ == '__main__':
    recreate_measurement_data()
