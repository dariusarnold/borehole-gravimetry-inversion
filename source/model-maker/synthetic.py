"""
Genereate snythetic density measurements
"""

import numpy as np

from DensityModel import DensityModel
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


def main():
    depths = np.array((25., 66., 143.))
    save_synthetic_data("test.dat", 2700, 3500, 40, 5, depths)


if __name__ == '__main__':
    main()
