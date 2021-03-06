"""
Generate synthetic density measurements
"""

import numpy as np

from pyGravInv.ModelMaker.DensityModel import DensityModel
from pyGravInv.ModelMaker.SyntheticGenerator import SyntheticGenerator
from pyGravInv.ModelMaker.DensityModel import make_model_fromfile


def save_synthetic_data(fname, density_background, density_spike, spike_top, spike_width, measurement_depths,
                        measurement_errors=None, save_dm=False):
    """
    Generate synthetic data for a density model with constant background density and one density spike.
    Data is saved to file
    :param measurement_errors: If a ndarray of floats is given, will be saved with the synthetic measurements. Hast
    to have same length as measurement_depths
    :param save_dm: If density model is discretized and saved as well
    :param fname: filename in which data is saved
    :param density_background: background density in kg/m³
    :param density_spike: density of spike in kg/m³
    :param spike_top: depth of the top of the spike layer in m
    :param spike_width: in m
    :param measurement_depths: iterable of depths at which synthetic data should be created
    :return: None
    """
    dm = DensityModel(density_background, spike_top, spike_width, density_spike)
    measurement_depths = np.array(measurement_depths)
    if measurement_errors is not None:
        measurement_errors = np.array(measurement_errors)
    if save_dm:
        path, extension = fname.split(".")
        dfname = path + ".dens"
        dm.save_to_file(dfname, 10000, 0, measurement_depths[-1])
    SyntheticGenerator.save_to_file(fname, dm, measurement_depths, measurement_errors)


def generate_synthetic_measurement():
    depths = np.linspace(0, 150, 15)
    meas_errors = np.array((0.1,)*len(depths))
    save_synthetic_data("ex94.dat", 2700, 3500, 40, 5, depths, meas_errors, save_dm=True)


def recreate_measurement_data():
    # read density model that was created as an inversion result
    dm = make_model_fromfile("input_dens_model.dens")
    depths = np.array((25., 66., 143.))
    # evaluate this density model at the same points
    SyntheticGenerator.save_to_file("output_meas_data.dat", dm, depths)


if __name__ == '__main__':
    #recreate_measurement_data()
    generate_synthetic_measurement()
