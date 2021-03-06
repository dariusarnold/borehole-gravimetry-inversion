import numpy as np


class DensityModel:
    """
    Class that describes a density model with constant background density and a
    jump of increased density at a certain depth with a certain thickness.
    Density model starts at depth zero and has a positive depth axis downward.
    """
    def __init__(self, background_density, spike_position_top, spike_width, spike_density):
        """
        Init model with background density and one density spike.
        Spike is inclusive, i.e. every depth value in [spike_position_top, spike_position_top+spike_width]
        will have the spike_density.
        :param background_density: value of the background density in kg/m³
        :param spike_position_top: position of the top of the density spike in m
        :param spike_width: width of the density spike in m
        :param spike_density: density of the spike in kg/m³
        """
        self.background_density = background_density
        self.spike_position_top = spike_position_top
        self.spike_position_bottom = spike_position_top + spike_width
        self.spike_density = spike_density

    @classmethod
    def make_model(cls, background_density, spike_position_top, spike_width, spike_density, eval_depths):
        """
        Create model and evaluate it at eval_depths, return evaluated/discretized density
        :param background_density: value of the background density in kg/m³
        :param spike_position_top: position of the top of the density spike in m
        :param spike_width: width of the density spike in m
        :param spike_density: density of the spike in kg/m³
        :param eval_depths: array of depth values where the model is to be evaluated
        :return: model evaluated at the given depth values
        """
        dm = cls(background_density, spike_position_top, spike_width, spike_density)
        dens = dm._eval_model(eval_depths)
        return dens

    def _eval_model(self, depths):
        """
        Evaluate the density model at a certain depth
        :param depths: np.array of depth values
        :return: np.array of the densities at the given depths
        """
        return np.vectorize(self._eval)(depths)

    def _eval(self, depth):
        """
        Eval the model at a constant depth
        :param depth: depth in m
        :return: density at this position in kg/m³
        """
        if self.spike_position_top <= depth <= self.spike_position_bottom:
            return self.spike_density
        else:
            return self.background_density

    def save_to_file(self, fname, discretization_steps, depth_top, depth_bottom):
        """
        Save a discretized version of the density model to the file.
        :param fname: Filename
        :param discretization_steps: number of steps used for discretization
        :param depth_top: Top of evaluation interval
        :param depth_bottom: Bottom of evaluation interval
        """
        depths = np.linspace(depth_top, depth_bottom, discretization_steps)
        dens = self._eval_model(depths)
        data = (depths, dens)
        np.savetxt(fname, np.transpose(data), delimiter=" ", fmt="%.2f")



class DiscretizedDensityModel:
    """
    Class that holds the result of DensityModel evaluated at discrete depths
    """
    def __init__(self, depths, densities):
        """
        :param depths: np.array of all depths points in m
        :param densities: np.array of the densities in kg/m³ associated with the depth points
        """
        self.depths = depths
        self.densities = densities
        # depth_deltas: distance between to depth points in m
        self.depths_delta = np.diff(depths)

    def _eval_model(self, depths):
        """
        Evaluate the density model at a certain depth
        :param depths: np.array of depth values
        :return: np.array of the densities at the given depths
        """
        return np.vectorize(self._eval)(depths)


    def _eval(self, depth):
        """
        Eval the model at a certain depth, returns the density value from the closest point
        :param depth: depth in m
        :return: density at this position in kg/m³
        """
        # index of the closest depth value
        index = np.abs(self.depths-depth).argmin()
        return self.densities[index]


def make_model_fromfile(fname):
    """Make discrete density model from """
    depth, densities = np.loadtxt(fname, unpack=True)
    # convert density reading from file (g/cm³) to kg/m³
    densities *= 1000
    return DiscretizedDensityModel(depth, densities)


def make_model_arange(background_density, spike_position_top, spike_width, spike_density, eval_depth_from, eval_depth_upto, eval_stepsize):
    """
    Make density model from a stepsize
    :param background_density: value of the background density in kg/m³
    :param spike_position_top: position of the top of the density spike in m
    :param spike_width: width of the density spike in m
    :param spike_density: density of the spike in kg/m³
    :param eval_depth_from: Density of model is evaluated from this depth onwards to eval_depth_upto
    :param eval_depth_upto: Density of model is evaluated from eval_depth up to this value
    :param eval_stepsize: Depth interval is evaluated using this step size
    :return: DiscretizedDensityModel holding the results, depths and associated density values
    """
    depths = np.arange(eval_depth_from, eval_depth_upto, eval_stepsize)
    dens = DensityModel.make_model(background_density, spike_position_top, spike_width, spike_density, depths)
    return DiscretizedDensityModel(depths, dens)


def make_model_linspace(background_density, spike_position_top, spike_width, spike_density, eval_depth_from, eval_depth_upto, eval_num_steps):
    """
    Make density model from a number of steps
    :param background_density: value of the background density in kg/m³
    :param spike_position_top: position of the top of the density spike in m
    :param spike_width: width of the density spike in m
    :param spike_density: density of the spike in kg/m³
    :param eval_depth_from: Density of model is evaluated from this depth onwards to eval_depth_upto
    :param eval_depth_upto: Depth is evaluated from 0 to this param
    :param eval_num_steps: Depth interval is evaluated using this number of steps
    :return: DiscretizedDensityModel holding the results, depths and associated density values
    """
    depths = np.linspace(eval_depth_from, eval_depth_upto, eval_num_steps)
    dens = DensityModel.make_model(background_density, spike_position_top, spike_width, spike_density, depths)
    return DiscretizedDensityModel(depths, dens)


def make_model_non_uniform(background_density, spike_position_top, spike_width, spike_density, eval_depths):
    """
    Make density model for given depth values, which do not need to have a uniform distance
    :param background_density: value of the background density in kg/m³
    :param spike_position_top: position of the top of the density spike in m
    :param spike_width: width of the density spike in m
    :param spike_density: density of the spike in kg/m³
    :param eval_depths: np.array of depths at which the model is to be evaluated
    :return:
    """
    dens = DensityModel.make_model(background_density, spike_position_top, spike_width, spike_density, eval_depths)
    return DiscretizedDensityModel(eval_depths, dens)
