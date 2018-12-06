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
        self.spke_width = spike_width
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
        :return: density at this position
        """
        if self.spike_position_top <= depth <= self.spike_position_bottom:
            return self.spike_density
        else:
            return self.background_density


class DiscretizedDensityModel:
    """
    Class that holds the result of DensityModel evaluated at discrete depths
    """
    def __init__(self, depths, densities):
        """
        :param depths: np.array of all depths points
        :param densities: np.array of the densities associated with the depth points
        """
        self.depths = depths
        self.densities = densities
        # depth_deltas: distance between to depth points in m
        self.depths_delta = np.diff(depths)


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
