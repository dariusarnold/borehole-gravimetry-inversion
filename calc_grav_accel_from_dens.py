#!/usr/bin/env python3

import numpy as np
import itertools
import matplotlib.pyplot as plt


pi = np.pi
G = 6.6740E-11      # gravitational constant
re = 6371000        # radius earth


class Layer:
    def __init__(self, r_upper, r_lower, density_u, density_l):
        if r_lower > r_lower:
            raise ValueError("Lower radius larger than upper radius")
        self.radius_upper = r_upper
        self.radius_lower = r_lower
        self.volume = volume_sphere(r_upper) - volume_sphere(r_lower)
        self.density_upper = density_u
        self.density_lower = density_l

    @property
    def mass(self):
        # Mass calculated from average density
        return self.volume*0.5*(self.density_upper+self.density_lower)

    def __repr__(self):
        return "Layer({r_up}, {r_low}, {d_up}, {d_low}".format(r_up=self.radius_upper,
                                                               r_low=self.radius_lower,
                                                               d_up=self.density_upper,
                                                               d_low=self.density_lower)


def load_density_prem_model(filepath):
    #radius,depth,density,Vpv,Vph,Vsv,Vsh,eta,Q-mu,Q-kappa
    radius, depth, density, *trash  = np.loadtxt("PREM_1s.csv", skiprows=1, delimiter=",", unpack=True)
    # convert density from g/cm³ to kg/m³
    density = [d*1000. for d in density]
    # convert radius from km to m
    radius = [r*1000. for r in radius]
    # invert lists to get the deepest values first and the go outwards
    density.reverse()
    radius.reverse()
    # plot density PREM model read from file
    radius_inverted = [-(r-6371000.) for r in radius]
    plt.plot(density, radius_inverted)
    plt.gca().invert_yaxis()
    plt.xlabel("Density (kg/m³)")
    plt.ylabel("Radius (m)")
    plt.savefig("PREM_dens.pdf")
    plt.cla()
    layers = []
    for i in range(len(radius)-1):
        r_up = radius[i]
        r_low = radius[i+1]
        d_up = density[i]
        d_low = density[i+1]
        layers.append(Layer(r_up, r_low, d_up, d_low))
    return layers, depth


def volume_sphere(radius):
    return 4/3*pi*radius**3


def calc_grav_accel(layers):
    mass_layers = [l.mass for l in layers]
    mass_layers.insert(0, 0.)
    accel = []
    for total_mass, layer in zip(itertools.accumulate(mass_layers), layers):
        a = -G/layer.radius_upper**2*total_mass
        if np.isnan(a):
            a = 0.
        accel.append(a)
    # convert accel from m/s² to mGal, meaning 10^-5 m/s²
    accel = [a*1E5 for a in accel]
    return accel



def main(filepath):
    layers, depth = load_density_prem_model(filepath)
    accel = calc_grav_accel(layers)
    # to convert from core upwards depth axis to surface downwards depth axis
    radii = [-(l.radius_upper - 6371000.) for l in layers]
    # in a borehole depth is measured with 0 at the top and increasing towards the bottom
    # our values are starting from the core at 0 km going towards the crust at 6371 km
    accel.reverse()
    radii.reverse()
    # correct data by surface value to get relative values
    accel = [a-9.81*1E5 for a in accel]
    plt.plot(radii, accel, '.')
    plt.show()
    np.savetxt("prem_accel.dat", np.c_[np.array(radii), np.array(accel)], delimiter=' ', fmt="%.6f %.6f")


if __name__ == '__main__':
    main("PREM_1s.csv")