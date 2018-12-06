"""
Genereate snythetic density measurements
"""

import matplotlib.pyplot as plt
from DensityModel import make_model_arange


def main():
    depths, dens = make_model_arange(1000, 10, 5, 2000, 30, 1)
    plt.plot(depths, dens, ".")
    plt.show()


if __name__ == '__main__':
    main()
