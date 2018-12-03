import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess


input_files = ("grav3.dat", "grav15.dat")

progname = "Inversion_with_Errors"
progname = os.path.join("cmake-build-debug", progname)
if os.name == 'nt': progname += ".exe"

call_string = "{progname} {{input_path}} {discretization_steps} {norm_id} {{nu}}"
call_string = call_string.format(progname=progname, discretization_steps=10000, norm_id=0)

for fname in input_files:
    nus = (0.01, 0.1, 1., 10., 100.)
    misfits = []
    norms = []
    for nu in nus:
        _call_string = call_string.format(input_path=fname, nu=nu)
        result = subprocess.check_output(_call_string, shell=True)
        result = result.split()
        misfit, norm = result[2], result[4]
        misfits.append(float(misfit))
        norms.append(float(norm))
    # plot misfit over norm, attach eta
    plt.figure()
    plt.scatter(norms, misfits)
    # add nu markers to points
    ax = plt.gca()
    for i, nu in enumerate(nus):
        ax.annotate(str(nu), xy=(norms[i], misfits[i]))
    plt.title("Misfit over norm for {}".format(fname))
    plt.xlabel(r"$||m||$")
    plt.ylabel(r"$\frac{\chi^2}{N}$")
    # make room for ylabel
    plt.tight_layout()
    plt.show()
