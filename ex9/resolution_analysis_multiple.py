import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

from source.ModelMaker.synthetic import save_synthetic_data


def call_inversion(filepath, nu=None):
    progname = "Inversion_with_Errors"
    if os.name == 'nt': progname += ".exe"
    progpath = os.path.join("..", "cmake-build-debug", progname)
    call_string = "{programm_path} {input_path} {discretization_steps} {norm_id} {nu}"
    call_string = call_string.format(programm_path=progpath,
                                     input_path=filepath,
                                     discretization_steps=10000,
                                     # this is W12 Seminorm
                                     norm_id=1,
                                     nu="" if nu is None else nu)
    try:
        result = subprocess.check_output(call_string, shell=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
    result = result.split()
    print(result)
    nu = float(result[1])
    misfit = float(result[4])
    norm = float(result[6])
    return nu, misfit, norm

def replace_extension(filepath, new_extension):
    root, ext = os.path.splitext(filepath)
    return f"{root}.{new_extension}"

def main():
    # parameters of synthetic model
    dens_bg = 2700
    dens_spike = 3500
    spike_width = 5
    depths_spike_top = np.arange(0, 151, 5)
    meas_depths = np.arange(0, 150, 15)
    meas_errors = np.array([0.1]*len(meas_depths))
    # parameters of inversion
    nu = 6.58E-4
    JUST_READ=False


    inversion_results = []
    for depth_spike in depths_spike_top:
        fname = "a4/res_ana_{}.dat".format(depth_spike)
        if not JUST_READ:
            save_synthetic_data(fname, dens_bg, dens_spike, depth_spike, spike_width, meas_depths, meas_errors)
            call_inversion(fname, nu)
        fname = replace_extension(fname, "dens")
        depth, dens = np.loadtxt(fname, unpack=True)
        inversion_results.append(dens)


    inversion_results = np.array(inversion_results)
    # this has to be done so the data is vertical instead of horizontal
    inversion_results = inversion_results.transpose()


    # plot results
    x_axis_extent = [depths_spike_top[0], depths_spike_top[-1]]
    y_axis_extent = [depth[-1], depth[0]]
    plt.imshow(inversion_results, aspect="auto", cmap="Greys", origin="upper", extent=(*x_axis_extent, *y_axis_extent))
    cbar = plt.colorbar()
    cbar.set_label("Density (g/m³)")
    plt.xlabel("Depth of spike (m)")
    plt.ylabel("Depth borehole (m)")
    plt.title("Inverted models")
    plt.show()


if __name__ == '__main__':
    main()