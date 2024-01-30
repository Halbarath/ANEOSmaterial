#!/usr/bin/env python3
"""
Plot the pressure table generated with writePhaseDirect.
"""
from __future__ import print_function
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import argparse

def load_axes(axes_file):
    data = np.loadtxt(axes_file)

    num_rho = int(data[0])
    num_T = int(data[1])

    rho = data[2:num_rho+2]
    T = data[num_rho+2:]

    rho_min = np.min(rho)
    rho_max = np.max(rho)

    T_min = np.min(T)
    T_max = np.max(T)

    print("num_rho= {:} num_T= {:}".format(num_rho, num_T))
    print("rho_min= {:.3E} rho_min= {:.3E} T_min= {:.3E} T_max= {:.3E}".format(rho_min, rho_max, T_min, T_max))

    if np.size(rho) != num_rho or np.size(T) != num_T:
        print("Error reading file.")
        exit(1)

    return rho, T

def main():
    # Command line parameters
    parser = argparse.ArgumentParser(description="Plot the pressure of an ANEOS material")

    # Required positional argument
    parser.add_argument('file', type=str, help="Input file generated with writePressureDirect")

    # Optional arguments
    parser.add_argument('--axes', type=str, help="File containing the density and temperature axes", default="axes.in")

    args = parser.parse_args()

    file = args.file
    axes_file = args.axes

    print(axes_file)

    # Read axes file
    rho, T = load_axes(axes_file)

    press = np.loadtxt(file)

    # Limit the table to rho >= 1e-4
    index_rho_min = np.min(np.where(rho >= 1e-4))
    print("index_rho_min= {:}".format(index_rho_min))

    rho = rho[index_rho_min:]
    press = press[:,index_rho_min:]

    press_min = np.min(press)
    press_max = np.max(press)
  
    # Plot the pressure
    fig, ax = plt.subplots(1,1)

    plt.contour(np.log10(rho), np.log10(T), np.log10(press), levels=np.linspace(-20, np.log10(press_max), 20), colors='black', linewidths=0.5)
    plt.contourf(np.log10(rho), np.log10(T), np.log10(press), levels=np.linspace(-20, np.log10(press_max), 20), extend='min')
    plt.colorbar()

    plt.xlabel("Log(Density) [g cm$^{-3}$]")
    plt.ylabel("Log(Temperature) [K]")

    output = "{:}.press.png".format(file)
    plt.savefig(output, dpi=300, bbox_inches='tight')

    exit(0)


if __name__ == '__main__':
    main()

