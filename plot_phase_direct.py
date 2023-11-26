#!/usr/bin/env python3
"""
Plot the pressure table generated with writePhaseDirect.
"""
from __future__ import print_function
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import argparse


def main():
    # Command line parameters
    parser = argparse.ArgumentParser(description="Plot the different phases of an ANEOS input file")

    # Required positional argument
    parser.add_argument('file', type=str, help="Input file generated with writePhaseDirect")

    # Optional arguments
    parser.add_argument('--axes', type=str, help="File containing the density and temperature axes", default=None)

    args = parser.parse_args()

    file = args.file
    axes_file = args.axes

    if axes_file is not None:
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

    phase = np.loadtxt(file)
 
    phase_min = np.min(phase)
    phase_max = np.max(phase)
 
    num_phase = np.size(np.unique(phase))
    
    print("num_phase = {:}".format(num_phase))
    print("Phases: {:}".format(np.unique(phase)))
 
    # Plot the phases
    fig, ax = plt.subplots(1,1)

    if axes_file is None:
        plt.imshow(phase, cmap=plt.cm.get_cmap('viridis', num_phase), origin='lower', aspect='auto')
    else:
        plt.imshow(phase, cmap=plt.cm.get_cmap('viridis', num_phase), origin='lower', extent=[np.log10(rho_min), np.log10(rho_max), np.log10(T_min), np.log10(T_max)], aspect='auto')
    
    plt.colorbar()
    
    plt.clim(phase_min-0.5, phase_max+0.5)

    plt.xlabel("Density [g cm$^{-3}$]")
    plt.ylabel("Temperature [K]")

    output = "{:}.phase.png".format(file)
    plt.savefig(output, dpi=300, bbox_inches='tight')

    exit(0)


if __name__ == '__main__':
    main()

