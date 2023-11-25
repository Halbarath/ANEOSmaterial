#!/usr/bin/env python3
"""
Plot the pressure table generated with writePressureTable.
"""
from __future__ import print_function
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys


def main():
    fig, ax = plt.subplots(1,1)

    if len(sys.argv) != 2:
        print("Usage: plot_phase_table.py <phase.txt>")
        exit(1)

    data = np.loadtxt(sys.argv[1])

    rho = data[:,0]
    T = data[:,1]
    phase = data[:,2]

    rho_min = np.min(rho)
    rho_max = np.max(rho)

    T_min = np.min(T)
    T_max = np.max(T)
    
    phase_min = np.min(phase)
    phase_max = np.max(phase)

    num_rho = np.size(np.where(rho == rho_max))
    num_T = np.size(np.where(T == T_max))

    num_phase = np.size(np.unique(phase))
    
    print("num_rho = {:} num_T = {:} num_phase = {:}".format(num_rho, num_T, num_phase))
    print("Phases: {:}".format(np.unique(phase)))

    # rho and T axis
    rho = np.unique(rho)
    T = np.unique(T)

    phase = np.split(phase, num_rho)
   
    # Plot the pressure
    plt.imshow(phase, cmap=plt.cm.get_cmap('viridis', num_phase), origin='lower', extent=[rho_min, rho_max, T_min, T_max], aspect='auto')
    plt.colorbar()
    
    plt.clim(phase_min-0.5, phase_max+0.5)

    plt.xlabel("Density [g cm$^{-3}$]")
    plt.ylabel("Temperature [K]")

    plt.savefig("phase.png", dpi=300, bbox_inches='tight')

    exit(0)


if __name__ == '__main__':
    main()

