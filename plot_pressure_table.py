#!/usr/bin/env python3
"""
Plot the pressure table generated with writePressureTable.
"""
from __future__ import print_function
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def main():
    fig, ax = plt.subplots(1,1)

    data = np.loadtxt("pressure.txt")

    rho = data[:,0]
    T = data[:,1]
    P = data[:,2]

    rho_min = np.min(rho)
    rho_max = np.max(rho)

    T_min = np.min(T)
    T_max = np.max(T)
    
    num_rho = np.size(np.where(rho == rho_max))
    num_T = np.size(np.where(T == T_max))

    # rho and T axis
    rho = np.unique(rho)
    T = np.unique(T)

    log_press = np.log(np.split(P, num_rho))
   
    # Plot the pressure
    plt.imshow(log_press, cmap=plt.cm.viridis, origin='lower', extent=[rho_min, rho_max, T_min, T_max], aspect='auto')
    plt.colorbar()
    
    plt.xlabel("Density [g cm$^{-3}$]")
    plt.ylabel("Temperature [K]")

    plt.savefig("pressure.png", dpi=300, bbox_inches='tight')


    exit(0)


if __name__ == '__main__':
    main()

