#!/usr/bin/env python3
#   This file is part of ANEOSmaterial.
#   Copyright (c) 2020-2021 Thomas Meier & Christian Reinhardt
#
#   ANEOSmaterial is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   ANEOSmaterial is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ANEOSmaterial.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np

nRho = 1401
nT = 1601

minRho = 1e-4
maxRho = 1000
minT = 1
maxT = 1e8

# Add one grid point at 1e-25 to have larger density range covered
rhoAxis = np.concatenate(([1e-25], np.logspace(np.log10(minRho), np.log10(maxRho), nRho)))
TAxis = np.logspace(np.log10(minT), np.log10(maxT), nT)

nT = len(TAxis)
nRho = len(rhoAxis)
with open('axes.in', 'w') as fid:
    fid.write(f'{nRho}\n')
    fid.write(f'{nT}\n')
    for rho in rhoAxis:
        fid.write(f'{rho:.15e}\n')
    for T in TAxis:
        fid.write(f'{T:.15e}\n')