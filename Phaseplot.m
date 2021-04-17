%   This file is part of ANEOSmaterial.
%   Copyright (c) 2020-2021 Thomas Meier & Christian Reinhardt
% 
%   ANEOSmaterial is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   ANEOSmaterial is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with ANEOSmaterial.  If not, see <http://www.gnu.org/licenses/>.

clear
close all
clc

A=textread('axes.in');

rho = A(3:A(1)+2);
T = A(A(1)+3:end);

% B=textread('phase_ice.txt');
B=textread('phase_dunite.txt');

[X,Y]=meshgrid(rho,T);

surface(X,Y,B)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
shading flat
colorbar

xlabel('Density [g/cm^3]')
ylabel('Temperature [K]')
title('Phase indicator for Dunite')
% print2pdf('phaseDUNITE')
hold on
% plot3(1.11,273.15,1e10,'x')
% plot3(1.11,373.15,1e10,'x')
% plot3(7.85,1811,1e10,'x')
% plot3(7.85,3273,1e10,'x')
axis tight
exportgraphics(gcf,'phase_dunite.png','Resolution',600)