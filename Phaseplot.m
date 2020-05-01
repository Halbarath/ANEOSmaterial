clear
close all
clc

A=textread('axes.in');

rho = A(3:A(1)+2);
T = A(A(1)+3:end);

% B=textread('phase_ice.txt');
B=textread('phase_iron.txt');

[X,Y]=meshgrid(rho,T);

surface(X,Y,B)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
shading flat
colorbar

xlabel('Density [g/cm3]')
ylabel('Temperature [K]')
title('Phase indicator for ICE')
% print2pdf('phaseDUNITE')
hold on
% plot3(1.11,273.15,1e10,'x')
% plot3(1.11,373.15,1e10,'x')
plot3(7.85,1811,1e10,'x')
plot3(7.85,3273,1e10,'x')