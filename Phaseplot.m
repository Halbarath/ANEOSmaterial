clear
close all
clc

A=textread('axes.in');

rho = A(3:A(1)+2);
T = A(A(1)+3:end);

B=textread('phase.txt');

[X,Y]=meshgrid(rho,T);

surface(X,Y,B)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
shading flat
colorbar

xlabel('Density [g/cm3]')
ylabel('Temperature [K]')
title('Phase indicator for DUNITE')
print2pdf('phaseDUNITE')