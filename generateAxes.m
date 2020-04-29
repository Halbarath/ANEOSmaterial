clear
close all
clc

nRho=1401;
nT=1201;

minRho=1e-4;
maxRho=1000;
minT=1;
maxT=1e6;

rhoAxis=logspace(log10(minRho),log10(maxRho),nRho);
% rhoAxis1 = logspace(log10(minRho),log10(3),201);
% rhoAxis2 = logspace(log10(3),log10(7),2001);
% rhoAxis3 = logspace(log10(7),log10(100),51);
% rhoAxis= [rhoAxis1, rhoAxis2(2:end), rhoAxis3(2:end)];
TAxis=logspace(log10(minT),log10(maxT),nT);

% TAxis1 = logspace(log10(minT),log10(1000),351);
% TAxis2 = logspace(log10(1000),log10(10000),601);
% TAxis3 = logspace(log10(10000),log10(1e6),101);
% TAxis = [TAxis1, TAxis2(2:end), TAxis3(2:end)];

nT=length(TAxis);
nRho=length(rhoAxis);
fid = fopen('axes.in','w');
fprintf(fid,'%d\n',nRho);
fprintf(fid,'%d\n',nT);
for i=1:nRho
    fprintf(fid,'%.15e\n',rhoAxis(i));
end
for i=1:nT
    fprintf(fid,'%.15e\n',TAxis(i));
end
fclose(fid);