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

nRho=1401;
nT=1601;

minRho=1e-4;
maxRho=1000;
minT=1;
maxT=1e8;

% Add one grid point at 1e-25 to have larger density range covered
rhoAxis=[1e-25 logspace(log10(minRho),log10(maxRho),nRho)];
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