close all
% Load data

set(0, 'DefaultAxesFontSize', 24);
set(0, 'defaulttextfontsize', 24);
set(0, 'defaulttextinterpreter','tex');
parg=[];
for i = 1:20
    eval(sprintf('parg=[parg;load(''parS%s.dat'')];', num2str(i)));
end
mina=-1;maxa=1;
minT1=0;maxT1=2;
minT2=0;maxT2=2;
minT3=0;maxT3=2;
minKa=5;
maxKa=12;
% Load relevant sets.
% Load relevant sets.
a=10.^(mina+(maxa-mina)*parg(:,1));
T1=10.^(minT1+(maxT1-minT1)*parg(:,2));
T2=10.^(minT2+(maxT2-minT2)*parg(:,3));
T3=10.^(minT3+(maxT3-minT3)*parg(:,4));
X0g=[a,T1,T2,T3];
Kd_aT1=1./(10.^(minKa+(maxKa-minKa)*parg(:,5))/10^9);
Kd_aT2=Ka_aT1;
Kd_bT2=1./(10.^(minKa+(maxKa-minKa)*parg(:,6))/10^9);
Kd_bT3=Ka_bT2;
Kd_inhT1=10.^(minKa+(maxKa-minKa)*parg(:,8))/10^9;
Kd_inhT3=1./((1+(5-1).*parg(:,7)).*Ka_inhT1);
Kd_inhT1=1./Kd_inhT1;

% Group dissociation constants
K=[Kd_aT1,Kd_bT2,Kd_inhT3,Kd_inhT1];
S=ones(size(K,1),1);
% And plot scatterplots.
K2=plotCorrelations(K,S,-3,4,'g',1);
% Write figure
% export_fig -r600 -painters -nocrop -transparent testKd.pdf

% Group concentrations
K=[a,T1,T2,T3];
S=ones(size(K,1),1);
S(7)=S(7);
% And plot scatterplots.
K3=plotCorrelations(K,S,-1,2,'g',0);
% Write figure
% export_fig -r600 -painters -nocrop -transparent test.pdf

%
% Rik van Roekel, 2014
%
% Licensing:
%   Copyright (C) 2011-2015 Rik van Roekel. All rights
%   reserved.
%
%   Contact: rikvanroekel@gmail.com
%
%   This file is part of the Automated DNA Design Software (ADDS).
%   
%   ADDS is free software: you can redistribute it 
%   and/or modify it under the terms of the GNU General 
%   Public License as published by the Free Software 
%   Foundation, either version 3 of the License, or (at 
%   your option) any later version.
%
%   ADDS is distributed in the hope that it will be
%   useful, but WITHOUT ANY WARRANTY; without even the 
%   implied warranty of MERCHANTABILITY or FITNESS FOR A 
%   PARTICULAR PURPOSE.  See the GNU General Public 
%   License for more details (http://www.gnu.org/licenses/).