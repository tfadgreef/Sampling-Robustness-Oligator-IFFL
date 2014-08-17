% File to plot responses from sampling.
% 
% 
%
clear all
close all
set(0, 'DefaultAxesFontSize', 24);
set(0, 'defaulttextfontsize', 24);
set(0, 'defaulttextinterpreter','tex');
mina=-1;
maxa=1;
minT1=0;
maxT1=2;
minT2=0;
maxT2=2;
minT3=0;
maxT3=2;
minKa=5;
maxKa=12;
%% Parameters
% Load relevant parameter samples.
parS=[];
for i=1:20
    parS=[parS;load(['parS',num2str(i),'.dat'])];
end
% Indicate which sets should be plotted (here, the first 100).
for n=1:100
    % Calculate parameter values.
    samples=parS;
    a=10^(mina+(maxa-mina)*samples(n,1));
    T1=10^(minT1+(maxT1-minT1)*samples(n,2));
    T2=10^(minT2+(maxT2-minT2)*samples(n,3));
    T3=10^(minT3+(maxT3-minT3)*samples(n,4));
    Ka_aT1=10^(minKa+(maxKa-minKa)*samples(n,5))/10^9;
    Ka_aT2=Ka_aT1;
    Ka_bT2=10^(minKa+(maxKa-minKa)*samples(n,6))/10^9;
    Ka_bT3=Ka_bT2;
    Ka_inhT1=10^(minKa+(maxKa-minKa)*samples(n,8))/10^9;
    Ka_inhT3=(1+(5-1)*samples(n,7))*Ka_inhT1;
    % Set parameter, initial condition and time vector.
    par = [1200; 0.06; Ka_bT2; 80; 80; 3.5; 44; 440; 440; 44; ...
       44; 44; 150; 44; Ka_inhT1; Ka_inhT3; 300; 300; 30; 30; ...
       30; Ka_aT1; Ka_aT2; 1200; 1200; 0.06; 1200; 0.06; 0.06; 40; ...
       1200; 300; 80; Ka_bT3; 0.06/Ka_aT1; 0.06/Ka_aT2; 0.06/Ka_bT2; ...
        0.06/Ka_bT3; 0.06/Ka_inhT3; 0.06/Ka_inhT1; 0.06*Ka_aT1/Ka_inhT1];
    
    x0  = [T1; a; T2; T3; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; ...
        0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
    tspan=[0:7000];
    %% Solver
    % Call ode15s (delete/comment if you use mex-file).
%     options = odeset('AbsTol', 1e-11, 'RelTol', 1e-11);
%     [t,x]=ode15s(@(t,x)OligatorFunction(t,x,par,0),tspan,x0,options);
    % Call MEX file (delete/comment if you use ode15s,m and also delete
    % transpose.
    [t,x]=Oligator(tspan, x0, par, 0, [1e-11, 1e-11, 10]);
    t=t';
    x=x';
    %% Figures
    % x(1) = [T1]
    % x(2) = [a]
    % x(3) = [T2]
    % x(4) = [T3]
    % x(5) = [a.T1]
    % x(6) = [T1.a]
    % x(7) = [a.T2]
    % x(8) = [a.T1.a]
    % x(9) = [aa.T1]
    % x(10) = [ab.T2]
    % x(11) = [a.T2.b]
    % x(12) = [b]
    % x(13) = [T2.b]
    % x(14) = [b.T3]
    % x(15) = [bi.T3]
    % x(16) = [b.T3.i]
    % x(17) = [i]
    % x(18) = [T3.i]
    % x(19) = [i.T1]

    figure
    plot(t,x(:,5)+x(:,6)+x(:,7)+x(:,8)+x(:,9)+x(:,10)+x(:,11)+x(:,13)+x(:,14)+x(:,15)+x(:,16)+x(:,18)+x(:,19),'r','LineWidth',3)
    axis([0 2000 0 round(max(x(:,5)+x(:,6)+x(:,7)+x(:,8)+x(:,9)+x(:,10)+x(:,11)+x(:,13)+x(:,14)+x(:,15)+x(:,16)+x(:,18)+x(:,19)))])
    set(gca,'LineWidth',2)
    hold on
    plot(t,x(:,17),'g','LineWidth',3)
    plot(t,x(:,2),'k','LineWidth',3)
    plot(t,x(:,12),'b','LineWidth',3)
    
    xlabel('Time /min')
    ylabel('Concentration /nM')
    legend('[dsDNA]','[inh]','[\alpha]','[\beta]')

end
% Rik van Roekel, 2014
% After Montagne et al. (2011), Mol.Syst.Biol. 7, 466
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

