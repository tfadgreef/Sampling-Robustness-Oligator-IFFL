set(0, 'DefaultAxesFontSize', 24);
set(0, 'defaulttextfontsize', 24);
set(0, 'defaulttextinterpreter','tex');
% Script to perform local perturbations on all parameters in the Oligator.
% After Hafner et al. (2009) PloS Comput. Biol. 5, e1000534
% Original scripts can be found on:
% http://www.ieu.uzh.ch/wagner/publications-software.html
randn('state',2)
load lastgood.mat

Rc=1.9858775e-3; % [kcal/mol]
T=42+273.15;
% Check what X is by removing from lastgood the samples that do not exhibit
% deltaG(4) < deltaG(3). The size that results, is X.
Rho=zeros(X,20);
for i =1:20
    eval(['load rho',num2str(i)])
    Rho(:,i)=rho;
end
A=[];

% Choose a parameter set that exhibits high robustness. (Here the first,
% the loop should be fully completed to the set you need).
for i = 1:1
    
    % Extract position of the maximum:
    column=find(max(Rho)==max(max(Rho)));
    indmax=find(Rho(:,column(1))==max(Rho(:,column(1))));
    r=Rho(indmax(1),column(1));
    A=[A;indmax(1),column(1)];
    Rho(indmax(1),column(1))=0;
    conc=lastgood(:,1:5);
    Klog=lastgood(:,6:10);
    Klog2=log10((10.^Klog).*10^9);
    deltaG=-Rc*T*Klog2*log(10)-Rc*T*log(55);
    %     Extract relevant parameter sets, as deltaG > 0 make no sense physically,
    %     when dealing with polymerized DNA from a template and frankly,
    %     deltaG > -10 makes no sense practically.
    validdG=deltaG(deltaG(:,3)<deltaG(:,5),:);
    validconc=conc(deltaG(:,3)<deltaG(:,5),:);
end
conc=10.^validconc(indmax(1),:);
deltaG=validdG(indmax(1),:);
deltaG=[deltaG(1) deltaG(2) deltaG(2) deltaG(2) deltaG(3) deltaG(4) deltaG(5)];
conc_p=repmat(conc,1000,1);
Rholocal=zeros(1,7);
for i = 1:7
    % Perturb each parameteer 1000 times, and the determine rho for each
    % parameter. This measure is used as a weithing factor in DNA
    % optimization.
    dG_p=repmat(deltaG,1000,1);

    dG_p(:,i)=dG_p(:,i)+2*randn(1000,1);
    
    sets=[conc_p,exp(-(dG_p+Rc*T*log(55))/(Rc*T))/10^9];
    for j=1:size(sets,1)
        % Function to classify sets. Explicit bounds are no longer used, since
        % this process was converged.
        [test, S, P] = checksetsfunctionsplit(sets(j,:));
        if test
            Rholocal(i)=Rholocal(i)+1;
            display(' ')
            display([num2str(i),'; ',num2str(j)])
            display(' ')
        end
    end
end
save Rholocal.mat Rholocal
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