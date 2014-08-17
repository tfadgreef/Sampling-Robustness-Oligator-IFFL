clear,close all
% Script to perform glocal analysis check on the sampling results of the
% Oligator.
% After Hafner et al. (2009) PloS Comput. Biol. 5, e1000534
% Original scripts can be found on:
% http://www.ieu.uzh.ch/wagner/publications-software.html
randn('state',4)
% Load data
curdir=pwd;
cd ..
parg=[];
for i=1:20
    parg=[parg;load(['parS',num2str(i),'.dat'])];
end
cd(curdir)
% Number of samples for Gaussian noise.
nbsets=1e6;
set(0, 'DefaultAxesFontSize', 24);
set(0, 'defaulttextfontsize', 24);
set(0, 'defaulttextinterpreter','tex');
%%
mina=-1;maxa=1;
minT1=0;maxT1=2;
minT2=0;maxT2=2;
minT3=0;maxT3=2;
minKa=5;
maxKa=12;
% Calculate parameter values.
a=10.^(mina+(maxa-mina)*parg(:,1));
T1=10.^(minT1+(maxT1-minT1)*parg(:,2));
T2=10.^(minT2+(maxT2-minT2)*parg(:,3));
T3=10.^(minT3+(maxT3-minT3)*parg(:,4));
X0g=[a,T1,T2,T3];
Ka_aT1=10.^(minKa+(maxKa-minKa)*parg(:,5))/10^9;
Ka_aT2=Ka_aT1;
Ka_bT2=10.^(minKa+(maxKa-minKa)*parg(:,6))/10^9;
Ka_bT3=Ka_bT2;
Ka_inhT1=10.^(minKa+(maxKa-minKa)*parg(:,8))/10^9;
Ka_inhT3=(1+(5-1).*parg(:,7)).*Ka_inhT1;
Kag=[Ka_aT1,Ka_bT2,Ka_inhT1,Ka_inhT3];

% These first sets are all 'good'.
K=log10(Kag);
X0g=log10(X0g);
newgoods=[X0g,K];
lk=size(newgoods,2);
% Initialize lambda bounds.
lambdamax=2.5;
lambdamin=1.5;
maxmajor=15;
allsets{1}=[];
Rc=1.9858775e-3; % [kcal/mol]
T=42+273.15;
%%
% Main loop.
for majori = 1:maxmajor
    display(num2str(majori))
    % Lambda drops linearly from 2.5 to 1.5 and is a scaling factor to
    % slightly increase the principle components in the sampling. These
    % values are chosen on the basis of literature (Hafner et al. (2009)
    % PLoS Comput. Biol. 5, e1000534).
    lambda = lambdamin + (lambdamax-lambdamin)*max([1-(majori-1)/5 0]);
    % Add Gaussian noise in the PCA directions.
    % Determine variances.
    vark = var(newgoods);
    vark = sort((vark~=0).*[1:length(vark)],'descend');
    vark = sort(vark(1:max((vark~=0).*[1:length(vark)])));
    % Use sets for which variance is not 0.
    logks = newgoods(:,vark);
    % Perform principal component analysis for those samples whose variance
    % is not 0.
    [axis2, coord, latent] = princomp(logks);
    axis = zeros(lk,lk);
    % Directions of the principal variances.
    axis(vark,vark) = axis2;
    % First values of the variances in a cell.
    variance{majori} = zeros(lk,1);
    % Standard deviations are used, but termed variance!
    variance{majori}(vark) = sqrt(latent);
    % Use those variance in the main loop ahead.
    vari = variance{majori};
    % Determine the mean of the log10 values of concentrations and equilibria.
    lastlogmean = -Inf*ones(1,lk);
    lastlogmean(vark) = mean(logks);
    allgoods{majori} = newgoods;
    clear newgoods
    % Terminate if maximum ratio of standard deviations deviate less than
    % 10%.
    if majori>1
        if max(abs((variance{majori}./variance{majori-1}) -1))< .1
            break
        end
    end
    % Deterine new sets.
    positions = (randn(nbsets,lk)) .* (lambda*ones(nbsets,1)*vari');
    sets = (axis*positions')';
    sets = sets + ones(nbsets,1)*lastlogmean(1:lk);
    sets = 10.^sets;
    allsets{majori+1} = sets;
    tic
    for i = 1:size(sets,1)
        % Function to classify sets. Explicit bounds are used for deltaG,
        % i.e. higher than -25 kcal/mol but lower than -8 kcal/mol.
        % Additionally, inhibitor binds stronger to T3 than to T1
        % (dG(4)<dG(3)).
        deltaG=-Rc*T*log(sets(i,5:8)*10^9)-Rc*T*log(55);
        if min(sets(i,1:4))>=0.1&&max(sets(i,1:4))<=100&&min(deltaG([1,2,3,4]))>=-25&&max(deltaG)<=-8&&deltaG(4)<deltaG(3)
            try
                [test] = checksetsfunction(sets(i,:));
            catch
                test=0;
            end
        else
            test=0;
        end
        if test
            display([num2str(majori),'; ' num2str(i)])
            if ~exist('newgoods')
                newgoods = log10(sets(i,:));
            else
                newgoods(end+1,:) = log10(sets(i,:));
            end
            if exist('newgoods') && size(newgoods,1)==1000
                break
            end
        end
    end
    toc
    if ~exist('newgoods')
        display 'No viable set found! Try a lower maximum value for lambda. Terminating...'
        break
    end
end
save allgoods.mat allgoods
save allsets.mat allsets
save allgoodcons.mat allgoodcond
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