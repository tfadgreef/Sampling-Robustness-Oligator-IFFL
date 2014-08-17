randn('state',7)
% Script for Monte Carlo integration of the results of the glocal
% robustness analysis on the Oligator.
% After Hafner et al. (2009) PloS Comput. Biol. 5, e1000534
% Original scripts can be found on:
% http://www.ieu.uzh.ch/wagner/publications-software.html

clear all
load allgoods.mat
nbsets=1e6;
% Take final set.
goodks=allgoods{end};
% Determine variances.
vark = var(goodks);
vark = sort((vark~=0).*[1:length(vark)],'descend');
vark = sort(vark(1:max((vark~=0).*[1:length(vark)])));
logks=goodks(:,vark);
lk=size(logks,2);
% Principal axes determination.
[axis2, coord, latent] = princomp(logks);
axis = zeros(lk,lk);
axis(vark,vark) = axis2;
variance = zeros(lk,1);
% Standard deviations, even though this variable is still called variance!
variance(vark) = sqrt(latent);
lastlogmean = -Inf*ones(1,lk);
lastlogmean(vark) = mean(logks);

% Determine hyperbox size.
maxcoord = zeros(1,lk);
maxcoord(vark) = max(coord);
mincoord = zeros(1,lk);
mincoord(vark) = min(coord);

V0 = prod(maxcoord(vark)-mincoord(vark));
reallk = length(vark);

% Generate the parameters sets.
% Sets are equvalently distribute along the major axis in the min/max
% interval.
positions = (rand(nbsets,lk)) .* ( ones(nbsets,1)*(maxcoord-mincoord)) + ( ones(nbsets,1)*mincoord);

sets = (axis*positions')';
sets = sets + ones(nbsets,1)*lastlogmean(1:lk);
sets = 10.^sets;

for i = 1:size(sets,1)
    % Function to classify sets. Explicit bounds are no longer used, since
    % this process was converged.
    try
        [test, S, P] = checksetsfunction(sets(i,:));
    catch
        test=0;
        S=0;
        P=0;
    end
    if test
        display(' ')
        display([num2str(i)])
        display(' ')
        if exist('newgoods')==0
            newgoods = log10(sets(i,:));
        else
            newgoods(end+1,:) = log10(sets(i,:));
        end
        if exist('goodconditions')==0
            goodconditions=[log10(S),log10(P)];
        else
            goodconditions(end+1,:)=[log10(S),log10(P)];
        end
    end
end
V = V0*size(newgoods,1)/size(sets,1);
% R is a measure for the overall, global robustness of the Oligator. R
% can be interpreted as the maximum order of magnitude that parameters may
% change while retaining sustained oscillations.
R =  V^(1/length(vark));
time=clock;
save V.mat V
save R.mat R
lastset=log10(sets);
save lastset.mat lastset
lastgood=newgoods;
save lastgood.mat lastgood
lastgoodconds=goodconditions;
save lastgoodconds.mat lastgoodconds
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
