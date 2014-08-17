function rho=local(ind1,ind2,batch)
% Function to perform glocal analysis check on the sampling results of the
% Oligator, the local part after MC integration.
% After Hafner et al. (2009) PloS Comput. Biol. 5, e1000534
% Original scripts can be found on:
% http://www.ieu.uzh.ch/wagner/publications-software.html
randn('state',13)
load lastgood.mat
load R.mat

Rc=1.9858775e-3; % [kcal/mol]
T=42+273.15;

% A dG deviation of `1 kcal/mol will be introduced for all association
% constants.
conc=lastgood(:,1:4);
Klog=lastgood(:,5:8);
Klog2=log10((10.^Klog).*10^9);
deltaG=-Rc*T*Klog2*log(10)-Rc*T*log(55);
% Extract relevant parameter sets, as deltaG > 0 make no sense physically,
% when dealing with polymerized DNA from a template and frankly,
% deltaG > -10 makes no sense practically.
validdG=deltaG(deltaG(:,3)>deltaG(:,4),:);
validconc=conc(deltaG(:,3)>deltaG(:,4),:);
rho=zeros(size(validdG,1),1);
for i = ind1:ind2
    conc_p=repmat(validconc(i,:),100,1);
    dG_p=zeros(100,size(deltaG,2));
    for j=1:100
        dG_p(j,:)=validdG(i,:)+1*randn(1,4);
        while dG_p(j,3)<=dG_p(j,4)
            dG_p(j,:)=validdG(i,:)+1*randn(1,4);
        end
    end
    sets=[10.^conc_p, (1/55*exp(-dG_p/(Rc*T)))/10^9];
    for j=1:100
        % Function to classify sets. Explicit bounds are no longer used, since
        % this process was converged.
        try
            [test] = checksetsfunction(sets(j,:));
        catch
            test=0;
        end
        if test
            rho(i)=rho(i)+1;
            display(' ')
            display([num2str(i),'; ',num2str(j)])
            display(' ')
        end
    end
    rho(i)=rho(i)/100;
    if mod(i,50)==0
        eval(['save rho',num2str(batch),'.mat rho'])
    end
end
eval(['save rho',num2str(batch),'.mat rho'])
end
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