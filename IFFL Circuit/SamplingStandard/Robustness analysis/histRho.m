set(0, 'DefaultAxesFontSize', 24);
set(0, 'defaulttextfontsize', 24);
set(0, 'defaulttextinterpreter','tex');
Rho=[];
% Load all local robustness (rho) values.
for i=1:20
    eval(['load rho',num2str(i),'.mat']);
    Rho=[Rho;rho];
    clear rho
end

% Plot histogram of distribution of local robustnesses. For visualization
% purposes, we omitted 0.
[N,X]=hist(nonzeros(Rho));
figure
bar(X,N/sum(N),'r')
xlabel('Fraction accepted sets')
ylabel('Relative frequency')
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