clear,close all
% Script to plot the probabilities of finding a certain network in a
% certain S-P region.
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

% Set font size and initialize data matrices.
set(0,'defaultaxesfontsize',24);
set(0,'defaulttextfontsize',24);
SPb=[];
SPg=[];
% Load all data.
for i = 1:20
    SPb=[SPb;load(['SPb',num2str(i),'.dat'])];
    SPg=[SPg;load(['SPg',num2str(i),'.dat'])];
    i
end
%%
% Load colormap (jet, but with a final value what is white instead of dark
% blue).
cmap=load('colormap.dat');
% Logarithmic increments of sensitivity and precision.
Sinc=-2.75:0.5:1.25;
Pinc=-1.25:0.5:2.75;

% figure
% Combine good and bad results.
S=[SPb(:,1);SPg(:,1)];
P=[SPb(:,2);SPg(:,2)];
% Anything under the diagonal should be projected above the the diagonal.
P(P<-S)=-S(P<-S);
% Relative histogram of sensitivities and precisions with bins defined by
% increments.
n=hist3([S(:,1),P(:,1)],{Sinc Pinc})/1e7;
% Transpose and take logarithm of relative histogram (=probabilities).
n1 = n';
n1=log10(n1);
% If we find minus infinity, set it to minimal value (for plotting).
n1(isinf(n1))=-7;
% Accommodate size.
n1( size(n,1) + 1 ,size(n,2) + 1 ) = 0;
% Define x and y-axes (log10(S) and log10(P), respectively).
sb1 = linspace(-3,1.5,size(n,1)+1);
pb1 = linspace(-1.5,3,size(n,1)+1);
% Pseudocolor plot of the relative histogram projection
figure
pcolor(sb1,pb1,n1)
% Set labels, use the correct colormap with the correct values.
xlabel('log_{10}(S)')
ylabel('log_{10}(P)')
colormap(cmap)
caxis([-7 0])
% Show colorbar and make axes square.
colorbar
axis square
set(gca,'LineWidth',3)