set(0, 'DefaultAxesFontSize', 24);
set(0, 'defaulttextfontsize', 24);
set(0, 'defaulttextinterpreter','tex');
clear all
% File to make boxplots of the sampling results, depicting desired behavior
% against sampling range.
% Number of parameters.
Npar=8;
par1=[];
par2=[];
% Load parameter samples (values between 0 and 1) corresponding to
% sustained oscillations in the Oligator.
for i = 1:20
    par1=[par1;load(['parS',num2str(i),'.dat'])];
end
close all
% Give sampling ranges.
parZ=[par1];
clear par1
mina = -1;
maxa = 1;
minT1 = 0;
maxT1 = 2;
minT2 = 0;
maxT2 = 2;
minT3 = 0;
maxT3 = 2;
minKa = 5;
maxKa = 12;

% Calculate the corresponding parameter values in nM and add them to the
% array par.
for i=1:1
    j=(i-1)*Npar;
    par(:,j+1)=(mina+(maxa-mina)*parZ(:,j+1));
    par(:,j+2)=(minT1+(maxT1-minT1)*parZ(:,j+2));
    par(:,j+3)=(minT2+(maxT2-minT2)*parZ(:,j+3));
    par(:,j+4)=(minT3+(maxT3-minT3)*parZ(:,j+4));
    par(:,j+5)=-(minKa+(maxKa-minKa)*parZ(:,j+5))+9;
    par(:,j+6)=-(minKa+(maxKa-minKa)*parZ(:,j+6))+9;
    par(:,j+8)=-(minKa+(maxKa-minKa)*parZ(:,j+8))+9;
    par(:,j+7)=-log10((1+(5-1)*parZ(:,j+7)).*10.^(-par(:,j+8)));
end

% Set variables for boxplot.
par1=par;
species1 = repmat({'[a]0','[T1]','[T2]','[T3]','Kda','Kdb','KdinhT3','KdinhT1'},1);

% Define spaces to group with species in the boxplot.
firstlast = [repmat({' '},1,8)];
% Show the boxplots.
figure
h=boxplot(par1,{species1,firstlast},'colors',repmat('r',1,8),'factorgap',[20 2],'labelverbosity','minor','widths',1,'plotstyle','compact','jitter',0.25);
set(h(1,:),'LineWidth',3)
set(h(2,:),'LineWidth',10)
set(h(3,:),'MarkerSize',20)
set(h(4,:),'MarkerSize',20)

set(gca,'XTickLabel','')
Xmax=get(gca,'XLim');

% Set ticks.
ticklocations=[1: floor(100*(Xmax(2)+0.2)/8)/100:Xmax(2)];
set(gca,'LineWidth',2,'XTick',ticklocations,'XTickLabel',{'[\alpha]_0','[T1]','[T2]','[T3]','K_d^{\alpha}','K_d^{\beta}','K_d^{inhT3}','K_d^{inhT1}'}) 
hold on
X=get(gca,'XTick');

% Rotate ticks.
XTickLabel = get(gca,'XTickLabel');
set(gca,'XTickLabel',' ');
hxLabel = get(gca,'XLabel');
set(hxLabel,'Units','data');
xLabelPosition = get(hxLabel,'Position');
y = xLabelPosition(2);
XTick = get(gca,'XTick');
y=repmat(y,length(XTick),1);
fs = get(gca,'fontsize');
hText = text(XTick,y,XTickLabel,'fontsize',fs);
set(hText,'Rotation',40,'HorizontalAlignment','center','VerticalAlignment','bottom');

% Draw the sampling ranges in gray.
x=[X(1) X(1)];
y1=[mina maxa];
y2=[minT1 maxT1];
y3=[minT2 maxT2];
y4=[minT3 maxT3];
y5=-[minKa maxKa]+9;
y6=y5;
y7=[-minKa -log10(5*10^maxKa)]+9;
y8=y5;
y9=y8;
plot(x,y1,'LineWidth',10,'color',[.5 .5 .5])

x=[X(2) X(2)];
plot(x,y2,'LineWidth',10,'color',[.5 .5 .5])

x=[X(3) X(3)];
plot(x,y3,'LineWidth',10,'color',[.5 .5 .5])

x=[X(4) X(4)];
plot(x,y4,'LineWidth',10,'color',[.5 .5 .5])

x=[X(5) X(5)];
plot(x,y5,'LineWidth',10,'color',[.5 .5 .5])

x=[X(6) X(6)];
plot(x,y6,'LineWidth',10,'color',[.5 .5 .5])

x=[X(7) X(7)];
plot(x,y7,'LineWidth',10,'color',[.5 .5 .5])

x=[X(8) X(8)];
plot(x,y8,'LineWidth',10,'color',[.5 .5 .5])

% Send ranges (gray bars) to the back.
chH = get(gca,'Children');
set(gca,'Children',[chH(end);chH(1:end-1)])
ylabel('log_{10}(ranges)')

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