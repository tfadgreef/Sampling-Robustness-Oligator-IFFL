function [test] = checksetsfunctionsplit(K)
% Function to perform glocal analysis check on the sampling results of the
% Oligator, for individual parameter robustnesses, ie. Ka_aT1  = Ka_aT2 and
% Ka_bT2 ~= Ka_bT3. It returns a Boolean (test) which is used to
% determine whether or not add the current parameterset K to a cell array.
% After Hafner et al. (2009) PloS Comput. Biol. 5, e1000534
% Original scripts can be found on:
% http://www.ieu.uzh.ch/wagner/publications-software.html

% Draw a parameter sample.
a = K(1);
T1 = K(2);
T2 = K(3);
T3 = K(4);
Ka_aT1 = K(5);
Ka_aT2 = K(6);
Ka_bT2 = K(7);
Ka_bT3 = K(8);
Ka_inhT3 = K(10);
Ka_inhT1 = K(9);

% Initiate parameter and initial condition vectors.
par = [1200; 0.06; Ka_bT2; 80; 80; 3.5; 44; 440; 440; 44; ...
       44; 44; 150; 44; Ka_inhT1; Ka_inhT3; 300; 300; 30; 30; ...
       30; Ka_aT1; Ka_aT2; 1200; 1200; 0.06; 1200; 0.06; 0.06; 40; ...
       1200; 300; 80; Ka_bT3; 0.06/Ka_aT1; 0.06/Ka_aT2; 0.06/Ka_bT2; ...
        0.06/Ka_bT3; 0.06/Ka_inhT3; 0.06/Ka_inhT1; 0.06*Ka_aT1/Ka_inhT1];


x0=[T1;a;T2;T3;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

% Initialize empty time vector and empty state vector to be filled
% during simulation.
tinit = 0;
tspan=[tinit tinit+7000];
try
    [t,x]=Oligator(tspan, x0, par, 0, [1e-11, 1e-11, 10]);
    t=t';
    x=x';
catch
    x=zeros(1e3,19);
    t=[1:1:1e3];
end
signal=x(:,5)+x(:,6)+x(:,7)+x(:,8)+x(:,9)+x(:,10)+x(:,11)+x(:,13)+x(:,14)+x(:,15)+x(:,16)+x(:,18)+x(:,19);
%plot(t,signal)
[x2min,x2max,base]=getOscillations(t,signal,0.05);

maxVal=max(signal);
timesMax=t(signal>0.9*max(signal));
timesMin=t(signal<0.1*max(signal));
if numel(timesMin)>0 && numel(timesMax)>0
    timesMin=timesMin(timesMin>timesMax(1));
    
        % If we see at least 1 nM response in dsDNA within 100 minutes,
        % then the response is
        % significant. If it is significant, it can be classified as
        % oscillating (sustained or damped), transient, or sustained, or
        % something else (significant but not sustained, not transient, not
        % oscillatory etc.) We are dealing with oscillations when we have
        % at least 15 peaks deviating 5% above the baseline (average).
        % If this 15th peak is larger than 97.5% times the 14th peak,
        % we are dealing with a sustained oscillation, otherwise,
        %with a damped oscillation.
    if max(signal)>=1 && timesMax(1)<100
        if length(x2max)>=15 && length(x2min)>=15
            if (signal(x2max(15))-base) >= 0.975*(signal(x2max(14))-base)
                % Sustained oscillations
                test=1;
            else
                % Damped oscillations
                test=0;
            end
        else
            test=0;
        end
        % If there is not a 5 muM change in Tr, we define the system as
        % non-reponsive.
    else
        test=0;
    end
else
    test=0;
end
end

function [indmin, indmax, ymean] = getOscillations(t,y, Tolosc)
Ndata = length(y);
cntmax=0;
cntmin=0;
indmax=[];
indmin=[];
sgn = sign(y(2)-y(1)); % initial slope
% Calculate the mean value
ymean=trapz(t,y)/t(end);

for i=2:Ndata-1
    if(sign(y(i+1)-y(i))>0)                         % POSITIVE SLOPE
        if(sgn<0)                                    % Minimum found
            cntmin = cntmin+1;
            sgn=1;
            if( abs((y(i)-ymean)/ymean)>=Tolosc )&&y(i)<ymean&&(y(i+1)-y(i))>1e-11 % Minimum significant?
                indmin = [indmin i];
            end
        end
        
    elseif(sign(y(i+1)-y(i))<0)                    % NEGATIVE SLOPE
        if(sgn>0)                                    % Maximum found
            cntmax = cntmax+1;
            sgn=-1;
            if( abs((y(i)-ymean)/ymean)>=Tolosc )&&y(i)>ymean&&(y(i+1)-y(i))<-1e-11  % Maximum significant?
                indmax = [indmax i];
            end
        end
    end
end
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