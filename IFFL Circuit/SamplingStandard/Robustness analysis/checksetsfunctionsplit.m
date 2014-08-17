function [test, S, P] = checksetsfunctionsplit(K)
% Function to perform glocal analysis check on the sampling results of the
% Oligator, for individual parameter robustnesses, ie. Ka_aT1  = Ka_aT2 and
% Ka_bT2 ~= Ka_bT3. It returns a Boolean (test) which is used to
% determine whether or not add the current parameterset K to a cell array.
% After Hafner et al. (2009) PloS Comput. Biol. 5, e1000534
% Original scripts can be found on:
% http://www.ieu.uzh.ch/wagner/publications-software.html

% Draw a parameter sample.
u0 = K(1);
T1 = K(2);
T2 = K(3);
T3 = K(4);
Ka_UT1 = K(5);
Ka_aT1 = K(6);
Ka_aT2 = K(7);
Ka_aT3 = K(8);
Ka_inhT2 = K(9);
Ka_YT3 = K(10);
Ka_inhT3 = K(11);

% Initiate parameter and initial condition vectors.
par = [0.06; 0.06; 150; 80; 80; 80; 440; 44; 54.4; 54.4; ...
    3.5; 44; 44; 300; 300; 440; 30; 30; 30; Ka_YT3; ...
    Ka_aT1; Ka_inhT2; 144; 1200; Ka_aT2; Ka_aT3; 1200; 0.06; 0.06; Ka_inhT3; ...
    40; 144; 0.06; 300; Ka_UT1; 1200; 0.06/Ka_UT1; 0.06/Ka_aT1; 0.06/Ka_aT2; 0.06/Ka_inhT2; ...
    0.06/Ka_aT3; 0.06/Ka_YT3; 0.06/Ka_inhT3; 0.06*Ka_aT3/Ka_inhT3; 0.06*Ka_YT3/Ka_inhT3];

x0=[T1;u0;0;T2;0;T3;0;0;0;0;0;0;0;0;0;0;0;0;0;0];


% Fold change is the ratio of new input minus old divided by old.
a=[0;1];
fold_change=[1;1+a(2)];

% Set inputs.
U0=u0;

% Initialize empty time vector and empty state vector to be filled
% during simulation.
t=[];
x=[];
tinit=0;
monitor_input=U0;
cont=1;
for j = 1:2
    % Initial state = free U plus a times additional U0.
    x0(2)=monitor_input+a(j)*U0;
    % Evaluate U0 and U1 (for first and second steady state).
    eval(['u',num2str(j),'=fold_change(j)*U0;']);
    % If we continue (i.e. a steady state is reached after 1000
    % minutes), then we simulate a second time too.
    if cont
        try
            tspan=[tinit tinit+1000];
            [ts,xs]=PIFFL(tspan, x0, par, 1, [1e-10, 1e-10, 10]);
            ts=ts';
            xs=xs';
            tinit=tspan(end);
            x0=xs(end,:);
            t=[t;ts];
            x=[x;xs];
            % If we have reached staady state after 1000 minutes, we
            % continue simulating, otherwise we stop; criterion is the
            % maximum percentage each state deviates at the end of the
            % simulation from its mean of the last 6 time points should be
            % less than 0.001%. Since we simulate 1000 minutes and results
            % are given each minute, this means that the system should not
            % have changed (much) for the last 6 indices.
            monitor_input=x(end,2);
            if max(abs(mean(x(end-5:end,:))-x(end,:))./x(end,:))*100<1e-3 %Steady state?
                cont=1;
            else
                cont=0;
            end
            monitor_input=x(end,2);
        catch
            t=[0:1000];
            x=zeros(1000,7);
            index=1000;
            index2=1000;
            cont=0;
            j=0;
        end
        
        if j==1 && cont
            % Index value of first steady state.
            index=find(t==t(end));
        elseif j==2 && cont
            % Index value of second steady state.
            index2=find(t==t(end));
        else
            % Index value of first part and second part are equal if we do
            % not reach steady state in either one of the parts, such that
            % parameter set will not be saved.
            index=find(t==t(end));
            index2=index;
        end
        % Evaluate final value of output.
        eval(['O',num2str(j),'=x(end,7);']);
    end
end

% If we did simulate twice a 1000 minutes, then current parameter set
% has a defined sensitivity and precision, and the test Boolean becomes one
% if the response is adaptive.
if (index~=index2) && (x(end,7)>=0.1)
    Op = max(x(index:index2,7));
    S = abs(((Op-O1)/O1)/((u2-u1)/u1));
    P = 1/abs(((O2-O1)/O1)/((u2-u1)/u1));
    if (log10(S)>=-0.5)&&(log10(P)>=1)
        test=1;
    else
        test=0;
    end
else
    S = 0;
    P = 0;
    test=0;
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