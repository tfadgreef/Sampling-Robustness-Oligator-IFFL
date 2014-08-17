function [parg,SPg,parb,SPb,parNo]=samplePIFFLFull(samples,minU,maxU,minT1,maxT1,minT2,maxT2,minT3,maxT3,minKa,maxKa,minE,maxE,n1,n2);
% Function to perform sampling on the IFFL circuit.
% Input: samples matrix, minimum and maximum concentrations of protected
% input ssDNA U, T1, T2, T3, exonuclease, 
% respectively, minumum and maximum association constant values,
% and the begin (n1) and end (n2) of the portion of samples to be
% calculated.

% Initiate output arrays: good parameters, sensitivity and presicion values
% corresponding to them, bad parameters, sensitivity and precision values
% and no response, respectively.
parg=[];
SPg=[];
parb=[];
SPb=[];
parNo=[];

tinit = 0;
a=[0;1];
% Fold change is the ratio of new input minus old divided by old.
fold_change=[1;1+a(2)];

for n = n1:n2
    % From the samples, calculate the parameters.
    u0=10^(minU+(maxU-minU)*samples(n,1));
    T1=10^(minT1+(maxT1-minT1)*samples(n,2));
    T2=10^(minT2+(maxT2-minT2)*samples(n,3));
    T3=10^(minT3+(maxT3-minT3)*samples(n,4));
    Ka_UT1=10^(minKa+(maxKa-minKa)*samples(n,5))/10^9;
    Ka_aT1=10^(minKa+(maxKa-minKa)*samples(n,6))/10^9;
    Ka_aT2=Ka_aT1;
    Ka_aT3=Ka_aT1;
    Ka_inhT3=10^(minKa+(maxKa-minKa)*samples(n,8))/10^9;
    Ka_inhT2=(1+(5-1)*samples(n,7))*Ka_inhT3;
    Ka_YT3=10^(minKa+(maxKa-minKa)*samples(n,9))/10^9;
    exoN=10^(minE+(maxE-minE)*samples(n,10));
    
    % Parameter and initial cindition vectors.
    par = [0.06; 0.06; 201.5; 80; 80; 80; 189.7; 44; 54.4; 54.4; ...
        3.5; 0.55; 44; 44; 121; 145.2; 189.7; 30; 30; 30; ...
        Ka_YT3; Ka_aT1; Ka_inhT2; 144; 1200; 4.3; 4.3; Ka_aT2; Ka_aT3; 1200; ...
        0.06; 0.06; Ka_inhT3; 30.27; 144; 0.06; 4.3; 121; Ka_UT1; 1200; ...
        0.06/Ka_UT1; 0.06/Ka_aT1; 0.06/Ka_aT2; 0.06/Ka_inhT2; 0.06/Ka_aT3; 0.06/Ka_YT3; 0.06/Ka_inhT3; 0.06*Ka_aT3/Ka_inhT3; 0.06*Ka_YT3/Ka_inhT3;0.025];
    
    x0=[T1;u0;0;T2;0;T3;0;0;0;0;0;0;0;0;0;0;0;0;0;0;exoN;0;0;0;0;0;0;0];
    t=[];
    x=[];
    i=1;
    tinit=0;
    % Set inputs.
    U0=u0;
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
                [ts,xs]=PIFFLFull(tspan, x0, par, 1, [1e-11, 1e-11, 10]);
                ts=ts';
                xs=xs';
                tinit=tspan(end);
                x0=xs(end,:);
                t=[t;ts];
                x=[x;xs];
                monitor_input=x(end,2);
                % If we have reached staady state after 1000 minutes, we
                % continue simulating, otherwise we stop; criterion is the
                % maximum percentage each state deviates at the end of the
                % simulation from its mean of the last 6 time points should be
                % less than 0.001%. Since we simulate 1000 minutes and results
                % are given each minute, this means that the system should not
                % have changed (much) for the last 6 indices.
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
    
    % If we did simulate twice a 1000 minutes, then current parameter set
    % has a defined sensitivity and precision, so parameter set is saved.
    % Otherwise, parameter set is useless and not saved.
    if (index~=index2) && (x(end,7)>=2)
        Op = max(x(index:index2,7));
        S = abs(((Op-O1)/O1)/((u2-u1)/u1));
        P = 1/abs(((O2-O1)/O1)/((u2-u1)/u1));
        if (log10(S)>=-0.5)&&(log10(P)>=1)
            parg = [parg;samples(n,:)];
            SPg = [SPg;log10(S) log10(P)];
        else
            parb = [parb;samples(n,:)];
            SPb = [SPb;log10(S) log10(P)];
        end
    else
        parNo = [parNo;samples(n,:)];
    end
    
    if mod(n,5000)==0
        n
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