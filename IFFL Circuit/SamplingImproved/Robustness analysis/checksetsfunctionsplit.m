function [test, S, P] = checksetsfunctionsplit(K)
% Function to perform glocal analysis check on the sampling results of the
% incoherent feedforward loop. It returns a Boolean (test) which is used to
% determine whether or not add the current parameterset K to a cell array.
% Copyright H.W.H. van Roekel, MSc. 2011.

% Invariant parameters.
par=[0.06;1;0.06;1;0.06;1;0.06;1;...
    0.06;1;1200;44;40;3.5;80;30;300;440;300;150];

% Draw a parameter sample.
u0 = K(1);
T1 = K(2);
T2 = K(3);
T3 = K(4);
exoN = K(5);
Ka_UT1 = K(6);
Ka_aT1 = K(7);
Ka_aT2 = K(8);
Ka_aT3 = K(9);
Ka_inhT2 = K(10);
Ka_YT3 = K(11);
Ka_inhT3 = K(12);

par = [0.06; 0.06; 201.5; 80; 80; 80; 189.7; 44; 54.4; 54.4; ...
    3.5; 0.55; 44; 44; 121; 145.2; 189.7; 30; 30; 30; ...
    Ka_YT3; Ka_aT1; Ka_inhT2; 144; 1200; 4.3; 4.3; Ka_aT2; Ka_aT3; 1200; ...
    0.06; 0.06; Ka_inhT3; 30.27; 144; 0.06; 4.3; 121; Ka_UT1; 1200; ...
    0.06/Ka_UT1; 0.06/Ka_aT1; 0.06/Ka_aT2; 0.06/Ka_inhT2; 0.06/Ka_aT3; 0.06/Ka_YT3; 0.06/Ka_inhT3; 0.06*Ka_aT3/Ka_inhT3; 0.06*Ka_YT3/Ka_inhT3;0.025];

x0=[T1;u0;0;T2;0;T3;0;0;0;0;0;0;0;0;0;0;0;0;0;0;exoN;0;0;0;0;0;0;0];


% Initiate options for solver and simulation.
a=[0;1];
fold_change=[1;1+a(2)];

% Set initial conditions (so turn to molar instead of log10(molar)).
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
            [ts,xs]=PIFFLFull(tspan, x0, par, 1, [1e-10, 1e-10, 10]);
            ts=ts';
            xs=xs';
            tinit=tspan(end);
            x0=xs(end,:);
            t=[t;ts];
            x=[x;xs];
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
% has a defined sensitivity and precision, so parameter set is saved.
% Otherwise, parameter set is useless and not saved.
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