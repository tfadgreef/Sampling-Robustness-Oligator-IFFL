% File to plot responses from sampling.
% 
% 
%
clear all
close all
textsize=24;
set(0, 'DefaultAxesFontSize', textsize);
set(0, 'defaulttextfontsize', textsize);
set(0, 'defaulttextinterpreter','tex');
minU=-1;
maxU=2;
minT1=-1;
maxT1=2;
minT2=-1;
maxT2=2;
minT3=-1;
maxT3=2;
minKa=5;
maxKa=12;
minE=1;
maxE=2;
% Load relevant parameter sets
parg=[];
for i=1:20
    i
    parg=[parg;load(['parg',num2str(i),'.dat'])];
end
% Indicate which sets should be plotted (here, the first 10).
for n=1:11
    % Calculate parameter values.
    samples=parg;
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
    
    % Set parameter and initial condition vector.
    par = [0.06; 0.06; 201.5; 37.296; 37.296; 37.296; 189.7; 38.15; 54.4; 54.4; ...
        3.5; 0.55; 38.15; 38.15; 121; 145.2; 189.7; 24; 24; 24; ...
        Ka_YT3; Ka_aT1; Ka_inhT2; 144; 908; 4.3; 4.3; Ka_aT2; Ka_aT3; 908; ...
        0.06; 0.06; Ka_inhT3; 30.27; 144; 0.06; 4.3; 121; Ka_UT1; 908; ...
        0.06/Ka_UT1; 0.06/Ka_aT1; 0.06/Ka_aT2; 0.06/Ka_inhT2; 0.06/Ka_aT3; 0.06/Ka_YT3; 0.06/Ka_inhT3; 0.06*Ka_aT3/Ka_inhT3; 0.06*Ka_YT3/Ka_inhT3;0.025];

    x0=[T1;u0;0;T2;0;T3;0;0;0;0;0;0;0;0;0;0;0;0;0;0;exoN;0;0;0;0;0;0;0];
    % Set inputs.
    x0init=x0;
    U0=u0;
    monitor_input=U0;
    It=[];
    a=[0;1];
    % Fold change is the ratio of new input minus old divided by old.
    fold_change=[1;1+a(2)];
    tinit=0;
    t=[];
    x=[];
    for j = 1:2
        % Initial state = free U plus a times additional U0.
        x0(2)=monitor_input+a(j)*U0;
        tspan=[tinit tinit+1000];
        % Call ode15s (delete/comment if you use mex-file).
        %     options = odeset('AbsTol', 1e-11, 'RelTol', 1e-11);
        %     [t,x]=ode15s(@(t,x)IFFLFunction(t,x,par,0),tspan,x0,options);
        % Call MEX file (delete/comment if you use ode15s,m and also delete
        % transpose.
        [ts,xs]=PIFFLFull(tspan, x0, par, 1, [1e-10, 1e-10, 10]);
        ts=ts';
        xs=xs';
        tinit=tspan(end);
        x0=xs(end,:);
        t=[t;ts];
        if j==1
            index1=length(t);
            Y1=xs(end,7);
        end
        if j==2
            index2=length(t);
            Y2=xs(end,7);
            Ymax=max(xs(:,7));
        end
        x=[x;xs];
        % If we have reached staady state after 1000 minutes, we
        % continue simulating, otherwise we stop.
        monitor_input=x(end,2);
        if j==1
            % Index value of first steady state.
            index=find(t==t(end));
        elseif j==2
            % Index value of second steady state.
            index2=find(t==t(end));
        end
    end
    % Evaluate U0 and U1 (for first and second steady state).
    U1=fold_change(1)*u0;
    U2=fold_change(2)*u0;
    It=[ones(index1,1)*U1;ones(index2-index1,1)*U2];
    S=abs((Ymax-Y1)/Y1)/abs((U2-U1)/U1);
    P=(abs((Y2-Y1)/Y1)/abs((U2-U1)/U1))^-1;
    % Generate the data to plot.
    figure
    tthres=find(t>tspan(1)-tspan(1)/10);
    tthres=tthres(1);
    t2=t(tthres);
    [AX,H1,H2] = plotyy(t((tthres:end)),It((tthres:end)),t(tthres:end),x((tthres:end),7));
    set(gca,'LineWidth',2)
    axis(AX(1),[t2(1) max(t) 0  2*ceil(It(end))])
    axis(AX(2),[t2(1) max(t) 0.75*x(end,7) 1.2*max(x(index-100:index2,7))])
    tickU=[0:2*ceil(It(end))/4:2*ceil(It(end))];
    tickY=[double(int64(0.75*x(end,7)*100))/100:double(int64((1.2*max(x(index-100:index2,7))-0.75*x(end,7))/4*100))/100:double(int64(1.2*max(x(index-100:index2,7))*100))/100];
    set(AX(1),'YTick',tickU)
    set(AX(2),'YTick',tickY)
    set(get(AX(1),'Ylabel'),'String','U [nM]','Color','b')
    set(get(AX(2),'Ylabel'),'String','Y [nM]','Color','r')
    set(AX(1),'YColor', 'b')
    set(AX(2),'YColor', 'r')
    set(H1,'color','b','LineWidth',3)
    set(H2,'color','r','LineWidth',3)
    xlabel('Time [min]')
    title(['log_{10}(S) = ',num2str(log10(S)),'; log_{10}(P) = ',num2str(log10(P)),],'FontSize',textsize)
end

