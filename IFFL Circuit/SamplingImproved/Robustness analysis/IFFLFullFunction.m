function dx = IFFLFullFunction(t,x,par,u)
% Full set of differential equations of the improved IFFL circuit.
% If you're going to solve this system with a Matlab ode solver, I
% recommend ode15s. But for faster results, I recommend you use compiled
% MEX files using numerical integrators from the SUNDIALS CVode package. I
% have included the .mexa64 and .mexw64 files, which provide solutions in
% Ubuntu and Windows on 64-bit systems, respectively.
% For the software (CVode wrapper) to compile it yourself, I refer you to
% dr. J. Vanlier: http://www.fdmold.uni-freiburg.de/~jvanlier/
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

% x is a column vector and it represents:
% x(1) = [T1] 
% x(2) = [U] 
% x(3) = [a] 
% x(4) = [T2] 
% x(5) = [g] 
% x(6) = [T3] 
% x(7) = [Y] 
% x(8) = [U.T1] 
% x(9) = [T1.a] 
% x(10) = [a.T2] 
% x(11) = [T2.g] 
% x(12) = [a.T3] 
% x(13) = [T3.Y] 
% x(14) = [g.T3] 
% x(15) = [U.T1.a] 
% x(16) = [Ua.T1] 
% x(17) = [a.T2.g] 
% x(18) = [ag.T2] 
% x(19) = [Y.T3.a] 
% x(20) = [Ya.T3]
% x(21) = [exoN]
% x(22) = [exoN.T1]
% x(23) = [exoN.T2]
% x(24) = [exoN.T3]
% x(25) = [exoN.U]
% x(26) = [exoN.a]
% x(27) = [exoN.inh]
% x(28) = [exoN.Y]

% par is a column vector and it represents:
% par(1) = kaTH_g 
% par(2) = ka_a 
% par(3) = Km_exoN_g 
% par(4) = Vm_nick_T1 
% par(5) = Vm_nick_T2 
% par(6) = Vm_nick_T3 
% par(7) = Km_exoN_a 
% par(8) = Km_pol_g 
% par(9) = Km_polSD_a 
% par(10) = Km_pol_a 
% par(11) = Km_polSD_g 
% par(12) = Km_exoN_U 
% par(13) = Km_polSD_Y 
% par(14) = Km_pol_Y 
% par(15) = Vm_exoN_a 
% par(16) = Vm_exoN_g 
% par(17) = Km_exoN_Y 
% par(18) = Km_nick_T3 
% par(19) = Km_nick_T2 
% par(20) = Km_nick_T1 
% par(21) = Ka_Y_T3 
% par(22) = Ka_a_T1 
% par(23) = Ka_g_T2 
% par(24) = Vm_pol_a 
% par(25) = Vm_pol_g 
% par(26) = Km_exoN_T2 
% par(27) = Km_exoN_T3 
% par(28) = Ka_a_T2 
% par(29) = Ka_a_T3 
% par(30) = Vm_polSD_Y 
% par(31) = ka_g 
% par(32) = ka_Y 
% par(33) = Ka_g_T3 
% par(34) = Vm_polSD_g 
% par(35) = Vm_polSD_a 
% par(36) = ka_U 
% par(37) = Km_exoN_T1 
% par(38) = Vm_exoN_Y 
% par(39) = Ka_U_T1 
% par(40) = Vm_pol_Y 
% par(41) = kd_U_T1 
% par(42) = kd_a_T1 
% par(43) = kd_a_T2 
% par(44) = kd_g_T2 
% par(45) = kd_a_T3 
% par(46) = kd_Y_T3 
% par(47) = kd_g_T3 
% par(48) = kd_g_a_T3 
% par(49) = kd_g_Y_T3 
% par(50) = kf_exoN

% Differential equations
dx(1) = par(41)*x(8) + par(42)*x(9) + -par(2)*x(1)*x(3) + -par(36)*x(1)*x(2)-par(50)*x(21)*x(1)+par(50)*par(37)*x(22);
dx(2) = par(41)*x(8) + -par(36)*x(2)*x(9) + par(41)*x(15) + -par(36)*x(1)*x(2)-par(50)*x(21)*x(2)+par(50)*par(12)*x(25);
dx(3) = par(42)*x(15) + par(42)*x(9) + -par(2)*x(3)*x(4) + par(35)*x(15)/(par(9)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10))) + -par(48)*x(14)*x(3) + -par(2)*x(1)*x(3) + par(43)*x(17) + -par(2)*x(3)*x(6) + -par(2)*x(11)*x(3) + par(45)*x(12) + -par(2)*x(13)*x(3) + par(43)*x(10) + -par(2)*x(3)*x(8) + par(45)*x(19) + par(1)*x(12)*x(5)-par(50)*x(21)*x(3)+(par(50)*par(7)-par(15)/50)*x(26);
dx(4) = -par(2)*x(3)*x(4) + par(43)*x(10) + par(44)*x(11) + -par(31)*x(4)*x(5)-par(50)*x(21)*x(4)+par(50)*par(26)*x(23);
dx(5) = -par(1)*x(12)*x(5) + -par(31)*x(10)*x(5) + par(44)*x(17) + par(34)*x(17)/(par(11)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10))) + -par(31)*x(5)*x(6) + par(49)*x(14)*x(7) + par(47)*x(14) + -par(1)*x(13)*x(5) + par(44)*x(11) + par(48)*x(14)*x(3) + -par(31)*x(4)*x(5)-par(50)*x(21)*x(5)+(par(50)*par(3)-par(16)/50)*x(27);
dx(6) = -par(31)*x(5)*x(6) + par(46)*x(13) + -par(32)*x(6)*x(7) + -par(2)*x(3)*x(6) + par(45)*x(12) + par(47)*x(14)-par(50)*x(21)*x(6)+par(50)*par(27)*x(24);
dx(7) = par(46)*x(19) + par(1)*x(13)*x(5) + -par(49)*x(14)*x(7) + par(46)*x(13) + -par(32)*x(6)*x(7) + -par(32)*x(12)*x(7) + par(30)*x(19)/(par(13)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10)))-par(50)*x(21)*x(7)+(par(50)*par(17)-par(38)/50)*x(28);
dx(8) = par(36)*x(1)*x(2) + par(42)*x(15) + -par(24)*x(8)/(par(10)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10))) + -par(2)*x(3)*x(8) + -par(41)*x(8);
dx(9) = par(2)*x(1)*x(3) + -par(36)*x(2)*x(9) + par(41)*x(15) + -par(42)*x(9);
dx(10)= -par(31)*x(10)*x(5) + par(44)*x(17) + -par(43)*x(10) + -par(25)*x(10)/(par(8)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10))) + par(2)*x(3)*x(4);
dx(11)= par(43)*x(17) + par(31)*x(4)*x(5) + -par(2)*x(11)*x(3) + -par(44)*x(11);
dx(12)= par(46)*x(19) + -par(1)*x(12)*x(5) + par(2)*x(3)*x(6) + -par(32)*x(12)*x(7) + par(48)*x(14)*x(3) + -par(40)*x(12)/(par(14)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10))) + -par(45)*x(12);
dx(13)= -par(46)*x(13) + par(32)*x(6)*x(7) + par(49)*x(14)*x(7) + -par(2)*x(13)*x(3) + -par(1)*x(13)*x(5) + par(45)*x(19);
dx(14)= par(31)*x(5)*x(6) + par(1)*x(13)*x(5) + -par(49)*x(14)*x(7) + -par(48)*x(14)*x(3) + -par(47)*x(14) + par(1)*x(12)*x(5);
dx(15)= -par(35)*x(15)/(par(9)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10))) + par(2)*x(3)*x(8) + par(36)*x(2)*x(9) + par(4)*x(16)/(par(20)*(1 + x(16)/par(20) + x(18)/par(19) + x(20)/par(18))) + -par(42)*x(15) + -par(41)*x(15);
dx(16)= par(35)*x(15)/(par(9)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10))) + -par(4)*x(16)/(par(20)*(1 + x(16)/par(20) + x(18)/par(19) + x(20)/par(18))) + par(24)*x(8)/(par(10)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10)));
dx(17)= par(2)*x(11)*x(3) + -par(44)*x(17) + -par(43)*x(17) + -par(34)*x(17)/(par(11)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10))) + par(5)*x(18)/(par(19)*(1 + x(16)/par(20) + x(18)/par(19) + x(20)/par(18))) + par(31)*x(10)*x(5);
dx(18)= -par(5)*x(18)/(par(19)*(1 + x(16)/par(20) + x(18)/par(19) + x(20)/par(18))) + par(34)*x(17)/(par(11)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10))) + par(25)*x(10)/(par(8)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10)));
dx(19)= -par(45)*x(19) + par(32)*x(12)*x(7) + par(6)*x(20)/(par(18)*(1 + x(16)/par(20) + x(18)/par(19) + x(20)/par(18))) + -par(46)*x(19) + -par(30)*x(19)/(par(13)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10))) + par(2)*x(13)*x(3);
dx(20)= par(40)*x(12)/(par(14)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10))) + -par(6)*x(20)/(par(18)*(1 + x(16)/par(20) + x(18)/par(19) + x(20)/par(18))) + par(30)*x(19)/(par(13)*(1 + x(15)/par(9) + x(10)/par(8) + x(12)/par(14) + x(19)/par(13) + x(17)/par(11) + x(8)/par(10)));
dx(21)=-par(50)*x(21)*(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7))+par(50)*par(37)*x(22)+par(50)*par(26)*x(23)+par(50)*par(27)*x(24)+par(50)*par(12)*x(25)+par(50)*par(7)*x(26)+par(50)*par(3)*x(27)+par(50)*par(17)*x(28);
dx(22)=par(50)*x(21)*x(1)-par(50)*par(37)*x(22);
dx(23)=par(50)*x(21)*x(4)-par(50)*par(26)*x(23);
dx(24)=par(50)*x(21)*x(6)-par(50)*par(27)*x(24);
dx(25)=par(50)*x(21)*x(2)-par(50)*par(12)*x(25);
dx(26)=par(50)*x(21)*x(3)-par(50)*par(7)*x(26);
dx(27)=par(50)*x(21)*x(5)-par(50)*par(3)*x(27);
dx(28)=par(50)*x(21)*x(7)-par(50)*par(17)*x(28);

