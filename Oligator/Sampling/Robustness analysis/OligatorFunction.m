function dx = OligatorFunction(t,x,par,u)
% Full set of differential equations of the Oligator
% If you're going to solve this system with a Matlab ode solver, I
% recommend ode15s. But for faster results, I recommend you use compiled
% MEX files using numerical integrators from the SUNDIALS CVode package. I
% have included the .mexa64 and .mexw64 files, which provide solutions in
% Ubuntu and Windows on 64-bit systems, respectively.
% For the software (CVode wrapper) to compile it yourself, I refer you to
% dr. J. Vanlier: http://www.fdmold.uni-freiburg.de/~jvanlier/
%
% Rik van Roekel, 2014
% After Montagne et al. (2011), Mol.Syst.Biol. 7, 466
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
% x(2) = [a] 
% x(3) = [T2] 
% x(4) = [T3] 
% x(5) = [a.T1] 
% x(6) = [T1.a] 
% x(7) = [a.T2] 
% x(8) = [a.T1.a] 
% x(9) = [aa.T1] 
% x(10) = [ab.T2] 
% x(11) = [a.T2.b] 
% x(12) = [b] 
% x(13) = [T2.b] 
% x(14) = [b.T3] 
% x(15) = [bi.T3] 
% x(16) = [b.T3.i] 
% x(17) = [i] 
% x(18) = [T3.i] 
% x(19) = [i.T1] 

% par is a column vector and it represents:
% par(1) = Vm_pol_b 
% par(2) = kaTH_i 
% par(3) = Ka_b_T2 
% par(4) = Vm_nick_T1 
% par(5) = Vm_nick_T2 
% par(6) = Km_polSD_i 
% par(7) = Km_pol_i 
% par(8) = Km_exoN_b 
% par(9) = Km_exoN_a 
% par(10) = Km_polSD_b 
% par(11) = Km_polSD_a 
% par(12) = Km_pol_a 
% par(13) = Km_exoN_i 
% par(14) = Km_pol_b 
% par(15) = Ka_i_T1 
% par(16) = Ka_i_T3 
% par(17) = Vm_exoN_b 
% par(18) = Vm_exoN_a 
% par(19) = Km_nick_T3 
% par(20) = Km_nick_T2 
% par(21) = Km_nick_T1 
% par(22) = Ka_a_T1 
% par(23) = Ka_a_T2 
% par(24) = Vm_pol_a 
% par(25) = Vm_pol_i 
% par(26) = ka_i 
% par(27) = Vm_polSD_b 
% par(28) = ka_a 
% par(29) = ka_b 
% par(30) = Vm_polSD_i 
% par(31) = Vm_polSD_a 
% par(32) = Vm_exoN_i 
% par(33) = Vm_nick_T3 
% par(34) = Ka_b_T3 
% par(35) = kd_a_T1 
% par(36) = kd_a_T2 
% par(37) = kd_b_T2 
% par(38) = kd_b_T3 
% par(39) = kd_i_T3 
% par(40) = kd_i_T1 
% par(41) = kd_i_a_T1 

% Differential equations
dx(1) = -2*par(28)*x(1)*x(2) + par(35)*x(5) + par(40)*x(19) + -par(26)*x(1)*x(17) + par(35)*x(6);
dx(2) = -2*par(28)*x(1)*x(2) + par(36)*x(7) + -par(18)*x(2)/(par(9)*(1 + x(2)/par(9) + x(12)/par(8) + x(17)/par(13))) + -par(28)*x(2)*x(3) + par(35)*x(5) + par(31)*x(8)/(par(11)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + par(2)*x(17)*x(5) + 2*par(35)*x(8) + -2*par(41)*x(19)*x(2) + -par(28)*x(2)*x(5) + par(2)*x(17)*x(6) + -par(28)*x(13)*x(2) + par(36)*x(11) + -par(28)*x(2)*x(6) + par(35)*x(6);
dx(3) = par(36)*x(7) + -par(28)*x(2)*x(3) + -par(29)*x(12)*x(3) + par(37)*x(13);
dx(4) = par(39)*x(18) + par(38)*x(14) + -par(26)*x(17)*x(4) + -par(29)*x(12)*x(4);
dx(5) = -par(2)*x(17)*x(5) + par(41)*x(19)*x(2) + -par(28)*x(2)*x(5) + -par(24)*x(5)/(par(12)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + par(28)*x(1)*x(2) + -par(35)*x(5) + par(35)*x(8);
dx(6) = -par(2)*x(17)*x(6) + -par(35)*x(6) + par(41)*x(19)*x(2) + par(28)*x(1)*x(2) + -par(28)*x(2)*x(6) + par(35)*x(8);
dx(7) = -par(29)*x(12)*x(7) + par(28)*x(2)*x(3) + -par(36)*x(7) + -par(1)*x(7)/(par(14)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + par(37)*x(11);
dx(8) = -2*par(35)*x(8) + par(28)*x(2)*x(6) + par(28)*x(2)*x(5) + -par(31)*x(8)/(par(11)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + par(4)*x(9)/(par(21)*(1 + x(9)/par(21) + x(10)/par(20) + x(15)/par(19)));
dx(9) = -par(4)*x(9)/(par(21)*(1 + x(9)/par(21) + x(10)/par(20) + x(15)/par(19))) + par(31)*x(8)/(par(11)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + par(24)*x(5)/(par(12)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10)));
dx(10)= par(1)*x(7)/(par(14)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + par(27)*x(11)/(par(10)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + -par(5)*x(10)/(par(20)*(1 + x(9)/par(21) + x(10)/par(20) + x(15)/par(19)));
dx(11)= par(29)*x(12)*x(7) + -par(36)*x(11) + -par(37)*x(11) + -par(27)*x(11)/(par(10)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + par(28)*x(13)*x(2) + par(5)*x(10)/(par(20)*(1 + x(9)/par(21) + x(10)/par(20) + x(15)/par(19)));
dx(12)= par(27)*x(11)/(par(10)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + par(38)*x(14) + -par(29)*x(12)*x(3) + -par(29)*x(12)*x(7) + -par(29)*x(12)*x(18) + par(37)*x(13) + -par(29)*x(12)*x(4) + -par(17)*x(12)/(par(8)*(1 + x(2)/par(9) + x(12)/par(8) + x(17)/par(13))) + par(37)*x(11) + par(38)*x(16);
dx(13)= par(29)*x(12)*x(3) + -par(37)*x(13) + -par(28)*x(13)*x(2) + par(36)*x(11);
dx(14)= -par(25)*x(14)/(par(7)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + par(39)*x(16) + -par(26)*x(14)*x(17) + par(29)*x(12)*x(4) + -par(38)*x(14);
dx(15)= par(25)*x(14)/(par(7)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + -par(33)*x(15)/(par(19)*(1 + x(9)/par(21) + x(10)/par(20) + x(15)/par(19))) + par(30)*x(16)/(par(6)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10)));
dx(16)= par(29)*x(12)*x(18) + -par(39)*x(16) + par(33)*x(15)/(par(19)*(1 + x(9)/par(21) + x(10)/par(20) + x(15)/par(19))) + -par(38)*x(16) + -par(30)*x(16)/(par(6)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + par(26)*x(14)*x(17);
dx(17)= par(39)*x(18) + par(39)*x(16) + 2*par(41)*x(19)*x(2) + -par(26)*x(14)*x(17) + -par(2)*x(17)*x(6) + -par(2)*x(17)*x(5) + -par(26)*x(17)*x(4) + par(30)*x(16)/(par(6)*(1 + x(14)/par(7) + x(16)/par(6) + x(7)/par(14) + x(5)/par(12) + x(8)/par(11) + x(11)/par(10))) + par(40)*x(19) + -par(26)*x(1)*x(17) + -par(32)*x(17)/(par(13)*(1 + x(2)/par(9) + x(12)/par(8) + x(17)/par(13)));
dx(18)= -par(29)*x(12)*x(18) + -par(39)*x(18) + par(26)*x(17)*x(4) + par(38)*x(16);
dx(19)= -par(40)*x(19) + par(2)*x(17)*x(5) + -2*par(41)*x(19)*x(2) + par(26)*x(1)*x(17) + par(2)*x(17)*x(6);
dx = dx(:);

