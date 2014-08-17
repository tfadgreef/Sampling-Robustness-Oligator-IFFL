% Compiling ODE to MEX file. For CVode wrapper,to compile it yourself,
% visit the home page of dr. J. Vanlier:
% http://www.fdmold.uni-freiburg.de/~jvanlier/
home=pwd;
cd 'C:\...\Documents\MATLAB\ODEMEXv11'
parser=pwd;
cd 'C:\...\Documents\MATLAB\ODEMEXv11\Parser\CVode'
setupCVode 
cd(parser)
addPaths
% 20 states
for i=1:20
    eval(['mStruct.s.x',num2str(i),'=',num2str(i),';']);
    eval(['mStruct.p.par',num2str(i),'=',num2str(i),';']);
end
% 45 parameters
for i=21:45
    eval(['mStruct.p.par',num2str(i),'=',num2str(i),';']);
end
mStruct.u.u1=1;
mStruct.c.c1=1;
cd(home)
% m-file containing ODEs
convertToC( mStruct, 'IFFL_MMFunction.m' );
compileC( 'PIFFL' );
quit
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