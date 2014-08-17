function K=plotCorrelations(K,S,A,B,string,c)

% K is matrix containing all dissociation constants: determine size.
N=size(K,2);

% Calculate deltaG values or just take the logarithm of the
% concentrations/dissociation constants. In any case, the variable with
% well be working is deltaG.

K=log10(K);

% Make the scatterplot
figure
[H,AX]=gplotmatrix(K,[],S,string,'o',1,'on','hist');

% Define the axes ticks
tickX=[round(A):round((B-A)/2.6):round(B)];
tickY=[round(A):round((B-A)/2.6):round(B)];

n=[];
for i=1:N
    eval(sprintf('n%s=hist(K(:,%s))', num2str(i),num2str(i)));
    eval(sprintf('n=[n,n%s]', num2str(i)));
end

% Change the plot to match axis limits.
for i = 1:N
    for j = 1:N % Change axes limits
        if i~=j
            axis(AX(i,j),[A B A B])
            for k = 1:length(string)
                set(H(i,j,k),'MarkerFaceColor',string(k)) % H containts plot info.
            end% Change marker settings
        else
            axis(AX(i,j),[A B A B])
            eval(sprintf('axis(AX(%s,j),[A B 0 max(n)])',num2str(N+1))) % Extra axis for histograms on the diagonal of the plot.
        end
        if j==1
            set(AX(i,j),'YTick',tickY)
        end
        if i==5
            set(AX(i,j),'XTick',tickX)
        end
    end
end
if c
    xlabel('log_{10}(K_d /nM)')
    ylabel('log_{10}(K_d /nM)')
else
    xlabel('log_{10}(Concentration /nM)')
    ylabel('log_{10}(Concentration /nM)')
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