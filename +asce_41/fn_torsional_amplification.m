function [ ] = fn_torsional_amplification( story, element )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% Calculate the actual torsional amplification multiplier
TAR = story.max_disp_x ./ story.ave_disp_x;
% if sum(TAR > 1.2) > 0
%     warning('Torsional Amplification Exceeds Limit. Forces and displacements caused by accidental torsion shall be amplified by a factor Ax')
% end

% For linear analysis calcuate the accidental torional amplification factor
% (this only applies to forces and displacements caused by accidental
% torsion)
for i = 1:length(TAR)
    A_x(i) = min([max([TAR(i)/1.2;1])^2,3]);
end

% For 2D models:
% Linear: forces and dispalcements shall be multipled by the max value of TAR calculated for the building
% Nonlinear: the ground motion shall be amplified by the max vale of TAR calculated for the building. 

%% NOTE
% It seems like a seperate analysis is needed for torsional considerations
% This amplification may happen after capacity and c factor calcs?

%% Amplify Forces and Displacements
element.Pmax = element.Pmax*max(TAR); 
element.Vmax = element.Vmax*max(TAR);
element.Mmax = element.Mmax*max(TAR);
story.max_disp_x = story.max_disp_x* max(TAR);
story.max_drift_x = story.max_drift_x*max(TAR);

% UPDATE TO WORK FOR Each DIRECTION
% story.max_disp_z = story.max_disp_z*max(TAR);
% story.max_drift_z = story.max_drift_z*max(TAR);

end

