function [ ] = fn_torsional_amplification( story, element )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% Calculate the actual torsional amplification multiplier
TAR_x = story.torsional_factor_x;
if sum(strcmp('torsional_factor_z',story.Properties.VariableNames)) > 0
    TAR_z = story.torsional_factor_z;
end

%% X direction
if sum(TAR_x > 1.2) > 0
    error('Torsional Amplification Exceeds Limit. Forces and displacements caused by accidental torsion shall be amplified by a factor Ax')
    
    % For linear analysis calcuate the accidental torional amplification factor
    % (this only applies to forces and displacements caused by accidental
    % torsion)
    for i = 1:length(TAR_x)
        A_x(i) = min([max([TAR_x(i)/1.2;1])^2,3]);
    end

    % Amplify Forces and Displacements
    element.Pmax = element.Pmax*max(A_x); 
    element.Vmax = element.Vmax*max(A_x);
    element.Mmax = element.Mmax*max(A_x);
    story.max_disp_x = story.max_disp_x .* A_x;
    story.max_drift_x = story.max_drift_x .* A_x;
end

%% Z direction
if exist(TAR_z,'var')
    if sum(TAR_z > 1.2) > 0
        error('Torsional Amplification Exceeds Limit. Forces and displacements caused by accidental torsion shall be amplified by a factor Az')

        % For linear analysis calcuate the accidental torional amplification factor
        % (this only applies to forces and displacements caused by accidental
        % torsion)
        for i = 1:length(TAR_z)
            A_z(i) = min([max([TAR_z(i)/1.2;1])^2,3]);
        end

        % Amplify Forces and Displacements
        element.Pmax = element.Pmax*max(A_z); 
        element.Vmax = element.Vmax*max(A_z);
        element.Mmax = element.Mmax*max(A_z);
        story.max_disp_z = story.max_disp_x .* A_z;
        story.max_drift_z = story.max_drift_x .* A_z;
    end
end

%% For 2D models:
% Linear: forces and dispalcements shall be multipled by the max value of TAR calculated for the building
% Nonlinear: the ground motion shall be amplified by the max value of TAR calculated for the building. 

%% NOTE
% It seems like a seperate analysis is needed for torsional considerations
% This amplification may happen after capacity and c factor calcs?
end

