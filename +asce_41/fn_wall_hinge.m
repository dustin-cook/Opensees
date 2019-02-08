function [ hinge ] = fn_wall_hinge( ele, ele_props, oop_tag )
% Description: Find wall hinge properties based on Table 10-19 of ASCE 41-17
% Created by: Dustin Cook
% Date Created: 9-21-18

%% Assmptions
% 1. Flexure controlled
% 2. Boundaries are not confined
% 3. Tension and Compression steel are equal
% 4. Shear term uses shear demand on wall
% 5. P term is axial load from gravity
% 6. Currently gravity felt by wall comes out to be very small
% 7. Assumes no coupling beams in the walls
% 8. Walls are deformation cotrolled and meet min tranverse reinforcement ratio of 0.0015

%% TODO
% - update gravity loads on walls
% - check if cracking moment exceeds yeild strength
% - shear v flexure controlled
% - splice checks
% - update strengh calc based on boundary rebar and ch 22 and ch 18 of ACI 318-14

%% Initial Setup
% Import Packages
import asce_41.*

% Define demands if they do not exist
if sum(strcmp('Pmax',ele.Properties.VariableNames)) == 0
    ele.Pmax = 0;
end
if sum(strcmp('Vmax',ele.Properties.VariableNames)) == 0
    ele.Vmax = 0;
end
if sum(strcmp('Vmax_oop',ele.Properties.VariableNames)) == 0
    ele.Vmax_oop = 0;
end

% Defined Critical Model
if oop_tag
    critical_mode = ele.critical_mode_oop;
else
    critical_mode = ele.critical_mode;
end

%% Load Wall Hinge Table 10-19 from ASCE 41-17
if strcmp(critical_mode,'flexure')
    hinge_table = readtable(['+asce_41' filesep 'wall_hinge_flexure.csv'],'ReadVariableNames',true);
    hinge_table.id = []; % Omit id 

    %% Filter Based on Confinement Boundary
    confined_boundary = 'N';
    hinge_filt = hinge_table(strcmp(hinge_table.confined_boundary,confined_boundary),:);

    %% Filter table based on strength term
    As = sum(str2double(strsplit(strrep(strrep(ele_props.As{1},']',''),'[',''))))/2;
    As_prime = As;
    strength_term = ((As - As_prime)*ele_props.fy_e + abs(ele.Pmax)) / (ele_props.w*ele_props.d*ele_props.fc_e); % The bigger the Axial load the more consertvative the hinge is, therefore the max axial load from the analysis is used.
    [ hinge_filt ] = fn_filter_asce41_table( hinge_filt, strength_term, 'strength_term', {'a_hinge','b_hinge','c_hinge','io','ls','cp'} );

    %% Filter table based on V ratio
    if oop_tag
        v_ratio = ele.Vmax_oop / (ele_props.w*ele_props.d*sqrt(ele_props.fc_e));
    else
        v_ratio = ele.Vmax / (ele_props.w*ele_props.d*sqrt(ele_props.fc_e));
    end
    [ hinge_filt ] = fn_filter_asce41_table( hinge_filt, v_ratio, 'v_ratio', {'a_hinge','b_hinge','c_hinge','io','ls','cp'} );

    %% Double Check only 1 row of the hinge table remains
    if height(hinge_filt) ~= 1
        error('Hinge table filtering failed to find unique result')
    end

    %% Save to table
    hinge = hinge_filt(:,4:9);
    
elseif strcmp(critical_mode,'shear')
    hinge_table = readtable(['+asce_41' filesep 'wall_hinge_shear.csv'],'ReadVariableNames',true);
    hinge_table.id = []; % Omit id 
    
    %% Filter table based on strength term
    As = sum(str2double(strsplit(strrep(strrep(ele_props.As{1},']',''),'[',''))))/2;
    As_prime = As;
    strength_term = ((As - As_prime)*ele_props.fy_e + abs(ele.Pmax)) / (ele_props.w*ele_props.d*ele_props.fc_e); % The bigger the Axial load the more consertvative the hinge is, therefore the max axial load from the analysis is used.
    if strength_term <= 0.05
        hinge_filt = hinge_table(strcmp(hinge_table.strength_term,'<=0.05'),:);
    else
        hinge_filt = hinge_table(strcmp(hinge_table.strength_term,'>0.05'),:);
    end
    
    %% Double Check only 1 row of the hinge table remains
    if height(hinge_filt) ~= 1
        error('Hinge table filtering failed to find unique result')
    end

    %% Save to table
    hinge = hinge_filt(:,2:9);
end

