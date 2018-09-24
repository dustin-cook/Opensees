function [ hinge ] = fn_wall_hinge( ele, ele_props )
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


%% Import Packages
import asce_41.*

%% Load Wall Hinge Table 10-19 from ASCE 41-17
hinge_table = readtable(['+asce_41' filesep 'wall_hinge.csv'],'ReadVariableNames',true);
hinge_table.id = []; % Omit id 

%% Filter Based on Confinement Boundary
confined_boundary = 'N';
hinge_filt = hinge_table(strcmp(hinge_table.confined_boundary,confined_boundary),:);
    
%% Filter table based on strength term
As = sum(str2double(strsplit(strrep(strrep(ele_props.As{1},']',''),'[',''))))/2;
As_prime = As;
strength_term = ((As - As_prime)*ele_props.fy_e + abs(ele.P_grav)) / (ele_props.w*ele_props.d*ele_props.fc_e);
[ hinge_filt ] = fn_filter_asce41_table( hinge_filt, strength_term, 'strength_term', {'a_hinge','b_hinge','c_hinge','io','ls','cp'} );

%% Filter table based on V ratio
v_ratio = ele.Vmax / (ele_props.w*ele_props.d*sqrt(ele_props.fc_e));
[ hinge_filt ] = fn_filter_asce41_table( hinge_filt, v_ratio, 'v_ratio', {'a_hinge','b_hinge','c_hinge','io','ls','cp'} );

%% Double Check only 1 row of the hinge table remains
if length(hinge_filt.a_hinge) ~= 1
    error('Hinge table filtering failed to find unique result')
end
    
%% Save to table
hinge = hinge_filt(:,4:9);
end

