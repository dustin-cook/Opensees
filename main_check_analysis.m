function [ ] = main_check_analysis( )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here


%% Torsion
% If the displacement multiplier ? caused by actual plus accidental 
% torsion at any level exceeds 1.5, two-dimensional models shall not be 
% permitted and three-dimensional models that account for the spatial 
% distribution of mass and stiffness shall be used.

% Increased forces and displacements caused by accidental torsion need not 
% be considered if either of the following conditions apply: (a) the 
% accidental torsional moment is less than 25% of the actual torsional 
% moment, or (b) the ratio of the displacement multiplier ? caused by the 
% actual plus accidental torsion and the displacement multiplier caused by 
% actual torsion is less than 1.1 at every floor.

%% Check convergence of model iterations

%% Calculate Beam Column Strength Ratios
% for i =1:length(joint.id)
%    beam1 = max([element.Mn_pos(element.node_2 == joint.x_neg(i)),element.Mn_neg(element.node_2 == joint.x_neg(i))]); % Maximum of beam pos and neg nominal bending strength
%    beam2 = max([element.Mn_pos(element.node_1 == joint.x_pos(i)),element.Mn_neg(element.node_1 == joint.x_pos(i))]); 
%    column1 = min([element.Mn_pos(element.node_2 == joint.y_neg(i)),element.Mn_neg(element.node_2 == joint.y_neg(i))]); % Minimum of column pos and negative nominal moment strength
%    column2 = min([element.Mn_pos(element.node_1 == joint.y_pos(i)),element.Mn_neg(element.node_1 == joint.y_pos(i))]); 
%    joint.beam_strength(i) = sum([beam1,beam2]);
%    joint.column_strength(i) = sum([column1,column2]);
%    joint.col_bm_ratio(i) = joint.column_strength(i)/joint.beam_strength(i);
% end

end

