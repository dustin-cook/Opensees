function [ node ] = fn_update_mass( node, nodes_on_slab, direction )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


slab_node_idx = ismember(node.id,nodes_on_slab);
        
% Find Cental Coordinate and nodes to each side of it
building_length = max(node.(direction)) - min(node.(direction));
node_center = building_length/2;

nodes_right = nodes_on_slab(node.(direction)(slab_node_idx) > node_center);
nodes_left = nodes_on_slab(node.(direction)(slab_node_idx) < node_center);
mass_as_is = sum(node.mass(slab_node_idx));

nodes_right_idx = ismember(node.id,nodes_right);
nodes_left_idx = ismember(node.id,nodes_left);

% Iteratevely Increase the Mass till a target offset is reached
offset = 0;
offset_target = 0.05; % as a percentage of the building length
threshold = 0.001;
mass_delta = 0;
while abs(offset_target-offset) > threshold
    mass_delta = mass_delta + 0.001;

    % Infinite loop happen safety device
    if mass_delta > 0.5
        error('Could Not FInd Accidental Torsion Adjustment, Recalibrate')
    end

    temp_mass = node.mass;
    temp_mass(nodes_right_idx) = node.mass(nodes_right_idx)*(1+mass_delta);
    temp_mass(nodes_left_idx) = node.mass(nodes_left_idx)*(1-mass_delta);

    % Make sure mass still adds up to what it did before
    mass_new = sum(temp_mass(slab_node_idx));
    mass_diff = round(abs(mass_new - mass_as_is),2);
    if mass_diff > 0
        error('Mass got thrown off trying to find accidental torsion')
    end

    % Find new offset
    offset = (sum(temp_mass(slab_node_idx).*(node.(direction)(slab_node_idx) - node_center))/sum(temp_mass(slab_node_idx)))/building_length;
end

% Set Found Masses
node.mass(nodes_right_idx) = node.mass(nodes_right_idx)*(1+mass_delta);
node.mass(nodes_left_idx) = node.mass(nodes_left_idx)*(1-mass_delta);
end

