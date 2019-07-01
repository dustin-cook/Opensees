function [ node, idx ] = fn_node_exist( node, x, y, z, add_mass )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isfield(node,'x') && isfield(node,'y') && isfield(node,'z')
    start_node_check = (node.x == x & node.y == y & node.z == z);
    if sum(start_node_check) == 0 % New Node
        node_id = max(node.id) + 1;
        node_table_index = length(node.id) + 1;
        node.id(node_table_index,1) = node_id;
        node.x(node_table_index,1) = x;
        node.y(node_table_index,1) = y;
        node.z(node_table_index,1) = z;
        if exist('add_mass','var') && add_mass == 1
            node.mass(node_table_index,1) = 1;
        else
            node.mass(node_table_index,1) = 0;
        end
        node.record_disp(node_table_index,1) = 0;
        node.record_accel(node_table_index,1) =  0;
        idx = node_table_index;
    else % Existing Node
        idx = start_node_check;
    end
else % First Node
    node.id(1,1) = 6000;
    node.x(1,1) = x;
    node.y(1,1) = y;
    node.z(1,1) = z;
    node.mass(1,1) = 0;
    node.record_disp(1,1) = 0;
    node.record_accel(1,1) =  0;
    idx = 1;
end

end

