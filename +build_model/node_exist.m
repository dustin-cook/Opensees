function [ node, id ] = node_exist( node, x, y, z )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isfield(node,'x') && isfield(node,'y') && isfield(node,'z')
    start_node_check = (node.x == x & node.y == y & node.z == z);
    if sum(start_node_check) == 0 % New Node
        node_id = max(node.id) + 1;
        node.id(node_id,1) = node_id;
        node.x(node_id,1) = x;
        node.y(node_id,1) = y;
        node.z(node_id,1) = z;
        id = node_id;
    else % Existing Node
        id = node.id(start_node_check);
    end
else % First Node
    node.id(1,1) = 1;
    node.x(1,1) = x;
    node.y(1,1) = y;
    node.z(1,1) = z;
    id = 1;
end

end

