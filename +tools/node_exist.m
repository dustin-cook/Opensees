function [ node, id ] = node_exist( node, x, y, z )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

start_node_check = (node.x == x & node.y == y & node.z == z);
if sum(start_node_check) == 0 % New Node
    node_id = max(node.id) + 1;
    node.id(node_id) = node_id;
    node.x(node_id) = x;
    node.y(node_id) = y;
    node.z(node_id) = z;
    id = node_id;
else % Existing Node
    id = node.id(start_node_check);
end

end

