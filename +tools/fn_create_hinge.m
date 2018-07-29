function [ node, element, hinge ] = fn_create_hinge( node, element, hinge, node_end, ele_id, hinge_id, foundation_nodes_id )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Define new nodes to connect to springs
new_node_id = node.id(end) + 1;
old_node_id = element.(node_end)(ele_id);

% Define new node properties
node.id(new_node_id) = new_node_id;
node.x(new_node_id) = node.x(old_node_id);
node.y(new_node_id) = node.y(old_node_id);
node.z(new_node_id) = node.z(old_node_id);
node.dead_load(new_node_id) = 0;
node.live_load(new_node_id) = 0;
node.mass(new_node_id) = 0;
node.trib_area_ration(new_node_id) = 0;
node.fix(new_node_id) = node.fix(old_node_id);

% connect element to new node
element.(node_end)(ele_id) = new_node_id;

% Define fixity of foundation nodes
if sum(old_node_id == foundation_nodes_id) > 0
    node.fix(new_node_id,:) = {[1 1 1 0 0 0]};
end


% Assign hinge properties
hinge.id(hinge_id) = hinge_id;
hinge.element_id(hinge_id) = element.id(ele_id);
hinge.node_1(hinge_id) = new_node_id;
hinge.node_2(hinge_id) = old_node_id;

end

