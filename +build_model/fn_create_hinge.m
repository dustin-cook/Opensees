function [ node, element, hinge ] = fn_create_hinge( node, element, hinge, node_end, ele_or_node_id, hinge_id, foundation_nodes_id, type, direction )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Define new nodes to connect to springs
new_node_id = node.id(end) + 1;
new_node_idx = length(node.id) + 1;
if strcmp(type,'foundation')
    old_node_id = ele_or_node_id;
else
    old_node_id = element.(node_end)(ele_or_node_id);
end
old_node_idx = find(node.id == old_node_id);

% Define new node properties
node.id(new_node_idx) = new_node_id;
node.x(new_node_idx) = node.x(old_node_idx);
node.y(new_node_idx) = node.y(old_node_idx);
node.z(new_node_idx) = node.z(old_node_idx);
node.dead_load(new_node_idx) = 0;
node.live_load(new_node_idx) = 0;
node.mass(new_node_idx) = 0;
node.record_disp(new_node_idx) = 0;
node.record_accel(new_node_idx) = 0;
node.story(new_node_idx) = 0;
node.primary_story(new_node_idx) = 0;
node.fix(new_node_idx) = node.fix(old_node_idx);
node.on_slab(new_node_idx) = 0;

% Define fixity of foundation nodes
if foundation_nodes_id(old_node_idx)
    if strcmp(type,'foundation')
        node.fix{new_node_idx} = '[111111]';
    else
        node.fix{new_node_idx} = '[000000]';
    end
end

% connect element to new node
if ~strcmp(type,'foundation')
    element.(node_end)(ele_or_node_id) = new_node_id;
end

% Assign hinge properties
hinge.id(hinge_id,1) = hinge_id;
hinge.type{hinge_id,1} = type;
hinge.direction{hinge_id,1} = direction;
hinge.node_1(hinge_id,1) = new_node_id;
hinge.node_2(hinge_id,1) = old_node_id;
if strcmp(type,'foundation')
    hinge.element_id(hinge_id,1) = 0;
else
    hinge.element_id(hinge_id,1) = element.id(ele_or_node_id);
end
end

