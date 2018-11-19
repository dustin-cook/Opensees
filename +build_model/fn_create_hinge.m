function [ node, element, hinge ] = fn_create_hinge( node, element, hinge, node_end, ele_or_node_id, hinge_id, foundation_nodes_id, type )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Define new nodes to connect to springs
new_node_id = node.id(end) + 1;
if strcmp(type,'foundation')
    old_node_id = ele_or_node_id;
else
    old_node_id = element.(node_end)(ele_or_node_id);
end

% Define new node properties
node.id(new_node_id) = new_node_id;
node.x(new_node_id) = node.x(old_node_id);
node.y(new_node_id) = node.y(old_node_id);
node.z(new_node_id) = node.z(old_node_id);
node.dead_load(new_node_id) = 0;
node.live_load(new_node_id) = 0;
node.mass(new_node_id) = 0;
node.story(new_node_id) = 0;
node.primary_story(new_node_id) = 0;
node.fix(new_node_id) = node.fix(old_node_id);
node.on_slab(new_node_id) = 0;

% Define fixity of foundation nodes
if (sum(old_node_id == foundation_nodes_id) > 0)
    if strcmp(type,'foundation')
        node.fix{new_node_id} = '[111111]';
    else
        node.fix{new_node_id} = '[000000]';
    end
end

% connect element to new node
if ~strcmp(type,'foundation')
    element.(node_end)(ele_or_node_id) = new_node_id;
end

% Assign hinge properties
hinge.id(hinge_id,1) = hinge_id;
hinge.type{hinge_id,1} = type;
hinge.node_1(hinge_id,1) = new_node_id;
hinge.node_2(hinge_id,1) = old_node_id;
if strcmp(type,'foundation')
    hinge.element_id(hinge_id,1) = 0;
else
    hinge.element_id(hinge_id,1) = element.id(ele_or_node_id);
end
end

