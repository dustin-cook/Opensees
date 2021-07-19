function [ node, element ] = fn_create_element( ele_type, ele_id, ele_props, element_id, nb, story_props, story_group, node, element, direction, trib_wt_1, trib_wt_2 )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

import build_model.fn_node_exist

element.id(ele_id,1) = ele_id;

% Element Global Position
start_bay = (nb-1) + find(story_props.(['bay_coor_' direction]){1} == story_group.([direction '_start']));
ele_prim_start = story_props.(['bay_coor_' direction]){1}(start_bay);
ele_y_start = story_props.y_start;
ele_y_end = story_props.y_start + story_props.story_ht;

% Assign Tributary weights
element.trib_wt_1(ele_id,1) = trib_wt_1;
element.trib_wt_2(ele_id,1) = trib_wt_2;

% Assign element type specific properties properties
ele = ele_props(ele_props.id == element_id,:);
if strcmp(ele_type,'col')
    ele_prim_end = ele_prim_start;
elseif strcmp(ele_type,'beam')
    ele_prim_end = story_props.(['bay_coor_' direction]){1}(start_bay + 1);
    ele_y_start = ele_y_end;
elseif strcmp(ele_type,'wall')
    ele_prim_end = ele_prim_start;
%     ele_prim_end = story_group.([direction '_start']) + story_props.(['bay_coor_' direction]){1}(start_bay + 1);
end

% General Element Properties
element.ele_id(ele_id,1) = ele.id;
element.direction{ele_id,1} = direction;
element.story(ele_id,1) = story_props.id;
element.type(ele_id,1) = ele.type;

% Define element locations
if strcmp(direction,'x')
    ele_x_start = ele_prim_start;
    ele_x_end = ele_prim_end;
    ele_z_start = story_group.z_start;
    ele_z_end = story_group.z_start;
elseif strcmp(direction,'z')
    ele_z_start = ele_prim_start;
    ele_z_end = ele_prim_end;
    ele_x_start = story_group.x_start;
    ele_x_end = story_group.x_start;
else
    error('Invalid Element Group Direction')
end
    
% Assign nodes to elements (create nodes if they do not already exist)
% if strcmp(ele_type,'wall')
%     [ node, id ] = fn_node_exist( node, ele_x_start, ele_y_start, ele_z_start );
%     element.node_1(ele_id,1) = node.id(id);
%     [ node, id ] = fn_node_exist( node, ele_x_end, ele_y_start, ele_z_end );
%     element.node_2(ele_id,1) = node.id(id);
%     [ node, id ] = fn_node_exist( node, ele_x_end, ele_y_end, ele_z_end );
%     element.node_3(ele_id,1) = node.id(id);
%     [ node, id ] = fn_node_exist( node, ele_x_start, ele_y_end, ele_z_start );
%     element.node_4(ele_id,1) = node.id(id);
% else
    [ node, id ] = fn_node_exist( node, ele_x_start, ele_y_start, ele_z_start );
    element.node_1(ele_id,1) = node.id(id);
    [ node, id ] = fn_node_exist( node, ele_x_end, ele_y_end, ele_z_end );
    element.node_2(ele_id,1) = node.id(id);
    element.node_3(ele_id,1) = 0;
    element.node_4(ele_id,1) = 0;
% end

end

