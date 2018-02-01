function [ node, element, story ] = fn_model_table( model )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Load data tables
story_table = readtable(['inputs' filesep 'story.csv'],'ReadVariableNames',true);
element_group_table = readtable(['inputs' filesep 'element_group.csv'],'ReadVariableNames',true);
element_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
% node_group_table = readtable(['inputs' filesep 'node_group.csv'],'ReadVariableNames',true);

node.id = 1;
node.x = 0;
node.y = 0;
node.z = 0;
node.weight = 0;

% Object Methods
story = story_table(story_table.model == model.id,:);
ele_id = 0;
node_id = 1;
for s = 1:length(story.id)
    element_group = element_group_table(element_group_table.element_group == story.element_group(s),:);
    element_group.y_start = element_group.y_start + story.y_offset(s);
    element_group.y_end = element_group.y_end + story.y_offset(s);
    for e = 1:length(element_group.id)
        ele_id = ele_id + 1;
        ele = element_table(element_table.id == element_group.element(e),:);
        element.id(ele_id) = ele_id;
        element.a(ele_id) = ele.a;
        element.e(ele_id) = ele.e;
        element.g(ele_id) = ele.g;
        element.j(ele_id) = ele.j;
        element.iy(ele_id) = ele.iy;
        element.iz(ele_id) = ele.iz;
        element.weight(ele_id) = element_group.weight(e);
        element.orientation(ele_id) = element_group.orientation(e);
        
        
        % Check to see if the starting node exists %COULD MAKE A FUNCTION
        start_node_check = (node.x == element_group.x_start(e) & node.y == element_group.y_start(e) & node.z == element_group.z_start(e));
        if sum(start_node_check) == 0 % New Node
            node_id = node_id + 1;
            node.id(node_id) = node_id;
            node.x(node_id) = element_group.x_start(e);
            node.y(node_id) = element_group.y_start(e);
            node.z(node_id) = element_group.z_start(e);
            node.weight(node_id) = element_group.weight(e)/2;
            element.node_start(ele_id) = node_id;
        else % Existing Node
            element.node_start(ele_id) = node.id(start_node_check);
            node.weight(start_node_check) = node.weight(start_node_check) + element_group.weight(e)/2;
        end

        % Check to see if the ending node exists %COULD MAKE A FUNCTION
        end_node_check = (node.x == element_group.x_end(e) & node.y == element_group.y_end(e) & node.z == element_group.z_end(e));
        if sum(end_node_check) == 0 % New Node
            node_id = node_id + 1;
            node.id(node_id) = node_id;
            node.x(node_id) = element_group.x_end(e);
            node.y(node_id) = element_group.y_end(e);
            node.z(node_id) = element_group.z_end(e);
            node.weight(node_id) = element_group.weight(e)/2;
            element.node_end(ele_id) = node_id;
        else % Existing Node
            element.node_end(ele_id) = node.id(end_node_check);
            node.weight(end_node_check) = node.weight(end_node_check) + element_group.weight(e)/2;
        end
    end
end

% Assign Nodal Fixity
node.fix = zeros(length(node.id),6);
foundation_nodes = (node.y == 0);
if strcmp(model.foundation,'fix')
    node.fix(foundation_nodes,:) = 1;
elseif strcmp(model.foundation,'fix')
    node.fix(foundation_nodes,:) = [1 1 1 0 0 0];
end

% Assign Node Wt and Lateral Force
node.mass = node.weight*100/386;
node.force = zeros(1,length(node.id));




end

