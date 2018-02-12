function [ node, element, story ] = fn_model_table( model )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Load data tables
story_table = readtable(['inputs' filesep model.name{1} filesep 'story.csv'],'ReadVariableNames',true);
story_group_table = readtable(['inputs' filesep model.name{1} filesep 'story_group.csv'],'ReadVariableNames',true);
grid_line_table = readtable(['inputs' filesep model.name{1} filesep 'grid_line.csv'],'ReadVariableNames',true);
element_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

node.id = 1;
node.x = 0;
node.y = 0;
node.z = 0;
node.weight = 0;

% Object Methods
story = story_table(story_table.model_id == model.id,:);
ele_id = 0;
node_id = 1;
for s = 1:length(story.id)
    story_group = story_group_table(story_group_table.story_group_id == story.story_group_id(s),:);
    for g = 1:length(story_group.id)
        grid_line = grid_line_table(grid_line_table.grid_line_id == story_group.grid_line_id(g),:);
        for e = 1:length(grid_line.id)
            
            % Element Properties
            ele_id = ele_id + 1;
            ele = element_table(element_table.id == grid_line.element_id(e),:);
            element.id(ele_id) = ele_id;
            element.a(ele_id) = ele.a;
            element.e(ele_id) = ele.e;
            element.g(ele_id) = ele.g;
            element.j(ele_id) = ele.j;
            element.iy(ele_id) = ele.iy;
            element.iz(ele_id) = ele.iz;
            element.weight(ele_id) = grid_line.weight(e);
            element.orientation(ele_id) = grid_line.orientation(e);
            element.story(ele_id) = s;
            element.depth(ele_id) = ele.d;
            
            % Element Global Position
            ele_y_start = grid_line.y_start(e)*story.story_ht(s) + story.y_offset(s);
            ele_y_end = grid_line.y_end(e)*story.story_ht(s) + story.y_offset(s);
            
            if story_group.orientation(g) == 1
                ele_x_start = grid_line.x_start(e) + story_group.x_start(g);
                ele_x_end = grid_line.x_end(e) + story_group.x_start(g);
                ele_z_start = story_group.z_start(g);
                ele_z_end = story_group.z_start(g);
            elseif story_group.orientation(g) == 3
                ele_x_start = story_group.x_start(g);
                ele_x_end = story_group.x_start(g);
                ele_z_start = grid_line.x_start(e) + story_group.z_start(g);
                ele_z_end = grid_line.x_end(e) + story_group.z_start(g);
            else
                error('Grid Line Oreintation Not Valid')
            end


            % Check to see if the starting node exists %COULD MAKE A FUNCTION
            start_node_check = (node.x == ele_x_start & node.y == ele_y_start & node.z == ele_z_start);
            if sum(start_node_check) == 0 % New Node
                node_id = node_id + 1;
                node.id(node_id) = node_id;
                node.x(node_id) = ele_x_start;
                node.y(node_id) = ele_y_start;
                node.z(node_id) = ele_z_start;
                node.weight(node_id) = grid_line.weight(e)/2;
                element.node_start(ele_id) = node_id;
            else % Existing Node
                element.node_start(ele_id) = node.id(start_node_check);
                node.weight(start_node_check) = node.weight(start_node_check) + grid_line.weight(e)/2;
            end

            % Check to see if the ending node exists %COULD MAKE A FUNCTION
            end_node_check = (node.x == ele_x_end & node.y == ele_y_end & node.z == ele_z_end);
            if sum(end_node_check) == 0 % New Node
                node_id = node_id + 1;
                node.id(node_id) = node_id;
                node.x(node_id) = ele_x_end;
                node.y(node_id) = ele_y_end;
                node.z(node_id) = ele_z_end;
                node.weight(node_id) = grid_line.weight(e)/2;
                element.node_end(ele_id) = node_id;
            else % Existing Node
                element.node_end(ele_id) = node.id(end_node_check);
                node.weight(end_node_check) = node.weight(end_node_check) + grid_line.weight(e)/2;
            end
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
node.mass = node.weight/386;
node.force = zeros(1,length(node.id));


% Assign Joints
for s = 1:length(story.id)
    story_node = node.id(node.y == (story.y_offset(s)+story.story_ht(s)));
    for n = 1:length(story_node)
        elements_at_node = element.id(((element.node_start == story_node(n)) | (element.node_end == story_node(n))) & (element.story == s));
        if length(elements_at_node) > 1
            for e = 1:length(elements_at_node)
                e_id = elements_at_node(e);
                if element.orientation(e_id) == 1 
                    col_d_x = element.depth(e_id);
                elseif element.orientation(e_id) == 2
                    bm_d_x = element.depth(e_id);
                elseif element.orientation(e_id) == 3
                    bm_d_z = element.depth(e_id);
                elseif element.orientation(e_id) == 4
                    col_d_z = element.depth(e_id);
                end
            end
            
            for i = 1:6
                node_id = node_id + 1;
                node.id(node_id) = node_id;
                node.x(node_id) = ele_x_end;
                node.y(node_id) = ele_y_end;
                node.z(node_id) = ele_z_end;
                node.weight(node_id) = grid_line.weight(e)/2;
                ele_id = ele_id + 1;
            end
        end
    end
end

end

