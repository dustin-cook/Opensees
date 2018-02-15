function [ node, element, story, joint ] = fn_model_table( model )
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
            element.orientation(ele_id) = grid_line.orientation(e);
            element.story(ele_id) = s;
            element.depth(ele_id) = ele.d;
            element.width(ele_id) = ele.w;
            
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
                element.node_start(ele_id) = node_id;
            else % Existing Node
                element.node_start(ele_id) = node.id(start_node_check);
            end

            % Check to see if the ending node exists %COULD MAKE A FUNCTION
            end_node_check = (node.x == ele_x_end & node.y == ele_y_end & node.z == ele_z_end);
            if sum(end_node_check) == 0 % New Node
                node_id = node_id + 1;
                node.id(node_id) = node_id;
                node.x(node_id) = ele_x_end;
                node.y(node_id) = ele_y_end;
                node.z(node_id) = ele_z_end;
                element.node_end(ele_id) = node_id;
            else % Existing Node
                element.node_end(ele_id) = node.id(end_node_check);
            end
        end
    end
end

%% Assign Joints
joint = [];
joint_id = 0;
for s = 1:length(story.id)
    story_node = node.id(node.y == (story.y_offset(s)+story.story_ht(s)));
    for n = 1:length(story_node)
        n_id = story_node(n);
        elements_at_node = element.id(((element.node_start == n_id) | (element.node_end == n_id)) & (element.story == s));
        if length(elements_at_node) > 1
            for e = 1:length(elements_at_node)
                e_id = elements_at_node(e);
                
                % Find Joint properties based on elements that frame in
                bm_d_x = 0;
                bm_d_z = 0;
                if element.orientation(e_id) == 1 
                    new_ele.col_y = e_id;
                    joint_dim_x = element.depth(e_id);
                    joint_dim_z = element.width(e_id);
                elseif element.orientation(e_id) == 2
                    bm_d_x(e) = element.depth(e_id);
                    if element.node_start(e_id) == n_id
                        new_ele.bm_x_pos = e_id;
                    else
                        new_ele.bm_x_neg = e_id;
                    end
                elseif element.orientation(e_id) == 3
                    bm_d_z(e) = element.depth(e_id);
                    if element.node_start(e_id) == n_id
                        new_ele.bm_z_pos = e_id;
                    else
                        new_ele.bm_z_neg = e_id;
                    end
                elseif element.orientation(e_id) == 4
                    new_ele.col_y = e_id;
                    joint_dim_x = element.width(e_id);
                    joint_dim_z = element.depth(e_id);
                end
            end
            joint_dim_y = max([bm_d_x,bm_d_z]);
            
            % Define New Nodes
            new_node.x = [node.x(n_id)-joint_dim_x/2,node.x(n_id)+joint_dim_x/2,node.x(n_id),node.x(n_id),node.x(n_id),node.x(n_id)];
            new_node.y = [node.y(n_id)-joint_dim_y/2,node.y(n_id)-joint_dim_y/2,node.y(n_id)-joint_dim_y/2,node.y(n_id)-joint_dim_y/2,node.y(n_id)-joint_dim_y,node.y(n_id)-joint_dim_y/2];
            new_node.z = [node.z(n_id),node.z(n_id),node.z(n_id)-joint_dim_z/2,node.z(n_id)+joint_dim_z/2,node.z(n_id),node.z(n_id)];
            new_node.id = (1:6) + node.id(end);
            
            % Change elements to connect to new nodes
            if isfield(new_ele,'bm_x_neg')
                element.node_end(new_ele.bm_x_neg) = new_node.id(1);
            end
            if isfield(new_ele,'bm_x_pos')
                element.node_start(new_ele.bm_x_pos) = new_node.id(2);
            end
            if isfield(new_ele,'bm_z_neg')
                element.node_end(new_ele.bm_z_neg) = new_node.id(3);
            end
            if isfield(new_ele,'bm_z_pos')
                element.node_start(new_ele.bm_z_pos) = new_node.id(4);
            end
            element.node_end(new_ele.col_y) = new_node.id(5);
            
            clear new_ele
           
            % Define Joint
            joint_id = joint_id + 1;
            joint.id(joint_id) = joint_id + element.id(end);
            joint.x_neg(joint_id) = new_node.id(1);
            joint.x_pos(joint_id) = new_node.id(2);
            joint.y_neg(joint_id) = new_node.id(5);
            joint.y_pos(joint_id) = n_id;
            joint.z_neg(joint_id) = new_node.id(3);
            joint.z_pos(joint_id) = new_node.id(4);
            joint.center(joint_id) = new_node.id(6);
            
            % Add new nodes to nodes list
            node.id = [node.id new_node.id];
            node.x = [node.x new_node.x];
            node.y = [node.y new_node.y];
            node.z = [node.z new_node.z];

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

end

