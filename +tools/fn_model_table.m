function [ node, element, story, joint, wall ] = fn_model_table( model )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


import tools.*

%% Load data tables
story_table = readtable(['inputs' filesep model.name{1} filesep 'story.csv'],'ReadVariableNames',true);
story_group_table = readtable(['inputs' filesep model.name{1} filesep 'story_group.csv'],'ReadVariableNames',true);
grid_line_table = readtable(['inputs' filesep model.name{1} filesep 'grid_line.csv'],'ReadVariableNames',true);
element_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

node.id = 1;
node.x = 0;
node.y = 0;
node.z = 0;

% Object Methods
story = story_table(story_table.model_id == model.id,:);
element = [];
ele_id = 0;
for s = 1:length(story.id)
    story_group = story_group_table(story_group_table.story_group_id == story.story_group_id(s),:);
    for g = 1:length(story_group.id)
        frame_line = grid_line_table((grid_line_table.grid_line_id == story_group.grid_line_id(g)) & ~strcmp(element_table.type(grid_line_table.element_id),'wall'),:);
        for e = 1:length(frame_line.id)
            
            % Element Properties
            ele_id = ele_id + 1;
            ele = element_table(element_table.id == frame_line.element_id(e),:);
            element.id(ele_id) = ele_id;
            element.a(ele_id) = ele.a;
            element.e(ele_id) = ele.e;
            element.g(ele_id) = ele.g;
            element.j(ele_id) = ele.j;
            element.iy(ele_id) = ele.iy;
            element.iz(ele_id) = ele.iz;
            element.orientation(ele_id) = frame_line.orientation(e);
            element.story(ele_id) = s;
            element.depth(ele_id) = ele.d;
            element.width(ele_id) = ele.w;
            
            % Element Global Position
            ele_y_start = frame_line.y_start(e)*story.story_ht(s) + story.y_offset(s);
            ele_y_end = frame_line.y_end(e)*story.story_ht(s) + story.y_offset(s);
            
            if story_group.orientation(g) == 1
                ele_x_start = frame_line.x_start(e) + story_group.x_start(g);
                ele_x_end = frame_line.x_end(e) + story_group.x_start(g);
                ele_z_start = story_group.z_start(g);
                ele_z_end = story_group.z_start(g);
            elseif story_group.orientation(g) == 3
                ele_x_start = story_group.x_start(g);
                ele_x_end = story_group.x_start(g);
                ele_z_start = frame_line.x_start(e) + story_group.z_start(g);
                ele_z_end = frame_line.x_end(e) + story_group.z_start(g);
            else
                error('Grid Line Oreintation Not Valid')
            end

            % Check to see if the element nodes exists and assign
            [ node, id ] = node_exist( node, ele_x_start, ele_y_start, ele_z_start );
            element.node_start(ele_id) = node.id(id);
            [ node, id ] = node_exist( node, ele_x_end, ele_y_end, ele_z_end );
            element.node_end(ele_id) = node.id(id);
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
            
            bm_d_x = 0;
            bm_d_z = 0;
            col_d_x = 0;
            col_d_z = 0;
            for e = 1:length(elements_at_node)
                e_id = elements_at_node(e);
                % Find Joint properties based on elements that frame in
                if element.orientation(e_id) == 1 
                    new_ele.col_y = e_id;
                    col_d_x(e) = element.depth(e_id);
                    col_d_z(e) = element.width(e_id);
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
                    col_d_x(e) = element.width(e_id);
                    col_d_z(e) = element.depth(e_id);
                end
            end
            joint_dim_y = max([bm_d_x,bm_d_z]);
            joint_dim_x = max(col_d_x);
            joint_dim_z = max(col_d_z);
            
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

%% Assign Walls
wall = [];
wall_id = 0;
for s = 1:length(story.id)
    story_group = story_group_table(story_group_table.story_group_id == story.story_group_id(s),:);
    for g = 1:length(story_group.id)
        wall_line = grid_line_table((grid_line_table.grid_line_id == story_group.grid_line_id(g)) & strcmp(element_table.type(grid_line_table.element_id),'wall'),:);
        for w = 1:length(wall_line.id)
            
            % Element Properties
            wall_id = wall_id + 1;
            ele = element_table(element_table.id == wall_line.element_id(w),:);
            wall.id(wall_id) = wall_id;
            wall.ele_id(wall_id) = ele.id;
            wall.e(wall_id) = ele.e;
            wall.poisson_ratio(wall_id) = ele.poisson_ratio;
            wall.thickness(wall_id) = ele.w;
            
            % Element Global Position
            wall_y_start = wall_line.y_start(w)*story.story_ht(s) + story.y_offset(s);
            wall_y_end = wall_line.y_end(w)*story.story_ht(s) + story.y_offset(s);
            
            if story_group.orientation(g) == 1
                wall_x_start = wall_line.x_start(w) + story_group.x_start(g);
                wall_x_end = wall_line.x_end(w) + story_group.x_start(g);
                wall_z_start = story_group.z_start(g);
                wall_z_end = story_group.z_start(g);
            elseif story_group.orientation(g) == 3
                wall_x_start = story_group.x_start(g);
                wall_x_end = story_group.x_start(g);
                wall_z_start = wall_line.x_start(w) + story_group.z_start(g);
                wall_z_end = wall_line.x_end(w) + story_group.z_start(g);
            else
                error('Grid Line Oreintation Not Valid')
            end

            % Check to see if the wall nodes exist and assign nodes
            [ node, id ] = node_exist( node, wall_x_start, wall_y_start, wall_z_start );
            wall.node_1(wall_id) = node.id(id); 
            [ node, id ] = node_exist( node, wall_x_end, wall_y_start, wall_z_end );
            wall.node_2(wall_id) = node.id(id); 
            [ node, id ] = node_exist( node, wall_x_end, wall_y_end, wall_z_end );
            wall.node_3(wall_id) = node.id(id); 
            [ node, id ] = node_exist( node, wall_x_start, wall_y_end, wall_z_start );
            wall.node_4(wall_id) = node.id(id); 
            
        end
    end
end

%% Assign Nodal Fixity
node.fix = zeros(length(node.id),6);
foundation_nodes = (node.y == 0);
if strcmp(model.foundation,'fix')
    node.fix(foundation_nodes,:) = 1;
elseif strcmp(model.foundation,'pin')
    node.fix(foundation_nodes,:) = [1 1 1 0 0 0];
end

%% Find first nodes in each story
for s = 1:length(story.id)
    filter = (node.x == 0 & node.y == (story.y_offset(s)+story.story_ht(s)) & node.z == 0);
    story.first_story_node(s) = node.id(filter);
end

%% Assign weight to floor levels
node.weight = zeros(1,length(node.id));

for i = 1:length(story.id)
    slab_ht = story.y_offset(i) + story.story_ht(i);
    node_filter = (node.y == slab_ht);% & (node.z ~= 450) & (node.x == 0);
    story.nodes_on_slab{i} = node.id(node_filter);
    
    % Find slab extreme values
    max_x = max(node.x(node_filter));
    min_x = min(node.x(node_filter));
    max_z = max(node.z(node_filter));
    min_z = min(node.z(node_filter));
    
    % assign each node trib area ratio
    for j = 1:length(story.nodes_on_slab{i})
        n_id = story.nodes_on_slab{i}(j);
        if (node.x(n_id) == max_x || node.x(n_id) == min_x) && (node.z(n_id) == max_z || node.z(n_id) == min_z) % Corner Slab Node
            node.trib_area_ratio(n_id) = 0.25;
        elseif (node.x(n_id) == max_x || node.x(n_id) == min_x) || (node.z(n_id) == max_z || node.z(n_id) == min_z) % Side Slab Node
            node.trib_area_ratio(n_id) = 0.5;
        else % Center Slab Node
            node.trib_area_ratio(n_id) = 1;
        end
    end
    
    % total trib ratio
    wt_per_node = story.story_wt(i)/sum(node.trib_area_ratio(node_filter));
    node.weight(node_filter) = node.trib_area_ratio(node_filter) * wt_per_node;

end

%% Define Mass
node.mass = node.weight/386;

%% Offset Mass for Accidental Torsion
% For X direction %%%% SHOULD CHANGE INTO FUNCTION AND RUN EACH DIR
% INDEPENDANTLY
for i = 1:length(story.id)
    % Find Cental Coordinate and nodes to each side of it
    building_length = max(node.x(story.nodes_on_slab{i})) - min(node.x(story.nodes_on_slab{i}));
    node_center = building_length/2;
    nodes_right = story.nodes_on_slab{i}(node.x(story.nodes_on_slab{i}) > node_center);
    nodes_left = story.nodes_on_slab{i}(node.x(story.nodes_on_slab{i}) < node_center);
    mass_as_is = sum(node.mass(story.nodes_on_slab{i}));

    % Iteratevely Increase the Mass till a target offset is reached
    offset = 0;
    offset_target = 0.05; % as a percentage of the building length
    threshold = 0.001;
    mass_delta = 0;
    while abs(offset_target-offset) > threshold
        mass_delta = mass_delta + 0.001;
        
        %Infinite loop happen safety device
        if mass_delta > 0.5
            error('Could Not FInd Accidental Torsion Adjustment, Recalibrate')
        end
        
        temp_mass = node.mass;
        temp_mass(nodes_right) = node.mass(nodes_right)*(1+mass_delta);
        temp_mass(nodes_left) = node.mass(nodes_left)*(1-mass_delta);
        
        % Make sure mass still adds up to what it did before
        mass_new = sum(temp_mass(story.nodes_on_slab{i}));
        mass_diff = round(abs(mass_new - mass_as_is),2);
        if mass_diff > 0
            error('Mass got thrown off trying to find accidental torsion')
        end
        
        % Find new offset
        offset = (sum(temp_mass(story.nodes_on_slab{i}).*(node.x(story.nodes_on_slab{i}) - node_center))/sum(temp_mass(story.nodes_on_slab{i})))/building_length;
    end
    
    % Set Found Masses
    node.mass(nodes_right) = node.mass(nodes_right)*(1+mass_delta);
    node.mass(nodes_left) = node.mass(nodes_left)*(1-mass_delta);
end

% For Z direction
for i = 1:length(story.id)
    % Find Cental Coordinate and nodes to each side of it
    building_length = max(node.z(story.nodes_on_slab{i})) - min(node.z(story.nodes_on_slab{i}));
    node_center = building_length/2;
    nodes_right = story.nodes_on_slab{i}(node.z(story.nodes_on_slab{i}) > node_center);
    nodes_left = story.nodes_on_slab{i}(node.z(story.nodes_on_slab{i}) < node_center);
    mass_as_is = sum(node.mass(story.nodes_on_slab{i}));

    % Iteratevely Increase the Mass till a target offset is reached
    offset = 0;
    offset_target = 0.05; % as a percentage of the building length
    threshold = 0.001;
    mass_delta = 0;
    while abs(offset_target-offset) > threshold
        mass_delta = mass_delta + 0.01;
        
        %Infinite loop happen safety device
        if mass_delta > 0.5
            error('Could Not FInd Accidental Torsion Adjustment, Recalibrate')
        end
        
        temp_mass = node.mass;
        temp_mass(nodes_right) = node.mass(nodes_right)*(1+mass_delta);
        temp_mass(nodes_left) = node.mass(nodes_left)*(1-mass_delta);
        
        % Make sure mass still adds up to what it did before
        mass_new = sum(temp_mass(story.nodes_on_slab{i}));
        mass_diff = round(abs(mass_new - mass_as_is),2);
        if mass_diff > 0
            error('Mass got thrown off trying to find accidental torsion')
        end
        
        % Find new offset
        offset = (sum(temp_mass(story.nodes_on_slab{i}).*(node.z(story.nodes_on_slab{i}) - node_center))/sum(temp_mass(story.nodes_on_slab{i})))/building_length;
    end
    
    % Set Found Masses
    node.mass(nodes_right) = node.mass(nodes_right)*(1+mass_delta);
    node.mass(nodes_left) = node.mass(nodes_left)*(1-mass_delta);
end

end

