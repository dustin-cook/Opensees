function [ node_table, ele_table, story, joint, hinge ] = fn_model_table( model, analysis )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% INITIAL SETUP
% Import Packages
import tools.*

% Load data tables
story_table = readtable(['inputs' filesep 'models' filesep model.name{1} filesep 'story.csv'],'ReadVariableNames',true);
story_group_table = readtable(['inputs' filesep 'models' filesep model.name{1} filesep 'story_group.csv'],'ReadVariableNames',true);
grid_line_table = readtable(['inputs' filesep 'models' filesep model.name{1} filesep 'grid_line.csv'],'ReadVariableNames',true);
ele_props = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

node.id = 1;
node.x = 0;
node.y = 0;
node.z = 0;

%% Assign Elements
story = story_table(story_table.model_id == model.id,:);
element = [];
ele_id = 0;
for s = 1:length(story.id)
    story_group = story_group_table(story_group_table.story_group_id == story.story_group_id(s),:);
    for g = 1:length(story_group.id)
        frame_line = grid_line_table((grid_line_table.grid_line_id == story_group.grid_line_id(g)) & ~strcmp(ele_props.type(grid_line_table.element_id),'wall'),:);
        for e = 1:length(frame_line.id)
            
            % Element Properties
            ele_id = ele_id + 1;
            ele = ele_props(ele_props.id == frame_line.element_id(e),:);
            element.id(ele_id,1) = ele_id;
            element.ele_id(ele_id,1) = ele.id;
            element.orientation(ele_id,1) = frame_line.orientation(e);
            element.story(ele_id,1) = s;
            element.type(ele_id,1) = ele.type;
            
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
            element.node_1(ele_id,1) = node.id(id);
            [ node, id ] = node_exist( node, ele_x_end, ele_y_end, ele_z_end );
            element.node_2(ele_id,1) = node.id(id);
            element.node_3(ele_id,1) = 0;
            element.node_4(ele_id,1) = 0;
            element.length(ele_id,1) = sqrt( (ele_x_end-ele_x_start)^2 + (ele_y_end-ele_y_start)^2 + (ele_z_end-ele_z_start)^2 );
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
        elements_at_node = element.id(((element.node_1 == n_id) | (element.node_2 == n_id)) & (element.story == s));
        if length(elements_at_node) > 1
            
            bm_d_x = 0;
            bm_d_z = 0;
            col_d_x = 0;
            col_d_z = 0;
            for e = 1:length(elements_at_node)
                e_id = elements_at_node(e);
                ele = ele_props(ele_props.id == element.ele_id(e_id),:);
                
                % Find Joint properties based on elements that frame in
                if element.orientation(e_id) == 1 
                    new_ele.col_y = e_id;
                    col_d_x(e) = ele.d;
                    col_d_z(e) = ele.w;
                elseif element.orientation(e_id,1) == 2
                    bm_d_x(e) = ele.d;
                    if element.node_1(e_id) == n_id
                        new_ele.bm_x_pos = e_id;
                    else
                        new_ele.bm_x_neg = e_id;
                    end
                elseif element.orientation(e_id) == 3
                    bm_d_z(e) = ele.d;
                    if element.node_1(e_id) == n_id
                        new_ele.bm_z_pos = e_id;
                    else
                        new_ele.bm_z_neg = e_id;
                    end
                elseif element.orientation(e_id) == 4
                    new_ele.col_y = e_id;
                    col_d_x(e) = ele.w;
                    col_d_z(e) = ele.d;
                end
            end
            joint_dim_y = max([bm_d_x,bm_d_z]);
            joint_dim_x = max(col_d_x);
%             joint_dim_z = max(col_d_z);
            
            % Define New Nodes
            new_node.x = [node.x(n_id)-joint_dim_x/2,node.x(n_id)+joint_dim_x/2,node.x(n_id),node.x(n_id)];
            new_node.y = [node.y(n_id)-joint_dim_y/2,node.y(n_id)-joint_dim_y/2,node.y(n_id)-joint_dim_y,node.y(n_id)-joint_dim_y/2];
%             new_node.z = [node.z(n_id),node.z(n_id),node.z(n_id)-joint_dim_z/2,node.z(n_id)+joint_dim_z/2,node.z(n_id),node.z(n_id)];
            new_node.z = [node.z(n_id),node.z(n_id),node.z(n_id),node.z(n_id)];
            new_node.id = (1:4) + node.id(end);
            
            % Change elements to connect to new nodes
            if isfield(new_ele,'bm_x_neg')
                element.node_2(new_ele.bm_x_neg,1) = new_node.id(1);
            end
            if isfield(new_ele,'bm_x_pos')
                element.node_1(new_ele.bm_x_pos,1) = new_node.id(2);
            end
%             if isfield(new_ele,'bm_z_neg')
%                 element.node_2(new_ele.bm_z_neg,1) = new_node.id(3);
%             end
%             if isfield(new_ele,'bm_z_pos')
%                 element.node_1(new_ele.bm_z_pos,1) = new_node.id(4);
%             end
            element.node_2(new_ele.col_y,1) = new_node.id(3);
            
            clear new_ele
           
            % Define Joint
            joint_id = joint_id + 1;
            joint.id(joint_id) = joint_id + 1000;
            joint.x_neg(joint_id) = new_node.id(1);
            joint.x_pos(joint_id) = new_node.id(2);
            joint.y_neg(joint_id) = new_node.id(3);
            joint.y_pos(joint_id) = n_id;
%             joint.z_neg(joint_id) = new_node.id(3);
%             joint.z_pos(joint_id) = new_node.id(4);
            joint.center(joint_id) = new_node.id(4);
            
            % Add new nodes to nodes list
            node.id = [node.id new_node.id];
            node.x = [node.x new_node.x];
            node.y = [node.y new_node.y];
            node.z = [node.z new_node.z];

        end
    end
end

%% Assign Walls
for s = 1:length(story.id)
    story_group = story_group_table(story_group_table.story_group_id == story.story_group_id(s),:);
    for g = 1:length(story_group.id)
        wall_line = grid_line_table((grid_line_table.grid_line_id == story_group.grid_line_id(g)) & strcmp(ele_props.type(grid_line_table.element_id),'wall'),:);
        for w = 1:length(wall_line.id)
            
            % Element Properties
            ele_id = ele_id + 1;
            ele = ele_props(ele_props.id == wall_line.element_id(w),:);
            element.id(ele_id,1) = ele_id;
            element.ele_id(ele_id,1) = ele.id;
            element.orientation(ele_id,1) = 0;
            element.story(ele_id,1) = s;
            element.type{ele_id,1} = 'wall';
            
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
            element.node_1(ele_id,1) = node.id(id); 
            [ node, id ] = node_exist( node, wall_x_end, wall_y_start, wall_z_end );
            element.node_2(ele_id,1) = node.id(id); 
            [ node, id ] = node_exist( node, wall_x_end, wall_y_end, wall_z_end );
            element.node_3(ele_id,1) = node.id(id); 
            [ node, id ] = node_exist( node, wall_x_start, wall_y_end, wall_z_start );
            element.node_4(ele_id,1) = node.id(id); 
            element.length(ele_id,1) = sqrt( (wall_x_end-wall_x_start)^2 + (wall_y_end-wall_y_start)^2 + (wall_z_start-wall_z_start)^2 );
        end
    end
end

%% Find first nodes in each story
for s = 1:length(story.id)
    filter = (node.x == 0 & node.y == (story.y_offset(s)+story.story_ht(s)) & node.z == 0);
    story.first_story_node(s) = node.id(filter);
end

%% Assign weight to floor levels
node.dead_load = zeros(1,length(node.id));
node.live_load = zeros(1,length(node.id));
node.trib_area_ratio = zeros(1,length(node.id));

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
    node.dead_load(node_filter) = node.trib_area_ratio(node_filter) * story.story_dead_load(i) / sum(node.trib_area_ratio(node_filter));
    node.live_load(node_filter) = node.trib_area_ratio(node_filter) * story.story_live_load(i) / sum(node.trib_area_ratio(node_filter));

end

%% Define Mass
node.mass = node.dead_load/386;

%% Offset Mass for Accidental Torsion
% For X direction %%%% SHOULD CHANGE INTO FUNCTION AND RUN EACH DIR
% INDEPENDANTLY
if analysis.accidental_torsion == 1
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

%% Define Nodal Fixity
node.fix = cell(length(node.id),1);
node.fix(1:end) = {[0 0 0 0 0 0]}; % All nodes
foundation_nodes_id = node.id(node.y == 0);
node.fix(foundation_nodes_id,:) = {[1 1 1 1 1 1]}; % foundation nodes

%% Create Nonlinear Rotational Springs at ends of all beams and columns
hinge = [];
hinge_id = 0;
if analysis.nonlinear ~= 0
    % Define Hinges
    for i = 1:length(element.id)
        if strcmp(element.type{i},'column') || strcmp(element.type{i},'beam') % For all columns and beams
            hinge_id = hinge_id+1;
            % Define hinge at start of element
            [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_1', i, hinge_id, foundation_nodes_id ); 
            hinge_id = hinge_id+1;
            % Define hinge at end of element
            [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_2', i, hinge_id, foundation_nodes_id );
        end
    end
end


%% Remove any Z dimension Nodes for 2D Analysis
clear new_node
if strcmp(model.dimension,'2D')
    non_z_nodes = (node.z == 0);
    new_node.id = node.id(non_z_nodes);
    new_node.x = node.x(non_z_nodes);
    new_node.y = node.y(non_z_nodes);
    new_node.dead_load = node.dead_load(non_z_nodes);
    new_node.live_load = node.live_load(non_z_nodes);
    new_node.mass = node.mass(non_z_nodes);
    new_node.fix = node.fix(non_z_nodes,:);
    node = new_node;
end

%% Reformat outputs to table Outputs Tables
if isfield(node,'z')
    node_table = table(node.id', node.x', node.y', node.z', node.dead_load', node.live_load', node.trib_area_ratio', node.mass', node.fix, 'VariableNames',{'id','x','y','z','dead_load','live_load','trib_area_ratio','mass','fix'});
else
    node_table = table(node.id', node.x', node.y', node.dead_load', node.live_load', node.mass', node.fix, 'VariableNames',{'id','x','y','dead_load','live_load','mass','fix'});
end
ele_table = struct2table(element);

end

