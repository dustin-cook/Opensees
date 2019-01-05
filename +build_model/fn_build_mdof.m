function [ ] = fn_build_mdof( model, ele_props_table, analysis, write_dir, read_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% INITIAL SETUP
% Import Packages
import build_model.*
import aci_318.fn_aci_moment_capacity

% Load data tables
story = readtable([read_dir filesep 'story.csv'],'ReadVariableNames',true);
story_group_table = readtable([read_dir filesep 'story_group.csv'],'ReadVariableNames',true);
element_group_table = readtable([read_dir filesep 'element_group.csv'],'ReadVariableNames',true);
additional_elements = readtable([read_dir filesep 'additional_elements.csv'],'ReadVariableNames',true);

% Story property calculations
for s = 1:height(story)
    story.y_start(s) = sum([0;story.story_ht(1:s-1)]);
    story.bay_length_x{s} = str2double(strsplit(strrep(strrep(story.bay_length_x{s},']',''),'[',''),','));
    story.bay_length_z{s} = str2double(strsplit(strrep(strrep(story.bay_length_z{s},']',''),'[',''),','));
    for x = 1:(length(story.bay_length_x{s})+1)
        bay_coor_x(x) = sum([0,story.bay_length_x{s}(1:x-1)]);
    end
    story.bay_coor_x{s} = bay_coor_x;
    for z = 1:(length(story.bay_length_z{s})+1)
        bay_coor_z(z) = sum([0,story.bay_length_z{s}(1:z-1)]);
    end
    story.bay_coor_z{s} = bay_coor_z;
end

%% Assign Elements
node = [];
element = [];
ele_id = 0;
for s = 1:height(story)
    story_props = story(s,:);
    story_group = story_group_table(story_group_table.story_group_id == story.story_group_id(s),:);
    for g = 1:height(story_group)
        element_group = element_group_table((element_group_table.id == story_group.element_group_id(g)),:);
        for e = 1:height(element_group)
            for nb = 1:element_group.num_bays
                % Element Properties Column
                if element_group.col_id ~= 0
                    ele_id = ele_id + 1;
                    [ node, element ] = fn_create_element( 'col', ele_id, ele_props_table, element_group, nb, story_props, story_group(g,:), node, element, story_group.direction{g} );
                end
                if element_group.beam_id ~= 0
                    ele_id = ele_id + 1;
                    [ node, element ] = fn_create_element( 'beam', ele_id, ele_props_table, element_group, nb, story_props, story_group(g,:), node, element, story_group.direction{g} );
                end
                if element_group.wall_id ~= 0
                    ele_id = ele_id + 1;
                    [ node, element ] = fn_create_element( 'wall', ele_id, ele_props_table, element_group, nb, story_props, story_group(g,:), node, element, story_group.direction{g} );
                end
            end
            % Last column in bay span
            if element_group.col_id ~= 0
                ele_id = ele_id + 1;
                [ node, element ] = fn_create_element( 'col', ele_id, ele_props_table, element_group, element_group.num_bays+1, story_props, story_group(g,:), node, element, story_group.direction{g} );
            end
        end
    end
    
    % Assign Gravity Load to Elements
    sum_trib_wt = sum(element.trib_wt(element.story == story.id(s)));
    element.dead_load(element.story == story.id(s),1) = element.trib_wt(element.story == story.id(s))*story.story_dead_load(s)/sum_trib_wt;
    element.live_load(element.story == story.id(s),1) = element.trib_wt(element.story == story.id(s))*story.story_live_load(s)/sum_trib_wt;
end

%% For each node created go through as say whats connected in
mf_id = 0;
mf_joint.id = [];
for i = 1:length(node.id)
    eles_at_node = element.id(((element.node_1 == node.id(i)) | (element.node_2 == node.id(i))));
    if length(eles_at_node) > 1 && ~strcmp(element.type(element.id == eles_at_node(1)),'wall')
        mf_id = mf_id + 1;
        mf_joint.id(mf_id,1) = mf_id;
        mf_joint.column_low(mf_id,1) = 0;
        mf_joint.column_high(mf_id,1) = 0;
        mf_joint.beam_left(mf_id,1) = 0;
        mf_joint.beam_right(mf_id,1) = 0;
        for e = 1:length(eles_at_node)
            if strcmp(element.type(element.id == eles_at_node(e)),'column')
                if element.node_1(element.id == eles_at_node(e)) == node.id(i)
                    mf_joint.column_high(mf_id,1) = eles_at_node(e);
                elseif element.node_2(element.id == eles_at_node(e)) == node.id(i)
                    mf_joint.column_low(mf_id,1) = eles_at_node(e);
                end
            elseif strcmp(element.type(element.id == eles_at_node(e)),'beam')
                if element.node_1(element.id == eles_at_node(e)) == node.id(i)
                    mf_joint.beam_right(mf_id,1) = eles_at_node(e);
                elseif element.node_2(element.id == eles_at_node(e)) == node.id(i)
                    mf_joint.beam_left(mf_id,1) = eles_at_node(e);
                end
            end
        end
    end
end

%% Calculate Total Element Gravity Load
element.gravity_load = element.dead_load*analysis.dead_load + element.live_load*analysis.live_load;

%% Assign Gravity Loads to Nodes
node.dead_load = zeros(length(node.id),1);
node.live_load = zeros(length(node.id),1);
for e = 1:length(element.id)
    if strcmp(element.type{e},'wall')
        node.dead_load(element.node_2(e)) = node.dead_load(element.node_2(e)) + element.dead_load(e);
        node.live_load(element.node_2(e)) = node.live_load(element.node_2(e)) + element.live_load(e);
%         node.dead_load(element.node_3(e)) = node.dead_load(element.node_3(e)) + element.dead_load(e)/2;
%         node.dead_load(element.node_4(e)) = node.dead_load(element.node_4(e)) + element.dead_load(e)/2;
%         node.live_load(element.node_3(e)) = node.live_load(element.node_3(e)) + element.live_load(e)/2;
%         node.live_load(element.node_4(e)) = node.live_load(element.node_4(e)) + element.live_load(e)/2;
    else
        node.dead_load(element.node_1(e)) = node.dead_load(element.node_1(e)) + element.dead_load(e)/2;
        node.dead_load(element.node_2(e)) = node.dead_load(element.node_2(e)) + element.dead_load(e)/2;
        node.live_load(element.node_1(e)) = node.live_load(element.node_1(e)) + element.live_load(e)/2;
        node.live_load(element.node_2(e)) = node.live_load(element.node_2(e)) + element.live_load(e)/2;
    end
end

%% Define Mass
node.mass = node.dead_load/386;

%% Assign Joints
joint.id = [];
joint.x_neg = [];
joint.x_pos = [];
joint.y_neg = [];
joint.y_pos = [];
if strcmp(model.dimension,'3D')
    joint.z_neg = [];
    joint.z_pos = [];
end
joint_id = 0;
for s = 1:height(story)
    story_node = node.id(node.y == (story.y_start(s)+story.story_ht(s)));
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
                ele = ele_props_table(ele_props_table.id == element.ele_id(e_id),:);
                
                % Find Joint properties based on elements that frame in
                if strcmp(element.type{e_id},'column')
                    new_ele.col_y = e_id;
                    if strcmp(element.direction{e_id},'x')
                        col_d_x(e) = ele.d;
                        col_d_z(e) = ele.w;
                    elseif strcmp(element.direction{e_id},'z')
                        col_d_x(e) = ele.w;
                        col_d_z(e) = ele.d;
                    else
                        error('element direction not recognized')
                    end
                elseif strcmp(element.type{e_id},'beam')
                    if strcmp(element.direction{e_id},'x')
                        bm_d_x(e) = ele.d;
                        if element.node_1(e_id) == n_id
                            new_ele.bm_x_pos = e_id;
                        else
                            new_ele.bm_x_neg = e_id;
                        end
                    elseif strcmp(element.direction{e_id},'z')
                        bm_d_z(e) = ele.d;
                        if element.node_1(e_id) == n_id
                            new_ele.bm_z_pos = e_id;
                        else
                            new_ele.bm_z_neg = e_id;
                        end
                    else
                        error('element direction not recognized')
                    end
                end
            end
            joint_dim_y = max([bm_d_x,bm_d_z]);
            joint_dim_x = max(col_d_x);
            joint_dim_z = max(col_d_z);
            
            % Define New Nodes
            if strcmp(model.dimension,'2D')
                new_node.id = [1;2;3] + node.id(end);
                new_node.x = [node.x(n_id)-joint_dim_x/2;node.x(n_id)+joint_dim_x/2;node.x(n_id)];
                new_node.y = [node.y(n_id)-joint_dim_y/2;node.y(n_id)-joint_dim_y/2;node.y(n_id)-joint_dim_y];
                new_node.z = [node.z(n_id);node.z(n_id);node.z(n_id)];
                new_node.dead_load = [0; 0; 0];
                new_node.live_load = [0; 0; 0];
                new_node.mass = [0; 0; 0];
            elseif strcmp(model.dimension,'3D')
                new_node.id = [1;2;3;4;5] + node.id(end);
                new_node.x = [node.x(n_id)-joint_dim_x/2;node.x(n_id)+joint_dim_x/2;node.x(n_id);node.x(n_id);node.x(n_id)];
                new_node.y = [node.y(n_id)-joint_dim_y/2;node.y(n_id)-joint_dim_y/2;node.y(n_id)-joint_dim_y;node.y(n_id)-joint_dim_y/2;node.y(n_id)-joint_dim_y/2];
                new_node.z = [node.z(n_id);node.z(n_id);node.z(n_id);node.z(n_id)-joint_dim_z/2;node.z(n_id)+joint_dim_z/2]; 
                new_node.dead_load = [0; 0; 0; 0; 0];
                new_node.live_load = [0; 0; 0; 0; 0];
                new_node.mass = [0; 0; 0; 0; 0];
            end
            
            % Change elements to connect to new nodes
            if isfield(new_ele,'bm_x_neg') %&& element.ele_id(new_ele.bm_x_neg,1) ~= 16 && element.ele_id(new_ele.bm_x_neg,1) ~= 17
                element.node_2(new_ele.bm_x_neg,1) = new_node.id(1);
            end
            if isfield(new_ele,'bm_x_pos') %&& element.ele_id(new_ele.bm_x_pos,1) ~= 16 && element.ele_id(new_ele.bm_x_pos,1) ~= 17
                element.node_1(new_ele.bm_x_pos,1) = new_node.id(2);
            end
            if isfield(new_ele,'bm_z_neg')
                element.node_2(new_ele.bm_z_neg,1) = new_node.id(4);
            end
            if isfield(new_ele,'bm_z_pos')
                element.node_1(new_ele.bm_z_pos,1) = new_node.id(5);
            end
            element.node_2(new_ele.col_y,1) = new_node.id(3);
            clear new_ele
           
            % Define Joint
            joint_id = joint_id + 1;
            joint.id(joint_id,1) = joint_id + 1000;
            joint.x_neg(joint_id,1) = new_node.id(1);
            joint.x_pos(joint_id,1) = new_node.id(2);
            joint.y_neg(joint_id,1) = new_node.id(3);
            joint.y_pos(joint_id,1) = n_id;
            if strcmp(model.dimension,'3D')
                joint.z_neg(joint_id,1) = new_node.id(4);
                joint.z_pos(joint_id,1) = new_node.id(5);
            end
%             joint.center(joint_id,1) = new_node.id(4);
            
            % Add new nodes to nodes list
            node.id = [node.id; new_node.id];
            node.x = [node.x; new_node.x];
            node.y = [node.y; new_node.y];
            node.z = [node.z; new_node.z];
            node.dead_load = [node.dead_load; new_node.dead_load];
            node.live_load = [node.live_load; new_node.live_load];
            node.mass = [node.mass; new_node.mass];
        end
    end
end

%% Assign Additional Elements
for ae = 1:height(additional_elements)
    % Element Properties
    ele_id = ele_id + 1;
    ele = ele_props_table(ele_props_table.id == additional_elements.ele_id(ae),:);
    element.id(ele_id,1) = ele_id;
    element.trib_wt(ele_id,1) = 0;
    element.ele_id(ele_id,1) = ele.id;
    element.direction{ele_id,1} = additional_elements.direction{ae};
    element.story(ele_id,1) = 0;
    element.type{ele_id,1} = ele.type;
    element.dead_load(ele_id,1) = 0;
    element.live_load(ele_id,1) = 0;
    element.gravity_load(ele_id,1) = 0;
    
    % Check to see if the element nodes exists and assign
    [ node, id ] = fn_node_exist( node, additional_elements.x_start(ae), additional_elements.y_start(ae), additional_elements.z_start(ae), 1 );
    element.node_1(ele_id,1) = node.id(id);
    [ node, id ] = fn_node_exist( node, additional_elements.x_end(ae), additional_elements.y_end(ae), additional_elements.z_end(ae), 1 );
    element.node_2(ele_id,1) = node.id(id);
    element.node_3(ele_id,1) = 0;
    element.node_4(ele_id,1) = 0;
end

%% Find Element Lengths 
for e = 1:length(element.id)
    ele_x_end = node.x(node.id == element.node_2(e));
    ele_x_start = node.x(node.id == element.node_1(e));
    ele_y_end = node.y(node.id == element.node_2(e));
    ele_y_start = node.y(node.id == element.node_1(e));
    ele_z_end = node.z(node.id == element.node_2(e));
    ele_z_start = node.z(node.id == element.node_1(e));
    element.length(e,1) = sqrt( (ele_x_end-ele_x_start)^2 + (ele_y_end-ele_y_start)^2 + (ele_z_end-ele_z_start)^2 );
end

%% Find first nodes in each story and Nodes on Slab
node.story = zeros(length(node.id),1);
node.primary_story = zeros(length(node.id),1);
for s = 1:height(story)
    slab_ht = story.y_start(s) + story.story_ht(s);
    node.story(node.y == slab_ht) = s;
    node.primary_story(node.x == model.primary_node_offset & node.z == 0 & node.y == slab_ht) = 1;
    nodes_on_slab{s} = node.id(node.y == slab_ht);
end

%% Define new nodes to connect rigid diaphram to
node.on_slab = zeros(length(node.id),1);
if analysis.rigid_diaphram
    count = 0;
    last_node = node.id(end);
    for s = 1:height(story)
        for i = 1:length(nodes_on_slab{s})
            count = count + 1;
            node_id = last_node + count;
            node.id(node_id,:) = node_id;
            node.x(node_id,:) = node.x(node.id == nodes_on_slab{s}(i),:);
            node.y(node_id,:) = node.y(node.id == nodes_on_slab{s}(i),:);
            node.z(node_id,:) = node.z(node.id == nodes_on_slab{s}(i),:);
            node.dead_load(node_id,:) = 0;
            node.live_load(node_id,:) = 0;
            node.mass(node_id,:) = 0;
            node.story(node_id,:) = 0;
            node.primary_story(node_id,:) = 0;
            node.on_slab(node_id,:) = s;
        end
    end
end

%% Offset Mass for Accidental Torsion
% For X direction %%%% SHOULD CHANGE INTO FUNCTION AND RUN EACH DIR
% INDEPENDANTLY
if analysis.accidental_torsion == 1
    for i = 1:height(story)
        % Find Cental Coordinate and nodes to each side of it
        building_length = max(node.x(nodes_on_slab{i})) - min(node.x(nodes_on_slab{i}));
        node_center = building_length/2;
        nodes_right = nodes_on_slab{i}(node.x(nodes_on_slab{i}) > node_center);
        nodes_left = nodes_on_slab{i}(node.x(nodes_on_slab{i}) < node_center);
        mass_as_is = sum(node.mass(nodes_on_slab{i}));

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
            mass_new = sum(temp_mass(nodes_on_slab{i}));
            mass_diff = round(abs(mass_new - mass_as_is),2);
            if mass_diff > 0
                error('Mass got thrown off trying to find accidental torsion')
            end

            % Find new offset
            offset = (sum(temp_mass(nodes_on_slab{i}).*(node.x(nodes_on_slab{i}) - node_center))/sum(temp_mass(nodes_on_slab{i})))/building_length;
        end

        % Set Found Masses
        node.mass(nodes_right) = node.mass(nodes_right)*(1+mass_delta);
        node.mass(nodes_left) = node.mass(nodes_left)*(1-mass_delta);
    end

    % For Z direction
    for i = 1:height(story)
        % Find Cental Coordinate and nodes to each side of it
        building_length = max(node.z(nodes_on_slab{i})) - min(node.z(nodes_on_slab{i}));
        node_center = building_length/2;
        nodes_right = nodes_on_slab{i}(node.z(nodes_on_slab{i}) > node_center);
        nodes_left = nodes_on_slab{i}(node.z(nodes_on_slab{i}) < node_center);
        mass_as_is = sum(node.mass(nodes_on_slab{i}));

        % Iteratevely Increase the Mass till a target offset is reached
        offset = 0;
        offset_target = 0.05; % as a percentage of the building length
        threshold = 0.001;
        mass_delta = 0;
        while abs(offset_target-offset) > threshold
            mass_delta = mass_delta + 0.01;

            %Infinite loop happen safety device
            if mass_delta > 0.5
                error('Could Not Find Accidental Torsion Adjustment, Recalibrate')
            end

            temp_mass = node.mass;
            temp_mass(nodes_right) = node.mass(nodes_right)*(1+mass_delta);
            temp_mass(nodes_left) = node.mass(nodes_left)*(1-mass_delta);

            % Make sure mass still adds up to what it did before
            mass_new = sum(temp_mass(nodes_on_slab{i}));
            mass_diff = round(abs(mass_new - mass_as_is),2);
            if mass_diff > 0
                error('Mass got thrown off trying to find accidental torsion')
            end

            % Find new offset
            offset = (sum(temp_mass(nodes_on_slab{i}).*(node.z(nodes_on_slab{i}) - node_center))/sum(temp_mass(nodes_on_slab{i})))/building_length;
        end

        % Set Found Masses
        node.mass(nodes_right) = node.mass(nodes_right)*(1+mass_delta);
        node.mass(nodes_left) = node.mass(nodes_left)*(1-mass_delta);
    end
end

%% Define Nodal Fixity
node.fix = cell(length(node.id),1);
node.fix(1:end) = {'[000000]'}; % All nodes
foundation_nodes_id = node.id(node.y == 0);
% foundation nodes
if model.foundation == 1
    node.fix(foundation_nodes_id,:) = {'[111111]'}; % Fixed
elseif model.foundation == 0 
    node.fix(foundation_nodes_id,:) = {'[111000]'}; % Pinned
elseif model.foundation == 2
    node.fix(foundation_nodes_id,:) = {'[000000]'}; % Partial Fixity (ie pile hinge)
end

%% Create Nonlinear Rotational Springs at ends of all beams and columns, and shear springs at the bottom of walls
hinge.id = [];
hinge.element_id = [];
hinge.node_1 = [];
hinge.node_2 = [];
hinge_id = 0;
if analysis.nonlinear ~= 0
    % Define Hinges
    for i = 1:length(element.id)
        if element.story(i) <= analysis.stories_nonlinear
            if strcmp(element.type{i},'column') || strcmp(element.type{i},'beam') % For all columns and beams
                hinge_id = hinge_id+1;
                % Define hinge at start of element
                [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_1', i, hinge_id, foundation_nodes_id, 'rotational' ); 
                hinge_id = hinge_id+1;
                % Define hinge at end of element
                [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_2', i, hinge_id, foundation_nodes_id, 'rotational' );
            elseif strcmp(element.type{i},'wall') % For walls
                % Check if wall is the bottom wall
                walls_with_same_node = element.id(strcmp(element.type,'wall') & element.node_2 == element.node_1(i));
                if isempty(walls_with_same_node)
                    hinge_id = hinge_id+1;
                    % Define hinge at start of element
                    [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_1', i, hinge_id, foundation_nodes_id, 'shear' ); 
                end
            end
        end
    end
    
    %% Calculate the intial (zero axial load) Moment Capacity of Each element in the model
    for e = 1:length(element.id)
        ele = ele_props_table(ele_props_table.id == element.ele_id(e),:);
        if strcmp(ele.type,'truss') || contains(ele.description,'rigid')
            element.Mn_pos(e,1) = 9999999999;
            element.Mn_neg(e,1) = 9999999999;
        else
            % Moment Capcity per ACI
            [ ~, element.Mn_pos(e,1) ] = fn_aci_moment_capacity( 'pos', ele.fc_e, ele.w, ele.d, ele.As, ele.As_d, ele.fy_e, ele.Es, 0, ele.slab_depth, ele.b_eff ); % change to be based on the gravity load instead?
            [ ~, element.Mn_neg(e,1) ] = fn_aci_moment_capacity( 'neg', ele.fc_e, ele.w, ele.d, ele.As, ele.As_d, ele.fy_e, ele.Es, 0, ele.slab_depth, ele.b_eff );
        end
    end
else
    % Define Shear Spings on Walls
    for i = 1:length(element.id)
        if strcmp(element.type{i},'wall')
            walls_with_same_node = element.id(strcmp(element.type,'wall') & element.node_2 == element.node_1(i));
            if isempty(walls_with_same_node)
                hinge_id = hinge_id+1;
                % Define hinge at start of element
                [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_1', i, hinge_id, foundation_nodes_id, 'shear' ); 
            end
        end
    end
end

%% Define Foundation Hinges
if model.foundation == 2 % partial fixity such as pile hinge
    for f_node = 1:length(foundation_nodes_id)
        hinge_id = hinge_id+1;
        % Define hinge at foundation
        [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'NA', foundation_nodes_id(f_node), hinge_id, foundation_nodes_id, 'foundation' ); 
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
    new_node.story = node.story(non_z_nodes,:);
    new_node.primary_story = node.primary_story(non_z_nodes,:);
    node = new_node;
end

%% Reformat outputs to table and write CSV's
writetable(story,[write_dir filesep 'story.csv'])
node_table = struct2table(node);
writetable(node_table,[write_dir filesep 'node.csv'])
ele_table = struct2table(element);
writetable(ele_table,[write_dir filesep 'element.csv'])
joint_table = struct2table(joint);
writetable(joint_table,[write_dir filesep 'joint.csv'])
mf_joint_table = struct2table(mf_joint);
writetable(mf_joint_table,[write_dir filesep 'mf_joint.csv'])
hinge_table = struct2table(hinge);
writetable(hinge_table,[write_dir filesep 'hinge.csv'])

end

