function [ ] = fn_build_mdof( model, ele_props_table, analysis, write_dir, read_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% INITIAL SETUP
% Import Packages
import build_model.*
import aci_318.*

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
    if ~isempty(story_group)
        for g = 1:height(story_group)
            element_group = element_group_table((element_group_table.id == story_group.element_group_id(g)),:);
            for e = 1:height(element_group)
                for nb = 1:element_group.num_bays
                    % Element Properties Column
                    if iscell(element_group.col_id)
                        element_group.col_id = str2double(strsplit(strrep(strrep(element_group.col_id{1},'[',''),']',''),','));
                    end
                    if iscell(element_group.beam_id)
                        element_group.beam_id = str2double(strsplit(strrep(strrep(element_group.beam_id{1},'[',''),']',''),','));
                    end
                    if iscell(element_group.wall_id)
                        element_group.wall_id = str2double(strsplit(strrep(strrep(element_group.wall_id{1},'[',''),']',''),','));
                    end
                    if iscell(element_group.trib_wt_1)
                        element_group.trib_wt_1 = str2double(strsplit(strrep(strrep(element_group.trib_wt_1{1},'[',''),']',''),','));
                    end
                    if iscell(element_group.trib_wt_2)
                        element_group.trib_wt_2 = str2double(strsplit(strrep(strrep(element_group.trib_wt_2{1},'[',''),']',''),','));
                    end

                    if element_group.col_id(nb) ~= 0
                        ele_id = ele_id + 1;
                        [ node, element ] = fn_create_element( 'col', ele_id, ele_props_table, element_group.col_id(nb), nb, story_props, story_group(g,:), node, element, story_group.direction{g}, element_group.trib_wt_1(nb), element_group.trib_wt_2(nb) );
                    end
                    if element_group.beam_id(nb) ~= 0
                        ele_id = ele_id + 1;
                        [ node, element ] = fn_create_element( 'beam', ele_id, ele_props_table, element_group.beam_id(nb), nb, story_props, story_group(g,:), node, element, story_group.direction{g}, element_group.trib_wt_1(nb), element_group.trib_wt_2(nb) );
                    end
                    if element_group.wall_id(nb) ~= 0
                        ele_id = ele_id + 1;
                        [ node, element ] = fn_create_element( 'wall', ele_id, ele_props_table, element_group.wall_id(nb), nb, story_props, story_group(g,:), node, element, story_group.direction{g}, element_group.trib_wt_1(nb), element_group.trib_wt_2(nb) );
                    end
                end
                % Last column in bay span
                if sum(element_group.col_id ~=0) ~= 0 && element_group.col_id(nb+1) ~= 0
                    ele_id = ele_id + 1;
                    [ node, element ] = fn_create_element( 'col', ele_id, ele_props_table, element_group.col_id(nb+1), element_group.num_bays+1, story_props, story_group(g,:), node, element, story_group.direction{g} );
                end
            end
        end

        % Assign Gravity Load to Elements
        sum_trib_wt = sum([element.trib_wt_1(element.story == story.id(s));element.trib_wt_2(element.story == story.id(s))]);
        element.dead_load_1(element.story == story.id(s),1) = element.trib_wt_1(element.story == story.id(s))*story.story_dead_load(s)/sum_trib_wt;
        element.live_load_1(element.story == story.id(s),1) = element.trib_wt_1(element.story == story.id(s))*story.story_live_load(s)/sum_trib_wt;
        element.dead_load_2(element.story == story.id(s),1) = element.trib_wt_2(element.story == story.id(s))*story.story_dead_load(s)/sum_trib_wt;
        element.live_load_2(element.story == story.id(s),1) = element.trib_wt_2(element.story == story.id(s))*story.story_live_load(s)/sum_trib_wt;
    end
end

%% Define Joints (For each node created go through as say whats connected in)
joint_id = 0;
joint.id = [];
for i = 1:length(node.id)
    eles_at_node = element.id(((element.node_1 == node.id(i)) | (element.node_2 == node.id(i))));
    if length(eles_at_node) > 1 && ~strcmp(element.type(element.id == eles_at_node(1)),'wall')
        joint_id = joint_id + 1;
        joint.id(joint_id,1) = joint_id;
        joint.column_low(joint_id,1) = 0;
        joint.column_high(joint_id,1) = 0;
        joint.beam_left(joint_id,1) = 0;
        joint.beam_right(joint_id,1) = 0;
        for e = 1:length(eles_at_node)
            if strcmp(element.type(element.id == eles_at_node(e)),'column')
                if element.node_1(element.id == eles_at_node(e)) == node.id(i)
                    joint.column_high(joint_id,1) = eles_at_node(e);
                elseif element.node_2(element.id == eles_at_node(e)) == node.id(i)
                    joint.column_low(joint_id,1) = eles_at_node(e);
                end
            elseif strcmp(element.type(element.id == eles_at_node(e)),'beam')
                if element.node_1(element.id == eles_at_node(e)) == node.id(i)
                    joint.beam_right(joint_id,1) = eles_at_node(e);
                elseif element.node_2(element.id == eles_at_node(e)) == node.id(i)
                    joint.beam_left(joint_id,1) = eles_at_node(e);
                end
            end
        end
    end
end

%% Calculate Total Element Gravity Load
element.dead_load = element.dead_load_1 + element.dead_load_2;
element.live_load = element.live_load_1 + element.live_load_2;
element.gravity_load = element.dead_load*analysis.dead_load + element.live_load*analysis.live_load;

%% Define Mass
node.mass = zeros(length(node.id),1);
for e = 1:length(element.id)
    node.mass(node.id == element.node_1(e)) = node.mass(node.id == element.node_1(e)) + element.dead_load_1(e)/386;
    node.mass(node.id == element.node_2(e)) = node.mass(node.id == element.node_2(e)) + element.dead_load_2(e)/386;
end

%% Set existing nodes as the ones that get recorded
node.record_disp = ones(length(node.id),1);
node.record_accel = ones(length(node.id),1);

%% Assign Joints
joint_id = 0;
for s = 1:height(story)
    this_story = story.id(s);
    node_y = story.y_start(s)+story.story_ht(s);
    if node_y > 0
        story_node = node.id(node.y == node_y);
        for n = 1:length(story_node)
            n_id = story_node(n);
            elements_at_node = element.id(((element.node_1 == n_id) | (element.node_2 == n_id)) & (element.story == this_story));
            if length(elements_at_node) > 1

                bm_d_x = 0;
                bm_d_z = 0;
                col_d_x = 0;
                col_d_z = 0;
                for e = 1:length(elements_at_node)
                    e_id = elements_at_node(e);
                    ele = ele_props_table(ele_props_table.id == element.ele_id(e_id),:);

                    % Find Joint properties based on elements that frame into the joint
                    if strcmp(element.type{e_id},'column')
                        new_ele.col_y = e_id;
                        if strcmp(element.direction{e_id},'x')
                            col_d_x(e) = ele.h;
                            col_d_z(e) = ele.w;
                        elseif strcmp(element.direction{e_id},'z')
                            col_d_x(e) = ele.w;
                            col_d_z(e) = ele.h;
                        else
                            error('element direction not recognized')
                        end
                    elseif strcmp(element.type{e_id},'beam')
                        if strcmp(element.direction{e_id},'x')
                            bm_d_x(e) = ele.h;
                            if element.node_1(e_id) == n_id
                                new_ele.bm_x_pos = e_id;
                            else
                                new_ele.bm_x_neg = e_id;
                            end
                        elseif strcmp(element.direction{e_id},'z')
                            bm_d_z(e) = ele.h;
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

                % Define Joints and joint Nodes for 2D representaton of the joint in the X direction Change elements to connect to new nodes
                joint_id = joint_id + 1;
                joint_node_id =  1 + node.id(end);
                new_node.id = joint_node_id;
                new_node.x = node.x(node.id == n_id);
                new_node.y = node.y(node.id == n_id)-joint_dim_y;
                new_node.z = node.z(node.id == n_id);
                new_node.mass = 0;
                new_node.record_disp = 0;
                new_node.record_accel = 0;
                element.node_2(new_ele.col_y,1) = new_node.id(1); % This breaks the joint 3D option
                joint.y_neg(joint_id,1) = new_node.id(1);
                joint.y_pos(joint_id,1) = n_id;  

                if isfield(new_ele,'bm_x_neg') %&& element.ele_id(new_ele.bm_x_neg,1) ~= 16 && element.ele_id(new_ele.bm_x_neg,1) ~= 17
                    joint_node_id =  joint_node_id + 1;
                    new_node.id = [new_node.id; joint_node_id];
                    new_node.x = [new_node.x; node.x(node.id == n_id)-joint_dim_x/2];
                    new_node.y = [new_node.y; node.y(node.id == n_id)-joint_dim_y/2];
                    new_node.z = [new_node.z; node.z(node.id == n_id);];
                    new_node.mass = [new_node.mass; 0];
                    new_node.record_disp = [new_node.record_disp; 0];
                    new_node.record_accel = [new_node.record_accel; 0];
                    element.node_2(new_ele.bm_x_neg,1) = new_node.id(end);
                    joint.x_neg(joint_id,1) = new_node.id(end);
                else
                    joint.x_neg(joint_id,1) = 0;
                end

                if isfield(new_ele,'bm_x_pos') %&& element.ele_id(new_ele.bm_x_pos,1) ~= 16 && element.ele_id(new_ele.bm_x_pos,1) ~= 17
                    joint_node_id =  joint_node_id + 1;
                    new_node.id = [new_node.id; joint_node_id];
                    new_node.x = [new_node.x; node.x(node.id == n_id)+joint_dim_x/2];
                    new_node.y = [new_node.y; node.y(node.id == n_id)-joint_dim_y/2];
                    new_node.z = [new_node.z; node.z(node.id == n_id)];
                    new_node.mass = [new_node.mass; 0];
                    new_node.record_disp = [new_node.record_disp; 0];
                    new_node.record_accel = [new_node.record_accel; 0];
                    element.node_1(new_ele.bm_x_pos,1) = new_node.id(end);
                    joint.x_pos(joint_id,1) = new_node.id(end);
                else
                    joint.x_pos(joint_id,1) = new_node.id(end);
                end           

                % Define Joint Classification according to ASCE 41-17 figure 10-3
                % ASSUMES THERE ARE NEVER ANY TRANSFER BEAMS (ie in the z
                % direction) AND THAT THE PRIMARY DIRECTION IS ALWAYS THE X.
                % ALSO THAT THERE WILL ALWAYS BE A COLUMN ABOVE UNLESS ITS THE
                % TOP STORY.
                if isfield(new_ele,'bm_x_pos') && isfield(new_ele,'bm_x_neg')
                    joint.class{joint_id,1} = 'b'; % interior joint without transfer beams
                elseif s == height(story) % top story
                    joint.class{joint_id,1} = 'e'; % knee joint with or without transfer beams
                else
                    joint.class{joint_id,1} = 'd'; % exterior joint without transfer beams
                end


                % Add new nodes to nodes list
                node.id = [node.id; new_node.id];
                node.x = [node.x; new_node.x];
                node.y = [node.y; new_node.y];
                node.z = [node.z; new_node.z];
                node.mass = [node.mass; new_node.mass];
                node.record_disp = [node.record_disp; new_node.record_disp];
                node.record_accel = [node.record_accel; new_node.record_accel];

                % Clear data
                clear new_ele
            end
        end
    end
end

%% Assign Additional Elements
element.elastic = zeros(length(element.id),1);
for ae = 1:height(additional_elements)
    % Element Properties
    ele_id = ele_id + 1;
    ele = ele_props_table(ele_props_table.id == additional_elements.ele_id(ae),:);
    element.id(ele_id,1) = ele_id;
    element.trib_wt_1(ele_id,1) = 0;
    element.trib_wt_2(ele_id,1) = 0;
    element.ele_id(ele_id,1) = ele.id;
    element.direction{ele_id,1} = additional_elements.direction{ae};
    element.story(ele_id,1) = additional_elements.story(ae);
    element.type{ele_id,1} = ele.type;
    element.dead_load_1(ele_id,1) = 0;
    element.live_load_1(ele_id,1) = 0;
    element.dead_load_2(ele_id,1) = 0;
    element.live_load_2(ele_id,1) = 0;
    element.dead_load(ele_id,1) = 0;
    element.live_load(ele_id,1) = 0;
    element.gravity_load(ele_id,1) = 0;
    element.elastic(ele_id,1) = 1;
    
    % Check to see if the element nodes exists and assign
    [ node, id ] = fn_node_exist( node, additional_elements.x_start(ae), additional_elements.y_start(ae), additional_elements.z_start(ae), 0 );
    element.node_1(ele_id,1) = node.id(id);
    [ node, id ] = fn_node_exist( node, additional_elements.x_end(ae), additional_elements.y_end(ae), additional_elements.z_end(ae), 0 );
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

%% Define if elements are rigid or not
for e = 1:length(element.id)
    ele_props = ele_props_table(ele_props_table.id == element.ele_id(e),:);
    if contains(ele_props.description,'rigid')
        element.rigid(e,1) = 1;
    else
        element.rigid(e,1) = 0;
    end
        
end

%% Find first nodes in each story and Nodes on Slab to connect to rigid diaphram
node.story = zeros(length(node.id),1);
node.on_slab = zeros(length(node.id),1);
node.on_slab(node.y == 0) = 1; % Set base nodes to be on the slab
node.primary_story = zeros(length(node.id),1);
for s = 1:height(story)
    slab_ht = story.y_start(s) + story.story_ht(s);
    node.story(node.y == slab_ht) = story.id(s);
    node.primary_story(node.x == model.primary_node_offset & node.z == 0 & node.y == slab_ht) = 1;
    node.on_slab(node.y == slab_ht) = 1;
end

%% Offset Mass for Accidental Torsion
% For X direction %%%% SHOULD CHANGE INTO FUNCTION AND RUN EACH DIR
% INDEPENDANTLY
if strcmp(model.dimension,'3D') && analysis.accidental_torsion == 1
    for s = 1:height(story)
        slab_ht = story.y_start(s) + story.story_ht(s);
        if slab_ht > 0
            % Direction x
            [ node ] = fn_update_mass( node, node.id(node.y == slab_ht), 'x' );
            % Direction z
            [ node ] = fn_update_mass( node, node.id(node.y == slab_ht), 'z' );
        end
    end
end

%% Define Nodal Fixity
node.fix = cell(length(node.id),1);
node.fix(1:end) = {'[000000]'}; % All nodes
foundation_nodes_filter = node.y == 0;
foundation_nodes_ids = node.id(foundation_nodes_filter);
wall_base_nodes_filter = node.y == 0 & ismember(node.id,element.node_1(strcmp(element.type,'wall')));
wall_foundation_nodes_ids = node.id(wall_base_nodes_filter);
% foundation nodes
if model.foundation == 1
    node.fix(foundation_nodes_filter,:) = {'[111111]'}; % Fixed
elseif model.foundation == 0 
    node.fix(foundation_nodes_filter,:) = {'[111000]'}; % Pinned
elseif model.foundation == 2
    node.fix(foundation_nodes_filter,:) = {'[000000]'}; % Partial Fixity (ie pile hinge)
end

%% Create Nonlinear Rotational Springs at ends of all beams and columns, and shear springs at the bottom of walls
hinge.id = [];
[ node, element, hinge ] = fn_define_hinge( analysis, model, hinge, element, node, foundation_nodes_filter );

%% Define Foundation Hinges
if model.foundation == 2 % partial fixity such as pile hinge
    for f = 1:length(foundation_nodes_ids)
        hinge_id = hinge.id(end)+1;
        if sum(ismember(wall_foundation_nodes_ids,foundation_nodes_ids(f))) > 0
            direction_str = 'wall'; % base of a wall foundation
        else
            direction_str = 'pile'; % base of a column foundation
        end
        % Define hinge at foundation
        [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'NA', foundation_nodes_ids(f), hinge_id, foundation_nodes_filter, 'foundation', direction_str ); 
    end
end

%% Convert outputs to tables
node = struct2table(node);
element = struct2table(element);
joint = struct2table(joint);
hinge = struct2table(hinge);

%% Define Center Nodes
node.center = zeros(length(node.id),1);
for s = 1:height(story)
    this_story = story.id(s);
    story_nodes = node(node.story == this_story & (node.record_accel == 1 | node.record_disp == 1),:);
    x_center = mean([min(story_nodes.x),max(story_nodes.x)]);
    z_center = mean([min(story_nodes.z),max(story_nodes.z)]);
    node_resultant_dist = sqrt((story_nodes.x-x_center).^2 + (story_nodes.z-z_center).^2);
    [~,closest_idx] = min(node_resultant_dist);
    center_node_id = story_nodes.id(closest_idx);
    node.center(node.id == center_node_id) = 1;
end

%% Reformat outputs to table and write CSV's
writetable(story,[write_dir filesep 'story.csv'])
writetable(node,[write_dir filesep 'node.csv'])
writetable(element,[write_dir filesep 'element.csv'])
writetable(joint,[write_dir filesep 'joint.csv'])
writetable(hinge,[write_dir filesep 'hinge.csv'])

end

