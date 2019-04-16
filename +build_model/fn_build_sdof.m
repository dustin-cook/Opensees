function [ ] = fn_build_sdof( model, analysis, write_dir )
% Funtion to build sdof element and node tables

%% INITIAL SETUP
% Import packages
import build_model.*

%% Build Story Table
story.id = 1;
story.story_ht = model.height;
story.y_start = 0;
story.story_dead_load = model.weight;


%% Build Node Table
node.id = [1;2];
node.x = [0;0];
node.y = [0; model.height];
node.dead_load = [0; model.weight];
node.live_load = [0;0];
node.mass = node.dead_load/386;
node.record_disp = [1;1];
node.record_accel = [1;1];
node.story = [0;1];
node.primary_story = [0;1];
node.fix = {'[111111]';'[000000]'};
node.on_slab = [0;0];


%% Build Element Table
element.id = 1;
element.ele_id = 0;
element.type = 'column';
element.node_1 = 1;
element.node_2 = 2;
element.length = model.height;
element.direction = 'x';
element.k = (model.weight/386)/(model.period/(2*pi))^2;
element.e = 1000000;
element.a = 1000000;
element.iz = element.k*(element.length^3)/(3*element.e);
element.Mn_pos_1 = model.moment_capacity;
element.Mn_neg_1 = model.moment_capacity;
element.Mp_pos_1 = model.moment_capacity;
element.Mp_neg_1 = model.moment_capacity;
element.a_hinge_1 = 0.01;
element.b_hinge_1 = 0.02;
element.c_hinge_1 = 0.2;
element.critical_mode_1 = 'flexure';

%% Assign Joints
joint.id = [];

%% Create Nonlinear Rotational Springs at ends of all beams and columns
hinge.id = [];
hinge.element_id = [];
hinge.node_1 = [];
hinge.node_2 = [];
if analysis.nonlinear ~= 0
    % Define hinge at start of element
    [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_1', 1, 1, 1, 'rotational', 'primary', 1 ); 
end

%% Reformat outputs to table and write CSV's
story_table = struct2table(story);
writetable(story_table,[write_dir filesep 'story.csv'])
node_table = struct2table(node);
writetable(node_table,[write_dir filesep 'node.csv'])
ele_table = struct2table(element);
writetable(ele_table,[write_dir filesep 'element.csv'])
joint_table = struct2table(joint);
writetable(joint_table,[write_dir filesep 'joint.csv'])
hinge_table = struct2table(hinge);
writetable(hinge_table,[write_dir filesep 'hinge.csv'])

end

