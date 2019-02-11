% Script to pull info from ATC 134 run and format into element collection
% database

clear all
close all
clc

%% User inputs
analysis_dir = ['outputs' filesep 'ICBS_model_3D_fixed' filesep 'NDP' filesep 'asce_41_data'];

%% Import Packages
import asce_41.fn_define_backbone_rot
import asce_41.fn_define_backbone_shear

%% Load Inputs
ele_inputs = readtable('element_collection_key.csv','ReadVariableNames',true);
ele_props_table = readtable('inputs/element.csv','ReadVariableNames',true);
load([analysis_dir filesep 'joint_analysis.mat'])
load([analysis_dir filesep 'element_analysis.mat'])
load([analysis_dir filesep 'story_analysis.mat'])
load(['ground_motions' filesep 'ICSB_recordings' filesep 'recorded_edp_profile.mat'])

% Formulate Beam Table
beams = ele_inputs(strcmp(ele_inputs.element_type,'beam'),:);
for i = 1:height(beams)
    ele = element(element.id == beams.id(i),:);
    if ~isempty(ele)
        this_story = story(story.id == ele.story,:);
        ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
        beams.ld_req(i) = ele.ld_req;
        beams.ld_avail{i} = ele.ld_avail;
        beams.rho_ratio(i) = (ele_props.row - ele_props.row_prime)/ele.row_bal;
        beams.shear_stress_ratio(i) = ele.Vmax/(ele_props.w*ele_props.d*sqrt(ele_props.fc_e));
        beams.shear_demand_ratio(i) = ele.Vmax/ele.Vn;
        if strcmp(ele.trans_rein_check,'NC')
            beams.trans_rein(i) = 0;
        elseif strcmp(ele.trans_rein_check,'C')
            beams.trans_rein(i) = 1;
        end
        beams.s(i) = ele_props.S;
        beams.hinge_splices(i) = 0;
        beams.K0(i) = (1/1000)*6*ele_props.e*ele_props.iz/ele.length;
        if strcmp(ele.rot_1_dir,'pos')
            beams.Mn_1{i} = (1/1000)*ele.Mn_pos;
            beams.Mp_1{i} = (1/1000)*ele.Mp_pos;
        elseif strcmp(ele.rot_1_dir,'neg')
            beams.Mn_1{i} = (1/1000)*ele.Mn_neg;
            beams.Mp_1{i} = (1/1000)*ele.Mp_neg;
        end
        if strcmp(ele.rot_2_dir,'pos')
            beams.Mn_2{i} = (1/1000)*ele.Mn_pos;
            beams.Mp_2{i} = (1/1000)*ele.Mp_pos;
        elseif strcmp(ele.rot_2_dir,'neg')
            beams.Mn_2{i} = (1/1000)*ele.Mn_neg;
            beams.Mp_2{i} = (1/1000)*ele.Mp_neg;
        end
        if ele.pass_aci_dev_length == 0
            beams.condition(i) = 3;
        elseif strcmp(ele.critical_mode,'shear')
            beams.condition(i) = 2;
        else   
            beams.condition(i) = 1;
        end 
        beams.max_story_drift_analysis(i) = this_story.(['max_drift_' ele.direction{1}]); % Need to fix elment story drift read such that it does not come in as NaN
%         beams.max_story_drift_analysis(i) = ele.(['drift_' ele.direction{1}]);
        if sum(strcmp('rot_1',ele.Properties.VariableNames)) > 0
            beams.max_element_rotation_or_drift_1(i) = ele.rot_1;
        end
        if sum(strcmp('rot_2',ele.Properties.VariableNames)) > 0
            beams.max_element_rotation_or_drift_2(i) = ele.rot_2;
        end
        beams.max_story_drift_record(i) = record_edp.max_disp.(ele.direction{1})(ele.story+1);
        beams.damage_image{i} = 'NA';
        beams.a(i) = ele.a_hinge;
        beams.b(i) = ele.b_hinge;
        beams.c(i) = ele.c_hinge;
        beams.io(i) = ele.io;
        beams.ls(i) = ele.ls;
        beams.cp(i) = ele.cp;
    end
end

% Formulate Column Table
columns = ele_inputs(strcmp(ele_inputs.element_type,'column'),:);
for i = 1:height(columns)
    ele = element(element.id == columns.id(i),:);
    if ~isempty(ele)
        this_story = story(story.id == ele.story,:);
        ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
        columns.fc_e(i) = (1/1000)*ele_props.fc_e;
        columns.fy_e(i) = (1/1000)*ele_props.fy_e;
        columns.fyt_e(i) = (1/1000)*ele_props.fy_e; % Assume they are the same
        columns.max_axial_load_ratio(i) = ele.Pmax/(ele_props.a*ele_props.fc_e);
        columns.gravity_axial_load_ratio(i) = ele.P_grav/(ele_props.a*ele_props.fc_e);
        columns.shear_flexure_yield_ratio(i) = ele.vye/ele.V0;
        columns.shear_demand_ratio(i) = ele.Vmax/ele.V0;
        columns.rho_t(i) = ele.rho_t;
        columns.rho_l(i) = ele_props.row;
        columns.s(i) = ele_props.S;
        columns.hoops_conform(i) = 1; % Assume the are anchored into the core based on 135 deg hooks
        columns.hinge_splices(i) = 0; % Plans say that splices are in center of column
        columns.K0(i) = 6*ele_props.e*ele_props.iz/ele.length*(1/1000);
        columns.Mn{i} = ele.Mn_pos*(1/1000);
        columns.Mp{i} = ele.Mp_pos*(1/1000);
        if ele.pass_aci_dev_length == 0
            columns.condition(i) = 2;
        else   
            columns.condition(i) = 1;
        end 
        columns.max_story_drift_analysis(i) = this_story.(['max_drift_' ele.direction{1}]);
%         columns.max_story_drift_analysis(i) = ele.(['drift_' ele.direction{1}]);
        if sum(strcmp('rot_1',ele.Properties.VariableNames)) > 0
            columns.max_element_rotation_or_drift_1(i) = ele.rot_1;
        end
        if sum(strcmp('rot_2',ele.Properties.VariableNames)) > 0
            columns.max_element_rotation_or_drift_2(i) = ele.rot_2;
        end
        columns.max_story_drift_record(i) = record_edp.max_disp.(ele.direction{1})(ele.story+1);
        columns.damage_image{i} = ['Fig - ', columns.member_id{i}];
        columns.a(i) = ele.a_hinge;
        columns.b(i) = ele.b_hinge;
        columns.c(i) = ele.c_hinge;
        columns.io(i) = ele.io;
        columns.ls(i) = ele.ls;
        columns.cp(i) = ele.cp;
    end
end

% Formulate Joint Table
joints = ele_inputs(strcmp(ele_inputs.element_type,'joint'),:);
joints.joint_id_1 = [];
joints.joint_id_2 = [];
for i = 1:height(joints)
    jnt = joint(joint.id == joints.id(i),:);
%     ele = element(element.id == joints.id(i),:);
%     ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
    if ~isempty(jnt)
        joints.axial_load_ratio(i) = jnt.Pmax/(jnt.a*jnt.fc_e);
        joints.shear_ratio(i) = jnt.Vmax/jnt.Vj;
        joints.scwb(i) = jnt.col_bm_ratio;
        if strcmp(jnt.trans_rien,'NC')
            joints.trans_rein(i) = 0; 
        else
            joints.trans_rein(i) = 1; 
        end
        [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', jnt.Mn, jnt.Mn, inf, inf, jnt.h, jnt.e, jnt.iz, jnt.a_hinge, jnt.b_hinge, jnt.c_hinge, NaN, 0.04 );
        joints.K0(i) = (1/1000)*1000*moment_vec_pos(1)/rot_vec_pos(1);
        joints.Mn{i} = moment_vec_pos(1)*(1/1000);
        joints.Mp{i} = moment_vec_pos(2)*(1/1000);
        if strcmp(jnt.class,'a') || strcmp(jnt.class,'b')
            joints.condition(i) = 1;
        else   
            joints.condition(i) = 2;
        end 
        joints.max_story_drift_analysis(i) = ele.(['drift_' 'x']); % Hardcode as X for now, all joints are in x direction
        joints.max_element_rotation_or_drift(i) = 0; % Not modeling joint nonlinearity
        joints.max_story_drift_record(i) = NaN; % Only worry about it for the columns
        joints.damage_image{i} = 'NA';
        joints.a(i) = jnt.a_hinge;
        joints.b(i) = jnt.b_hinge;
        joints.c(i) = jnt.c_hinge;
        joints.io(i) = jnt.io;
        joints.ls(i) = jnt.ls;
        joints.cp(i) = jnt.cp;
    end
end

% Formulate Wall Table
walls = ele_inputs(strcmp(ele_inputs.element_type,'wall'),:);
walls.joint_id_1 = [];
walls.joint_id_2 = [];
for i = 1:height(walls)
    ele = element(element.id == walls.id(i),:);
    if ~isempty(ele)
        this_story = story(story.id == ele.story,:);
        ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
        As = sum(str2double(strsplit(strrep(strrep(ele_props.As{1},']',''),'[',''))))/2; % Assume 1/2 tensiona and 1/2 compression just like in wall hinge selection
        As_prime = As;
        walls.axial_load_ratio(i) = ((As-As_prime)*ele_props.fy_e+ele.Pmax)/(ele_props.w*ele_props.d*ele_props.fc_e);
        if strcmp(ele.critical_mode,'shear')
            walls.condition(i) = 2;
        elseif strcmp(ele.critical_mode,'flexure') 
            walls.condition(i) = 1;
        end 
        walls.shear_stress_ratio(i) = ele.Vmax/(ele_props.w*ele_props.d*sqrt(ele_props.fc_e));
        walls.confined_boundary(i) = 0; % Assume they are unconfined (doens't matter for shear controlled walls)
        [ force_vec, disp_vec ] = fn_define_backbone_shear( ele.Vn, ele.length, ele_props.g, ele_props.av, ele ); % Assumes walls are shear controlled, probably should check explicitly
        walls.K0(i) = (1/1000)*force_vec(1)/disp_vec(1);
        walls.Mn{i} = force_vec(2)*(1/1000);
        walls.Mp{i} = force_vec(3)*(1/1000);
        if sum(strcmp(['drift_' ele.direction{1}],ele.Properties.VariableNames)) > 0
            walls.max_story_drift_analysis(i) = ele.(['drift_' ele.direction{1}]);
        else
            walls.max_story_drift_analysis(i) = 0;
        end
        if sum(strcmp('rot_1',ele.Properties.VariableNames)) > 0
            walls.max_element_rotation_or_drift(i) = ele.rot_1;
        end
        walls.max_story_drift_record(i) = record_edp.max_disp.(ele.direction{1})(ele.story+1);
        walls.damage_image{i} = 'NA';
        walls.a(i) = ele.a_hinge;
        walls.b(i) = ele.b_hinge;
        walls.c(i) = ele.c_hinge;
        walls.d(i) = ele.d_hinge;
        walls.e(i) = ele.e_hinge;
        walls.f(i) = ele.f_hinge;
        walls.g(i) = ele.g_hinge;
        walls.io(i) = ele.io;
        walls.ls(i) = ele.ls;
        walls.cp(i) = ele.cp;
    end
end

% Save Data
if ~isempty(beams)
    writetable(beams,[analysis_dir filesep 'raw_beams_database.csv'])
end
if ~isempty(columns)
    writetable(columns,[analysis_dir filesep 'raw_columns_database.csv'])
end
if ~isempty(joints)
    writetable(joints,[analysis_dir filesep 'raw_joints_database.csv'])
end
if ~isempty(walls)
    writetable(walls,[analysis_dir filesep 'raw_walls_database.csv'])
end