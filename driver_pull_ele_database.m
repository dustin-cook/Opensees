% Script to pull info from ATC 134 run and format into element collection
% database

clear all
close all
clc

%% User inputs
analysis_dir = ['outputs' filesep 'simple_frame_and_wall_3D' filesep 'test' filesep 'asce_41_data'];

%% Import Packages
import asce_41.fn_define_backbone_rot

%% Load Inputs
ele_inputs = readtable('element_collection_inputs.csv','ReadVariableNames',true);
ele_props_table = readtable('inputs/element.csv','ReadVariableNames',true);
load([analysis_dir filesep 'joint_analysis.mat'])
load([analysis_dir filesep 'element_analysis.mat'])
load(['ground_motions' filesep 'ICSB_recordings' filesep 'recorded_edp_profile.mat'])

% Formulate Beam Table
beams = ele_inputs(strcmp(ele_inputs.ele_type,'beam'),:);
for i = 1:height(beams)
    ele = element(element.id == beams.id(i),:);
    if ~isempty(ele)
        ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
        if ele.pass_aci_dev_length == 1
            beams.ld{i} = 'adequate development length';
            beams.lb{i} = 'adequate development length';
        else
            beams.ld{i} = 'insufficient development length';
            beams.lb{i} = 'insufficient development length';
        end
        beams.row_ratio(i) = (ele_props.row - ele_props.row_prime)/ele.row_bal;
        beams.V_ratio_1(i) = ele.Vmax/(ele_props.w*ele_props.d*sqrt(ele_props.fc_e));
        beams.V_ratio_2(i) = ele.Vmax/ele.Vn;
        if strcmp(ele.trans_rein_check,'NC')
            beams.trans_rein(i) = 0;
        elseif strcmp(ele.trans_rein_check,'C')
            beams.trans_rein(i) = 1;
        end
        beams.S(i) = ele_props.S;
        beams.splices(i) = 0;
        beams.k0(i) = (1/1000)*6*ele_props.e*ele_props.iz/ele.length;
        beams.Mn{i} = num2str((1/1000)*[ele.Mn_pos, ele.Mn_neg]);
        beams.Mp{i} = num2str((1/1000)*[ele.Mp_pos, ele.Mp_neg]);
        if ele.pass_aci_dev_length == 0
            beams.condition(i) = 3;
        elseif strcmp(ele.critical_mode,'shear')
            beams.condition(i) = 2;
        else   
            beams.condition(i) = 1;
        end 
        beams.analysis_drift{i} = num2str([ele.drift_x, ele.drift_z]);
        if isfield(ele,'rot_1')
            beams.analysis_rotation{i} = num2str([ele.rot_1, ele.rot_2]);
        end
        beams.recording_drift{i} = num2str([record_edp.max_disp.x(ele.story+1), record_edp.max_disp.z(ele.story+1)]);
        beams.image{i} = 'NA';
        beams.a(i) = ele.a_hinge;
        beams.b(i) = ele.b_hinge;
        beams.c(i) = ele.c_hinge;
        beams.io(i) = ele.io;
        beams.ls(i) = ele.ls;
        beams.cp(i) = ele.cp;
    end
end

% Formulate Column Table
columns = ele_inputs(strcmp(ele_inputs.ele_type,'column'),:);
for i = 1:height(columns)
    ele = element(element.id == columns.id(i),:);
    if ~isempty(ele)
        ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
        columns.fce(i) = (1/1000)*ele_props.fc_e;
        columns.fyle(i) = (1/1000)*ele_props.fy_e;
        columns.fyte(i) = (1/1000)*ele_props.fy_e; % Assume they are the same
        columns.axial_max(i) = ele.Pmax/(ele_props.a*ele_props.fc_e);
        columns.axial_grav(i) = ele.P_grav/(ele_props.a*ele_props.fc_e);
        columns.V_ratio_1(i) = ele.vye/ele.V0;
        columns.V_ratio_2(i) = ele.Vmax/ele.V0;
        columns.row_t_ratio(i) = ele.rho_t;
        columns.row_ratio(i) = ele_props.row;
        columns.S(i) = ele_props.S;
        columns.trans_rein(i) = 1; % Assume the are anchored into the core based on 135 deg hooks
        columns.splices(i) = 0; % Plans say that splices are in center of column
        columns.k0(i) = 6*ele_props.e*ele_props.iz/ele.length*(1/1000);
        columns.Mn{i} = num2str([ele.Mn_pos, ele.Mn_neg]*(1/1000));
        columns.Mp{i} = num2str([ele.Mp_pos, ele.Mp_neg]*(1/1000));
        if ele.pass_aci_dev_length == 0
            columns.condition(i) = 2;
        else   
            columns.condition(i) = 1;
        end 
        columns.analysis_drift{i} = num2str([ele.drift_x, ele.drift_z]);
        if isfield(ele,'rot_1')
            columns.analysis_rotation{i} = num2str([ele.rot_1, ele.rot_2]);
        end
        columns.recording_drift{i} = num2str([record_edp.max_disp.x(ele.story+1), record_edp.max_disp.z(ele.story+1)]);
        columns.image{i} = ['Fig - ', columns.member_id{i}];
        columns.a(i) = ele.a_hinge;
        columns.b(i) = ele.b_hinge;
        columns.c(i) = ele.c_hinge;
        columns.io(i) = ele.io;
        columns.ls(i) = ele.ls;
        columns.cp(i) = ele.cp;
    end
end

% Formulate Joint Table
joints = ele_inputs(strcmp(ele_inputs.ele_type,'joint'),:);
joints.joint_id_1 = [];
joints.joint_id_2 = [];
for i = 1:height(joints)
    jnt = joint(joint.id == joints.id(i),:);
%     ele = element(element.id == joints.id(i),:);
%     ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
    if ~isempty(ele)
        joints.axial_max(i) = jnt.Pmax/(jnt.a*jnt.fc_e);
        joints.V_ratio_1(i) = jnt.Vmax/jnt.Vj;
        joints.scwb(i) = jnt.col_bm_ratio;
        if strcmp(jnt.trans_rien,'NC')
            joints.trans_rein(i) = 0; 
        else
            joints.trans_rein(i) = 1; 
        end
        [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', jnt.Mn(i), jnt.Mn(i), inf, inf, jnt.h(i), jnt.e(i), jnt.iz(i), jnt.a_hinge(i), jnt.b_hinge(i), jnt.c_hinge(i), NaN, 0.04 );
        joints.k0(i) = (1/1000)*1000*moment_vec_pos(1)/rot_vec_pos(1);
        joints.Mn{i} = num2str([moment_vec_pos(1), moment_vec_pos(1)]*(1/1000));
        joints.Mp{i} = num2str([moment_vec_pos(2), moment_vec_pos(2)]*(1/1000));
        if strcmp(jnt.class,'a') || strcmp(jnt.class,'b')
            joints.condition(i) = 1;
        else   
            joints.condition(i) = 2;
        end 
        joints.analysis_drift{i} = num2str([jnt.drift_x, jnt.drift_z]);
        joints.analysis_rotation(i) = 0; % Not modeling joint nonlinearity
        joints.recording_drift(i) = NaN; % Only worry about it for the columns
        joints.image{i} = 'NA';
        joints.a(i) = jnt.a_hinge;
        joints.b(i) = jnt.b_hinge;
        joints.c(i) = jnt.c_hinge;
        joints.io(i) = jnt.io;
        joints.ls(i) = jnt.ls;
        joints.cp(i) = jnt.cp;
    end
end

% Save Data
writetable(beams,[analysis_dir filesep 'damage_data_beams.csv'])
writetable(columns,[analysis_dir filesep 'damage_data_columns.csv'])
writetable(joints,[analysis_dir filesep 'damage_data_joints.csv'])