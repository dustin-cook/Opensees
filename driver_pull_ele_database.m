% Script to pull info from ATC 134 run and format into element collection
% database

clear all
close all
clc

%% User inputs
analysis_dir = ['outputs' filesep 'simple_frame_and_wall_3D' filesep 'test' filesep 'asce_41_data'];

%% Load Inputs
ele_inputs = readtable('element_collection_inputs.csv','ReadVariableNames',true);
ele_props_table = readtable('inputs/element.csv','ReadVariableNames',true);
load([analysis_dir filesep 'joint_analysis.mat'])
load([analysis_dir filesep 'element_analysis.mat'])

% Formulate Beam Table
beams = ele_inputs(strcmp(ele_inputs.ele_type,'beam'),:);
for i = 1:height(beams)
    ele = element(element.id == beams.id(i),:);
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
    beams.analysis_drift(i) = 0; % LOOK UP EXPLICITY
    beams.analysis_rotation(i) = 0; % LOOK UP EXPLICITY
    beams.recording_drift(i) = 0; % LOOK UP EXPLICITY
    beams.image{i} = 'NA';
    beams.a(i) = ele.a_hinge;
    beams.b(i) = ele.b_hinge;
    beams.c(i) = ele.c_hinge;
    beams.io(i) = ele.io;
    beams.ls(i) = ele.ls;
    beams.cp(i) = ele.cp;
end

% Formulate Column Table
columns = ele_inputs(strcmp(ele_inputs.ele_type,'column'),:);
for i = 1:height(columns)
    ele = element(element.id == columns.id(i),:);
    ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
    columns.fce(i) = (1/1000)*ele_props.fc_e;
    columns.fyle(i) = (1/1000)*ele_props.fy_e;
    columns.fyte(i) = (1/1000)*ele_props.fy_e; % LOOK UP EXPLICITY FROM PLANS
    columns.axial_max(i) = ele.Pmax/(ele_props.a*ele_props.fc_e);
    columns.axial_grav(i) = ele.P_grav/(ele_props.a*ele_props.fc_e);
    columns.V_ratio_1(i) = ele.vye/ele.V0; % DOUBLE CHECK V0 is what we are looking for
    columns.V_ratio_2(i) = ele.Vmax/ele.V0; % DOUBLE CHECK V0 is what we are looking for
    columns.row_t_ratio(i) = ele_props.Av/(ele_props.S*ele_props.w); % LOOK UP EXPLICITY
    columns.row_ratio(i) = ele_props.row;
    columns.S(i) = ele_props.S;
    if strcmp(ele.critical_mode,'shear')
        columns.trans_rein(i) = 0; % LOOK UP EXPLICITY
    else
        columns.trans_rein(i) = 1; % LOOK UP EXPLICITY
    end
    columns.splices(i) = 0;
    columns.k0(i) = 6*ele_props.e*ele_props.iz/ele.length*(1/1000);
    columns.Mn{i} = num2str([ele.Mn_pos, ele.Mn_neg]*(1/1000));
    columns.Mp{i} = num2str([ele.Mp_pos, ele.Mp_neg]*(1/1000));
    if ele.pass_aci_dev_length == 0
        columns.condition(i) = 2;
    else   
        columns.condition(i) = 1;
    end 
    columns.analysis_drift(i) = 0; % LOOK UP EXPLICITY
    columns.analysis_rotation(i) = 0; % LOOK UP EXPLICITY
    columns.recording_drift(i) = 0; % LOOK UP EXPLICITY
    columns.image{i} = ['Fig - ', columns.member_id{i}];
    columns.a(i) = ele.a_hinge;
    columns.b(i) = ele.b_hinge;
    columns.c(i) = ele.c_hinge;
    columns.io(i) = ele.io;
    columns.ls(i) = ele.ls;
    columns.cp(i) = ele.cp;
end

% Formulate Joint Table
joints = ele_inputs(strcmp(ele_inputs.ele_type,'joint'),:);
joints.joint_id_1 = [];
joints.joint_id_2 = [];
for i = 1:height(joints)
    jnt = joint(joint.id == joints.id(i),:);
%     ele = element(element.id == joints.id(i),:);
%     ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
    joints.axial_max(i) = jnt.Pmax/(jnt.a*jnt.fc_e); % DOUBLE CHECK AG IS WHAT WE WANT IT TO BE
    joints.V_ratio_1(i) = jnt.Vmax/jnt.Vj;
    joints.scwb(i) = jnt.col_bm_ratio;
    if strcmp(jnt.trans_rien,'NC')
        joints.trans_rein(i) = 0; 
    else
        joints.trans_rein(i) = 1; 
    end
    joints.k0(i) = 6*jnt.e*jnt.iz/jnt.d*(1/1000); % Double Check this is what we want to present here
    joints.Mn{i} = num2str([jnt.Mn, jnt.Mn]*(1/1000));
    joints.Mp{i} = num2str([jnt.Mn, jnt.Mn]*(1/1000)); % Double Check this is what we want to present here
    if strcmp(jnt.class,'a') || strcmp(jnt.class,'b')
        joints.condition(i) = 1;
    else   
        joints.condition(i) = 2;
    end 
    joints.analysis_drift(i) = 0; % LOOK UP EXPLICITY
    joints.analysis_rotation(i) = 0; % LOOK UP EXPLICITY
    joints.recording_drift(i) = 0; % LOOK UP EXPLICITY
    joints.image{i} = 'NA';
    joints.a(i) = jnt.a_hinge;
    joints.b(i) = jnt.b_hinge;
    joints.c(i) = jnt.c_hinge;
    joints.io(i) = jnt.io;
    joints.ls(i) = jnt.ls;
    joints.cp(i) = jnt.cp;
end

% Save Data
writetable(beams,[analysis_dir filesep 'damage_data_beams.csv'])
writetable(columns,[analysis_dir filesep 'damage_data_columns.csv'])
writetable(joints,[analysis_dir filesep 'damage_data_joints.csv'])