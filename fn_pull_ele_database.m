function [] = fn_pull_ele_database(model, analysis, ele_props_table)
% Script to pull info from ATC 134 run and format into element collection
% database

%% Define data directories
analysis_dir = ['outputs' filesep model.name{1} filesep analysis.proceedure '_' analysis.id filesep 'asce_41_data'];
opensees_dir = ['outputs' filesep model.name{1} filesep analysis.proceedure '_' analysis.id filesep 'opensees_data'];
write_dir = ['outputs' filesep model.name{1} filesep analysis.proceedure '_' analysis.id filesep 'element_data'];

%% Import Packages
import asce_41.*

%% Load Inputs
ele_inputs = readtable('element_collection_key_story_1.csv','ReadVariableNames',true);
load([analysis_dir filesep 'joint_analysis.mat'])
load([analysis_dir filesep 'element_analysis.mat'])
load([analysis_dir filesep 'story_analysis.mat'])
load([analysis_dir filesep 'hinge_analysis.mat'])
load(['ground_motions' filesep 'ICSB_recordings' filesep 'recorded_edp_profile.mat'])

%% Build element table
prime_hinges = hinge(strcmp(hinge.direction,'primary'),:);
for i = 1:height(prime_hinges)
    element2write.id(i,:) = i;
    element2write.element_id(i,:) = prime_hinges.element_id(i);
    element2write.hinge_side(i,:) = prime_hinges.ele_side(i);
    element2write.element_type{i,:} = element.type{element.id == prime_hinges.element_id(i)};
    element2write.member_id{i,:} = [num2str(prime_hinges.element_id(i)) '_' num2str(prime_hinges.ele_side(i))];
    if prime_hinges.ele_side(i) == 1
        jnt_filter = joint.column_high == prime_hinges.element_id(i) | joint.beam_right == prime_hinges.element_id(i);
    elseif prime_hinges.ele_side(i) == 2
        jnt_filter = joint.column_low == prime_hinges.element_id(i) | joint.beam_left == prime_hinges.element_id(i);
    end
    if any(jnt_filter)
        element2write.joint_id(i,:) = joint.id(jnt_filter);
    else
        element2write.joint_id(i,:) = 1000 + prime_hinges.element_id(i);
    end
    element2write.system_type(i,:) = 1; % Assume Nonconforming for now
    element2write.story(i,:) = prime_hinges.story(i);
end
% Add joints to the mix
c = i;
for j = 1:height(joint)
%     col_low = element(element.id == joint.column_low(j),:);
%     if col_low.story == 1
        c = c + 1;
        element2write.id(c,:) = c;
        element2write.element_id(c,:) = joint.id(j);
        element2write.hinge_side(c,:) = 0;
        element2write.element_type{c,:} = 'joint';
        element2write.member_id{c,:} = ['jnt-' num2str(joint.id(j));];
        element2write.joint_id(c,:) = joint.id(j);
        element2write.system_type(c,:) = 1; % Assume Nonconforming for now
        element2write.story(c,:) = element.story(element.id == joint.column_low(j));
%     end
end
% Add column base joints
base_columns = element(strcmp(element.type,'column') & element.story == 1,:);
b = c;
for col = 1:height(base_columns)
    b = b + 1;
    element2write.id(b,:) = b;
    element2write.element_id(b,:) = 1000+base_columns.id(col);
    element2write.hinge_side(b,:) = 0;
    element2write.element_type{b,:} = 'joint';
    element2write.member_id{b,:} = ['jnt-base' num2str(1000+base_columns.id(col));];
    element2write.joint_id(b,:) = 1000+base_columns.id(col);
    element2write.system_type(b,:) = 1; % Assume Nonconforming for now
    element2write.story(b,:) = 0;
end
ele_inputs = struct2table(element2write);

%% Formulate Beam Table
beams = ele_inputs(strcmp(ele_inputs.element_type,'beam'),:);
for i = 1:height(beams)
    hin_side = num2str(beams.hinge_side(i));
    ele = element(element.id == beams.element_id(i),:);
    hin = hinge(hinge.element_id == ele.id & strcmp(hinge.direction,'primary') & hinge.ele_side == str2double(hin_side),:);
    load([opensees_dir filesep 'hinge_TH_' num2str(hin.id)])
    if ~isempty(ele)
        this_story = story(story.id == ele.story,:);
        ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
        beams.ld_req(i) = ele.ld_req;
        beams.ld_avail{i} = ele.ld_avail;
        if ele.ld_req < ele.ld_avail
            beams.development(i) = 1;
        else
            beams.development(i) = 0;
        end
        beams.rho_ratio(i) = (ele_props.row - ele_props.row_prime)/ele.row_bal;
        beams.shear_stress_ratio(i) = ele.(['Vmax_' hin_side])/(ele_props.w*ele_props.d_eff*sqrt(ele_props.fc_e));
        beams.shear_demand_ratio(i) = ele.(['Vmax_' hin_side])/ele.(['Vn_' hin_side]);
        if strcmp(ele.trans_rein_check,'NC')
            beams.trans_rein(i) = 0;
        elseif strcmp(ele.trans_rein_check,'C')
            beams.trans_rein(i) = 1;
        end
        beams.s(i) = ele_props.(['S_' hin_side]);
        beams.hinge_splices(i) = 0;
        beams.K0(i) = (1/1000)*6*ele_props.e*ele_props.iz/ele.length;
        if abs(max(hin_TH.deformation_TH)) >= abs(min(hin_TH.deformation_TH))% Assumes beam strength is same on each side of the beam
            beams.Mn(i) = (1/1000)*ele.(['Mn_pos_' hin_side]); % Positive Bendiong
            beams.Mp(i) = (1/1000)*ele.(['Mp_pos_' hin_side]);
        elseif abs(max(hin_TH.deformation_TH)) < abs(min(hin_TH.deformation_TH))
            beams.Mn(i) = (1/1000)*ele.(['Mn_neg_' hin_side]); % Negative Bending
            beams.Mp(i) = (1/1000)*ele.(['Mp_neg_' hin_side]);
        end
        if ele.pass_aci_dev_length == 0
            beams.condition(i) = 3;
        elseif strcmp(ele.(['critical_mode_' hin_side]),'shear')
            beams.condition(i) = 2;
        else   
            beams.condition(i) = 1;
        end
        beams.length(i) = ele.length;
        beams.max_element_rotation_or_drift(i) = max(abs(hin_TH.deformation_TH));
        beams.max_story_drift_analysis(i) = this_story.(['max_drift_' ele.direction{1}]); % Need to fix elment story drift read such that it does not come in as NaN
        beams.max_story_drift_record(i) = record_edp.max_disp.(ele.direction{1})(ele.story+1);
        beams.failure_mech_analysis{i} = 'Rotational Yeilding'; % Need to dynamically check this
        beams.failure_mech_recorded{i} = 'No Failure'; % hard coded from ICSB results
        beams.damage_image{i} = 'NA';
        beams.a(i) = ele.(['a_hinge_' hin_side]);
        beams.b(i) = ele.(['b_hinge_' hin_side]);
        beams.c(i) = ele.(['c_hinge_' hin_side]);
        beams.io(i) = ele.(['io_' hin_side]);
        beams.ls(i) = ele.(['ls_' hin_side]);
        beams.cp(i) = ele.(['cp_' hin_side]);
    end
end

% Formulate Column Table
columns = ele_inputs(strcmp(ele_inputs.element_type,'column'),:);
for i = 1:height(columns)
    hin_side = num2str(columns.hinge_side(i));
    ele = element(element.id == columns.element_id(i),:);
    hin = hinge(hinge.element_id == ele.id & strcmp(hinge.direction,'primary') & hinge.ele_side == str2double(hin_side),:);
    load([opensees_dir filesep 'hinge_TH_' num2str(hin.id)])
    if ~isempty(ele)
        this_story = story(story.id == ele.story,:);
        ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
        columns.fc_e(i) = (1/1000)*ele_props.fc_e;
        columns.fy_e(i) = (1/1000)*ele_props.fy_e;
        columns.fyt_e(i) = (1/1000)*ele_props.fy_e; % Assume they are the same
        columns.max_expected_axial_load_ratio(i) = ele.Pmax/(ele_props.a*ele_props.fc_e);% Use expected to related to testing data
        columns.max_nominal_axial_load_ratio(i) = ele.Pmax/(ele_props.a*ele_props.fc_n);% Use nominal capacity since axial is force controlled
        columns.gravity_axial_load_ratio(i) = ele.P_grav/(ele_props.a*ele_props.fc_e);
        columns.shear_flexure_yield_ratio(i) = ele.(['vye_' hin_side])/ele.(['V0_' hin_side]);
        columns.shear_demand_ratio(i) = ele.(['Vmax_' hin_side])/ele.(['V0_' hin_side]);
        columns.k(i) = ele.(['Vn_' hin_side])/ele.(['V0_' hin_side]);
        min_d_b = min(str2double(strsplit(strrep(strrep(ele_props.d_b{1},']',''),'[',''))));
        columns.shear_length_ratio(i) = ele_props.a/ele_props.d_eff;
        columns.rho_t(i) = ele.(['rho_t_' hin_side]);
        columns.rho_l(i) = ele.rho_l;
        columns.s(i) = ele_props.(['S_' hin_side]);
        columns.s_db(i) = ele_props.(['S_' hin_side])/min_d_b;
        columns.hoops_conform(i) = 0; % Even though you can assume rein is anchored into the core based on 135 deg hooks from the plans I am going to assume idadequate here to trigger shear behavior in the bottom story
        columns.hinge_splices(i) = 0; % Plans say that splices are in center of column
        columns.K0(i) = 6*ele_props.e*ele_props.iz/ele.length*(1/1000);
        columns.Mn(i) = ele.(['Mn_pos_' hin_side])*(1/1000); % Assuming columns are the same in both directions (and have the same capacity top and bottom)
        columns.Mp(i) = ele.(['Mp_pos_' hin_side])*(1/1000);
        if ele.pass_aci_dev_length == 0
            columns.condition(i) = 2;
        else   
            columns.condition(i) = 1;
        end 
        columns.length(i) = ele.length;
        columns.yield_rotation(i) = columns.Mn(i)/columns.K0(i);
        columns.max_element_hinge_rotation(i) = max(abs(hin_TH.deformation_TH));
        if columns.max_element_hinge_rotation(i) < columns.yield_rotation(i)/11
            columns.max_element_rotation_or_drift(i) = columns.max_element_hinge_rotation(i)*11;
        else
            columns.max_element_rotation_or_drift(i) = columns.max_element_hinge_rotation(i) + (10/11)*columns.yield_rotation(i);
        end
        columns.max_story_drift_analysis(i) = this_story.(['max_drift_' ele.direction{1}]);
        columns.max_story_drift_record(i) = record_edp.max_disp.(ele.direction{1})(ele.story+1)/story.story_ht(story.id == 1); % Update this to read from story
        columns.failure_mech_analysis{i} = 'ADD FAILURE MECH'; % Need to dynamically check this
        columns.failure_mech_recorded{i} = 'ADD FAILURE MECH'; % need to hard cod from ICSB documenta
        columns.damage_image{i} = ['Fig - ', columns.member_id{i}];
        columns.a(i) = ele.(['a_hinge_' hin_side]);
        columns.b(i) = ele.(['b_hinge_' hin_side]);
        columns.c(i) = ele.(['c_hinge_' hin_side]);
        columns.io(i) = ele.(['io_' hin_side]);
        columns.ls(i) = ele.(['ls_' hin_side]);
        columns.cp(i) = ele.(['cp_' hin_side]);
    end
end

% Formulate Joint Table
joints = ele_inputs(strcmp(ele_inputs.element_type,'joint'),:);
joints.joint_id = [];
for i = 1:height(joints)
    jnt = joint(joint.id == joints.element_id(i),:);
    if ~isempty(jnt) % this is a real joint
        joints.axial_load_ratio(i) = jnt.Pmax/(jnt.a*jnt.fc_e);
        joints.shear_ratio(i) = jnt.Vmax/jnt.Vj;
        joints.scwb(i) = jnt.col_bm_ratio;
        if strcmp(jnt.trans_rien,'NC')
            joints.trans_rein(i) = 0; 
        else
            joints.trans_rein(i) = 1; 
        end
        [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', jnt.Mn, jnt.Mn, inf, inf, jnt.h, jnt.e, jnt.iz, jnt.a_hinge, jnt.b_hinge, jnt.c_hinge, NaN, 0.04, 'shear' );
        joints.K0(i) = (1/1000)*1000*moment_vec_pos(1)/rot_vec_pos(1);
        joints.Mn(i) = moment_vec_pos(1)*(1/1000);
        joints.Mp(i) = moment_vec_pos(2)*(1/1000);
        if strcmp(jnt.class,'a') || strcmp(jnt.class,'b')
            joints.condition(i) = 1;
        else   
            joints.condition(i) = 2;
        end 
        joints.max_element_rotation_or_drift(i) = 0; % Not modeling joint nonlinearity
        joints.max_story_drift_analysis(i) = ele.(['drift_' 'x']); % Hardcode as X for now, all joints are in x direction
        joints.max_story_drift_record(i) = NaN; % Only worry about it for the columns
        joints.failure_mech_analysis{i} = 'No Failure'; % Need to dynamically check this
        joints.failure_mech_recorded{i} = 'No Failure'; % hard coded from ICSB documentation
        joints.damage_image{i} = 'NA';
        joints.a(i) = jnt.a_hinge;
        joints.b(i) = jnt.b_hinge;
        joints.c(i) = jnt.c_hinge;
        joints.io(i) = jnt.io;
        joints.ls(i) = jnt.ls;
        joints.cp(i) = jnt.cp;
    else % This is a fake joint (ie ground joint)
        ele = element(element.id == (joints.element_id(i)-1000),:);
        ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
        joints.axial_load_ratio(i) = ele.Pmax/(ele_props.a*ele_props.fc_e);
        joints.shear_ratio(i) = 0;
        joints.scwb(i) = 0;
        joints.trans_rein(i) = 1; % Not really a joint so its conforming
        joints.K0(i) = inf;
        joints.Mn(i) = inf;
        joints.Mp(i) = inf;
        joints.condition(i) = 1; % Dont think this matters
        joints.max_element_rotation_or_drift(i) = 0; 
        joints.max_story_drift_analysis(i) = 0; 
        joints.max_story_drift_record(i) = NaN; 
        joints.failure_mech_analysis{i} = 'No Failure'; 
        joints.failure_mech_recorded{i} = 'No Failure';
        joints.damage_image{i} = 'NA';
        joints.a(i) = 0;
        joints.b(i) = 0;
        joints.c(i) = 0;
        joints.io(i) = 0;
        joints.ls(i) = 0;
        joints.cp(i) = 0;
    end
end

% Formulate Wall Table
walls = ele_inputs(strcmp(ele_inputs.element_type,'wall'),:);
walls.joint_id = [];
for i = 1:height(walls)
    ele = element(element.id == walls.element_id(i),:);
    hin = hinge(hinge.element_id == ele.id & strcmp(hinge.direction,'primary') & hinge.ele_side == 1,:);
    load([opensees_dir filesep 'hinge_TH_' num2str(hin.id)])
    if ~isempty(ele)
        this_story = story(story.id == ele.story,:);
        ele_props = ele_props_table(ele_props_table.id == ele.ele_id,:);
        As = sum(str2double(strsplit(strrep(strrep(ele_props.As{1},']',''),'[',''))))/2; % Assume 1/2 tensiona and 1/2 compression just like in wall hinge selection
        As_prime = As;
        walls.axial_load_ratio(i) = ((As-As_prime)*ele_props.fy_e+ele.Pmax)/(ele_props.w*ele_props.d_eff*ele_props.fc_e);
        if strcmp(ele.critical_mode_1,'shear')
            walls.condition(i) = 2;
        elseif strcmp(ele.critical_mode_1,'flexure') 
            walls.condition(i) = 1;
        end 
        walls.length(i) = ele.length;
        walls.shear_stress_ratio(i) = ele.Vmax_1/(ele_props.w*ele_props.d_eff*sqrt(ele_props.fc_e));
        walls.confined_boundary(i) = 0; % Assume they are unconfined (doens't matter for shear controlled walls)
        [ force_vec, disp_vec ] = fn_define_backbone_shear( ele.Vn_1, ele.length, ele_props.g, ele_props.av, ele.c_hinge_1, ele.d_hinge_1, ele.e_hinge_1, ele.f_hinge_1, ele.g_hinge_1 ); % Assumes walls are shear controlled, probably should check explicitly
        walls.K0(i) = (1/1000)*force_vec(1)/disp_vec(1);
        walls.Mn(i) = force_vec(2)*(1/1000);
        walls.Mp(i) = force_vec(3)*(1/1000);
        walls.max_element_rotation_or_drift(i) = max(abs(hin_TH.deformation_TH))/ele.length;
        walls.max_story_drift_analysis(i) = max(abs(hin_TH.deformation_TH))/ele.length;
        walls.max_story_drift_record(i) = record_edp.max_disp.(ele.direction{1})(ele.story+1);
        walls.failure_mech_analysis{i} = 'No Failure'; % Need to dynamically check this
        walls.failure_mech_recorded{i} = 'No Failure'; % hard coded from ICSB documentation
        walls.damage_image{i} = 'NA';
        walls.a(i) = ele.a_hinge_1;
        walls.b(i) = ele.b_hinge_1;
        walls.c(i) = ele.c_hinge_1;
        walls.d(i) = ele.d_hinge_1;
        walls.e(i) = ele.e_hinge_1;
        walls.f(i) = ele.f_hinge_1;
        walls.g(i) = ele.g_hinge_1;
        walls.io(i) = ele.io_1;
        walls.ls(i) = ele.ls_1;
        walls.cp(i) = ele.cp_1;
    end
end

% Save Data
if ~exist(write_dir,'dir')
    mkdir(write_dir)
end
if ~isempty(beams)
    writetable(beams,[write_dir filesep 'raw_beams_database.csv'])
end
if ~isempty(columns)
    writetable(columns,[write_dir filesep 'raw_columns_database.csv'])
end
if ~isempty(joints)
    writetable(joints,[write_dir filesep 'raw_joints_database.csv'])
end
if ~isempty(walls)
    writetable(walls,[write_dir filesep 'raw_walls_database.csv'])
end

end