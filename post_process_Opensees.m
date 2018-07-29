% Post Process Opensees Results
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 6;
analysis.name = 'NL_10DL10LL';

%% Import Packages
import tools.*

%% Load Analysis Data
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
gm_seq_table = readtable(['inputs' filesep 'ground_motion_sequence.csv'],'ReadVariableNames',true);
gm_table = readtable(['inputs' filesep 'ground_motion.csv'],'ReadVariableNames',true);
ground_motion_seq = gm_seq_table(gm_seq_table.id == analysis.gm_id,:);
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

output_dir = ['outputs' filesep model.name{1} filesep analysis.name];

element_table = readtable([output_dir filesep 'element.csv'],'ReadVariableNames',true);
node = readtable([output_dir filesep 'node.csv'],'ReadVariableNames',true);

load([output_dir filesep 'analysis_data.mat']); % Get rid of this system in favor of loading specific data
clear element

% story.nodes_on_slab{1} = story.nodes_on_slab{1}(1:32);

% Filter table to remove rigid elements
element = element_table(ismember(element_table.type,{'beam','column'}) & element_table.ele_id ~= 16 & element_table.ele_id ~= 17,:); % Also removes slab beams (change ele type later)

% Load ground motion data
[ dirs_ran, ground_motion ] = fn_load_gm_data( ground_motion_seq, gm_table );

% Load element force data
for i = 1:length(element.id)
    ele_force_TH = dlmread([output_dir filesep ['element_force_' num2str(element.id(i)) '.txt']],' ');
    ele_force_max_abs = max(abs(ele_force_TH));
    ele_force_max = max(ele_force_TH);
    ele_force_min = min(ele_force_TH);
    if length(dirs_ran) == 1 % 2D
        element_TH.(['ele_' num2str(element.id(i))]).P_TH_1 = ele_force_TH(:,1)';
        element_TH.(['ele_' num2str(element.id(i))]).P_TH_2 = ele_force_TH(:,4)';
        element_TH.(['ele_' num2str(element.id(i))]).V_TH_1 = ele_force_TH(:,2)';
        element_TH.(['ele_' num2str(element.id(i))]).V_TH_2 = ele_force_TH(:,5)';
        element_TH.(['ele_' num2str(element.id(i))]).M_TH_1 = ele_force_TH(:,3)';
        element_TH.(['ele_' num2str(element.id(i))]).M_TH_2 = ele_force_TH(:,6)';
        element.Pmax(i) = max([ele_force_max(1),-ele_force_max(4)],[],2);
        element.Pmin(i) = min([ele_force_min(1),-ele_force_min(4)],[],2);
        element.P_grav(i) = ele_force_TH(1,1);
        element.Vmax(i) = max(abs([ele_force_max_abs(2),ele_force_max_abs(5)]),[],2);
        element.Mmax(i) = max(abs([ele_force_max_abs(3),ele_force_max_abs(6)]),[],2);
    elseif length(dirs_ran) == 3 % 3D
        element_TH.(['ele_' num2str(element.id(i))]).P_TH_1 = ele_force_TH(:,1)';
        element_TH.(['ele_' num2str(element.id(i))]).V_TH_1 = ele_force_TH(:,2)';
        element_TH.(['ele_' num2str(element.id(i))]).M_TH_1 = ele_force_TH(:,6)';
        element_TH.(['ele_' num2str(element.id(i))]).M_TH_2 = ele_force_TH(:,12)';
        element.Pmax(i) = max([ele_force_max(1),-ele_force_max(7)],[],2);
        element.Pmin(i) = min([ele_force_min(1),-ele_force_min(7)],[],2);
        element.P_grav(i) = ele_force_TH(1,1);
        element.Vmax(i) = max(abs([ele_force_max_abs(2),ele_force_max_abs(8),ele_force_max_abs(3),ele_force_max_abs(9)]),[],2);
        element.Mmax(i) = max(abs([ele_force_max_abs(4),ele_force_max_abs(10),ele_force_max_abs(5),ele_force_max_abs(11),ele_force_max_abs(6),ele_force_max_abs(12)]),[],2);
    end
end
    
% Load hinge moment and rotation
if analysis.nonlinear ~= 0
    hinge.rotation = max(abs(dlmread([output_dir filesep 'hinge_rotation_all.txt'],' ')));
end

% Load Period data
periods = dlmread([output_dir filesep 'period.txt']);

% Perform calcs For each direction
for i = 1:length(dirs_ran)
    %% Load and Read Outputs
    eq.(dirs_ran(i)) = load([ground_motion.(dirs_ran(i)).eq_dir{1} filesep ground_motion.(dirs_ran(i)).eq_name{1}]);

   % EDP response history at each node
    node.(['disp_' dirs_ran(i) '_TH']) = dlmread([output_dir filesep ['nodal_disp_' dirs_ran(i) '.txt']],' ')';
    node.(['accel_' dirs_ran(i) '_rel_TH']) = dlmread([output_dir filesep ['nodal_accel_' dirs_ran(i) '.txt']],' ')'/386; % Convert to G
    node.(['accel_' dirs_ran(i) '_abs_TH']) = node.(['accel_' dirs_ran(i) '_rel_TH'])/386 + ones(length(node.id),1)*eq.(dirs_ran(i))';
    
    % Max edp's at each node
    node.(['max_disp_' dirs_ran(i)]) = max(abs(node.(['disp_' dirs_ran(i) '_TH'])),[],2);
    node.(['max_accel_' dirs_ran(i) '_rel']) = max(abs(node.(['accel_' dirs_ran(i) '_rel_TH'])),[],2);
    node.(['max_accel_' dirs_ran(i) '_abs']) = max(abs(node.(['accel_' dirs_ran(i) '_abs_TH'])),[],2);
    
    % EDP Profiles
    [ story.(['max_disp_' dirs_ran(i)]) ] = fn_calc_max_repsonse_profile( node.(['max_disp_' dirs_ran(i)]), story, 0 );
    [ story.(['ave_disp_' dirs_ran(i)]) ] = fn_calc_max_repsonse_profile( node.(['max_disp_' dirs_ran(i)]), story, 1 );
    [ story.(['max_accel_' dirs_ran(i)]) ] = fn_calc_max_repsonse_profile( node.(['max_accel_' dirs_ran(i) '_abs']), story, 0 );
    [ story.(['max_drift_' dirs_ran(i)]) ] = fn_drift_profile( node.(['disp_' dirs_ran(i) '_TH']), story );
    
    % Load Mode shape data
    if strcmp(dirs_ran(i),'x')
        % Save periods
        model.(['T1_' dirs_ran(i)]) = periods(1);
        % Save mode shapes
        mode_shape_raw = dlmread([output_dir filesep ['mode_shape_1.txt']]);
        mode_shape_norm = mode_shape_raw(1:2:end)/mode_shape_raw(end-1); % Extract odd rows and normalize by roof
        story.(['mode_shape_x']) = mode_shape_norm';
    elseif strcmp(dirs_ran(i),'z')
        % Save periods
        model.(['T1_' dirs_ran(i)]) = periods(2);
        % Save mode shapes
        mode_shape_raw = dlmread([output_dir filesep ['mode_shape_2.txt']]);
        mode_shape_norm = mode_shape_raw(1:2:end)/mode_shape_raw(end-1); % Extract odd rows and normalize by roof
        story.(['mode_shape_z']) = mode_shape_norm';
    end
end

%% Save element Data
save([output_dir filesep 'element_TH.mat'],'element_TH')
save([output_dir filesep 'element_analysis.mat'],'element')

%% Save Data
save([output_dir filesep 'post_process_data'])
