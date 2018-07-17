% Post Process Opensees Results
clear
close all
rehash
clc

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 1;
analysis.name = '11DL11LL';

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

element = readtable([output_dir filesep 'element.csv'],'ReadVariableNames',true);
node = readtable([output_dir filesep 'node.csv'],'ReadVariableNames',true);

load([output_dir filesep 'analysis_data.mat']);

% Load ground motion data
[ dirs_ran, ground_motion ] = fn_load_gm_data( ground_motion_seq, gm_table );

% Load element force data
for i = 1:length(element.id)
    if ~strcmp(element.type{i},'wall')
        ele_force = max(abs(dlmread([output_dir filesep ['element_force_' num2str(element.id(i)) '.txt']],' ')));
        if length(dirs_ran) == 1 % 2D
            element.Pmax(i) = max(abs([ele_force(1),ele_force(4)]),[],2);
            element.Vmax(i) = max(abs([ele_force(2),ele_force(5)]),[],2);
            element.Mmax(i) = max(abs([ele_force(3),ele_force(6)]),[],2);
        elseif length(dirs_ran) == 2 % 3D
            element.Pmax(i) = max(abs([ele_force(1),ele_force(7)]),[],2);
            element.Vmax(i) = max(abs([ele_force(2),ele_force(8),ele_force(3),ele_force(9)]),[],2);
            element.Mmax(i) = max(abs([ele_force(4),ele_force(10),ele_force(5),ele_force(11),ele_force(6),ele_force(12)]),[],2);
        end
    else
        test = 5;
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
    
    % Save periods
    model.(['T1_' dirs_ran(i)]) = periods(i);
    
    % Load Mode shape data
    mode_shape_raw = dlmread([output_dir filesep ['mode_shape_' num2str(i) '.txt']]);
    mode_shape_norm = mode_shape_raw(1:2:end)/mode_shape_raw(end-1); % Extract odd rows and normalize by roof
    story.(['mode_shape_' dirs_ran(i)]) = mode_shape_norm';
end

%% Save element table
writetable(element,[output_dir filesep 'element.csv'])

%% Save Data
save([output_dir filesep 'post_process_data'])
