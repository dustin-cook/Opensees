function [ ] = main_post_process_opensees( analysis, model, story, node, element, ground_motion, output_dir )
% Main function that load raw opensees recorder data and transoforms it into something more readily usable 

%% Import Packages
import opensees.post_process.*

%% Load in Analysis data
% Time Step Used in Analysis
time_step_reduce = dlmread([output_dir filesep 'final_time_step_reduction.txt']);
% Load Period data
periods = dlmread([output_dir filesep 'period.txt']);
% Load element force data
element_force_recorders = dlmread([output_dir filesep 'element_force.txt'],' ');
% if exist([output_dir filesep 'element_force_beams_and_columns.txt'],'file')
%     beam_column_force_TH = dlmread([output_dir filesep 'element_force_beams_and_columns.txt'],' ');
% end
% if exist([output_dir filesep 'element_force_walls.txt'],'file')
%     wall_force_TH = dlmread([output_dir filesep 'element_force_walls.txt'],' ');
% end

% Ground mottion data
dirs_ran = fieldnames(ground_motion);

% Element Forces
for i = 1:length(element.id)
    if length(dirs_ran) == 1 % 2D
        ele_force_TH = element_force_recorders(:,((i-1)*6+1):(i*6));
        ele_force_max_abs = max(abs(ele_force_TH));
        ele_force_max = max(ele_force_TH);
        ele_force_min = min(ele_force_TH);
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
        ele_force_TH = element_force_recorders(:,((i-1)*12+1):(i*12));
        ele_force_max_abs = max(abs(ele_force_TH));
        ele_force_max = max(ele_force_TH);
        ele_force_min = min(ele_force_TH);
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
    
% clear raw opesees data
clear element_force_recorders

% Load hinge moment and rotation
if analysis.nonlinear ~= 0
    hinge.rotation = max(abs(dlmread([output_dir filesep 'hinge_rotation_all.txt'],' ')));
end

% Perform calcs For each direction
for i = 1:length(dirs_ran)
    %% Load and Read Outputs
    eq.(dirs_ran{i}) = load([ground_motion.(dirs_ran{i}).eq_dir{1} filesep ground_motion.(dirs_ran{i}).eq_name{1}]);
    % Scale EQ (linear interpolation) based on time step used
    eq_length = ground_motion.(dirs_ran{i}).eq_length;
    eq_dt = ground_motion.(dirs_ran{i}).eq_dt;
    eq_timespace = linspace(eq_dt,eq_length*eq_dt,eq_length);
    analysis_timespace = linspace(eq_dt/time_step_reduce,eq_length*eq_dt,eq_length*(time_step_reduce));
    eq.(dirs_ran{i}) = interp1(eq_timespace,eq.(dirs_ran{i}),analysis_timespace);

   % EDP response history at each node
    node.(['disp_' dirs_ran{i} '_TH']) = dlmread([output_dir filesep ['nodal_disp_' dirs_ran{i} '.txt']],' ')';
    node.(['accel_' dirs_ran{i} '_rel_TH']) = dlmread([output_dir filesep ['nodal_accel_' dirs_ran{i} '.txt']],' ')'/386; % Convert to G
    node.(['accel_' dirs_ran{i} '_abs_TH']) = node.(['accel_' dirs_ran{i} '_rel_TH'])/386 + ones(length(node.id),1)*eq.(dirs_ran{i});
    
    % Max edp's at each node
    node.(['max_disp_' dirs_ran{i}]) = max(abs(node.(['disp_' dirs_ran{i} '_TH'])),[],2);
    node.(['max_accel_' dirs_ran{i} '_rel']) = max(abs(node.(['accel_' dirs_ran{i} '_rel_TH'])),[],2);
    node.(['max_accel_' dirs_ran{i} '_abs']) = max(abs(node.(['accel_' dirs_ran{i} '_abs_TH'])),[],2);
    
    % EDP Profiles
    [ story.(['max_disp_' dirs_ran{i}]) ] = fn_calc_max_repsonse_profile( node.(['max_disp_' dirs_ran{i}]), story, node, 0 );
    [ story.(['ave_disp_' dirs_ran{i}]) ] = fn_calc_max_repsonse_profile( node.(['max_disp_' dirs_ran{i}]), story, node, 1 );
    [ story.(['max_accel_' dirs_ran{i}]) ] = fn_calc_max_repsonse_profile( node.(['max_accel_' dirs_ran{i} '_abs']), story, node, 0 );
    [ story.(['max_drift_' dirs_ran{i}]) ] = fn_drift_profile( node.(['disp_' dirs_ran{i} '_TH']), story, node );
    
    % Load Mode shape data
    if strcmp(dirs_ran{i},'x')
        % Save periods
        model.(['T1_' dirs_ran{i}]) = periods(1);
        % Save mode shapes
        mode_shape_raw = dlmread([output_dir filesep ['mode_shape_1.txt']]);
        mode_shape_norm = mode_shape_raw(1:2:end)/mode_shape_raw(end-1); % Extract odd rows and normalize by roof
        story.(['mode_shape_x']) = mode_shape_norm';
    elseif strcmp(dirs_ran{i},'z')
        % Save periods
        model.(['T1_' dirs_ran{i}]) = periods(2);
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

end

