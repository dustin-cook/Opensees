function [ ] = main_post_process_opensees( analysis, model, story, node, element, hinge, ground_motion, output_dir )
% Main function that load raw opensees recorder data and transoforms it into something more readily usable 

%% Import Packages
import opensees.post_process.*

%% Load in Analysis data
if analysis.type == 1 % dynamic analysis
    % Ground mottion data
    dirs_ran = fieldnames(ground_motion);
    
    % Load element force data
    element_force_recorders = dlmread([output_dir filesep 'element_force.txt'],' ');
    time_step_vector = element_force_recorders(:,1)';
else % pushover analysis
    dirs_ran = {analysis.pushover_direction};
end

%% Element Forces
% Force component IDs
comp_names = {'P_TH_1','V_TH_1','M_TH_1','M_TH_2'};
if analysis.full_recorders == 1
    if length(dirs_ran) == 1 % 2D
        num_comps = 6;
        comp_keys = [1,2,3,6];
    elseif length(dirs_ran) == 3 % 3D
        num_comps = 12;
        comp_keys = [1,2,6,12];
    end
else
    num_comps = 4;
    comp_keys = [1,2,3,4];
end

%% Loop through elements and save data
if analysis.type == 1 % dynamic analysis
    for i = 1:length(element.id)
        ele_force_TH = element_force_recorders(:,((i-1)*num_comps+2):(i*num_comps+1));
        ele_force_max_abs = max(abs(ele_force_TH));
        ele_force_max = max(ele_force_TH);
        ele_force_min = min(ele_force_TH);
        for j = 1:length(comp_names)
            element_TH.(['ele_' num2str(element.id(i))]).(comp_names{j}) = ele_force_TH(:,comp_keys(j))';
        end
        element.P_grav(i) = ele_force_TH(1,1);
        element.Pmax(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).P_TH_1));
        element.Pmin(i) = min(abs(element_TH.(['ele_' num2str(element.id(i))]).P_TH_1));
        element.Vmax(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).V_TH_1));
        element.Mmax(i) = max(abs([element_TH.(['ele_' num2str(element.id(i))]).M_TH_1,element_TH.(['ele_' num2str(element.id(i))]).M_TH_1]));
    end
end
    
% clear raw opesees data
clear element_force_recorders

%% Load hinge moment and rotation
if analysis.nonlinear ~= 0
    for i = 1:height(hinge)
        rotation_TH = dlmread([output_dir filesep 'hinge_rotation_all.txt'],' ');
        hinge.rotation_TH{i} = rotation_TH(:,i+1)';
        moment_TH = dlmread([output_dir filesep 'hinge_moment_all.txt'],' ');
        hinge.moment_TH{i} = moment_TH(:,i+1)';
    end
end

%% Perform calcs For each direction
for i = 1:length(dirs_ran)
    %% Load and Read Outputs
    if analysis.type == 1 % Dynamic Analysis
        eq.(dirs_ran{i}) = load([ground_motion.(dirs_ran{i}).eq_dir{1} filesep ground_motion.(dirs_ran{i}).eq_name{1}]);
        % Scale EQ (linear interpolation) based on time step used
        eq_length = ground_motion.(dirs_ran{i}).eq_length;
        eq_dt = ground_motion.(dirs_ran{i}).eq_dt;
        eq_timespace = linspace(eq_dt,eq_length*eq_dt,eq_length);
        eq.(dirs_ran{i}) = interp1(eq_timespace,eq.(dirs_ran{i}),time_step_vector);
    end
    
   % EDP response history at each node
   node_disp_raw = dlmread([output_dir filesep ['nodal_disp_' dirs_ran{i} '.txt']],' ')';
   node.(['disp_' dirs_ran{i} '_TH']) = node_disp_raw(2:end,:);
   if analysis.type == 1 % Dynamic Analysis
       node_accel_raw = dlmread([output_dir filesep ['nodal_accel_' dirs_ran{i} '.txt']],' ')';
       node.(['accel_' dirs_ran{i} '_rel_TH']) = node_disp_raw(2:end,:)/386; % Convert to G
       node.(['accel_' dirs_ran{i} '_abs_TH']) = node.(['accel_' dirs_ran{i} '_rel_TH'])/386 + ones(length(node.id),1)*eq.(dirs_ran{i})(1:length(node.(['accel_' dirs_ran{i} '_rel_TH'])));
   elseif analysis.type == 2 % Pushover Analysis
       node_reac_raw = dlmread([output_dir filesep ['nodal_reaction_' dirs_ran{i} '.txt']],' ')';
       node.(['reaction_' dirs_ran{i} '_TH']) = node_disp_raw(2:end,:);
   end
    
    % Max edp's at each node
    node.(['max_disp_' dirs_ran{i}]) = max(abs(node.(['disp_' dirs_ran{i} '_TH'])),[],2);
    if analysis.type == 1 % Dynamic Analysis
        node.(['max_accel_' dirs_ran{i} '_rel']) = max(abs(node.(['accel_' dirs_ran{i} '_rel_TH'])),[],2);
        node.(['max_accel_' dirs_ran{i} '_abs']) = max(abs(node.(['accel_' dirs_ran{i} '_abs_TH'])),[],2);
    end
    
    % EDP Profiles
    [ story.(['max_disp_' dirs_ran{i}]) ] = fn_calc_max_repsonse_profile( node.(['max_disp_' dirs_ran{i}]), story, node, 0 );
    [ story.(['ave_disp_' dirs_ran{i}]) ] = fn_calc_max_repsonse_profile( node.(['max_disp_' dirs_ran{i}]), story, node, 1 );
    if analysis.type == 1 % Dynamic Analysis
        [ story.(['max_accel_' dirs_ran{i}]) ] = fn_calc_max_repsonse_profile( node.(['max_accel_' dirs_ran{i} '_abs']), story, node, 0 );
    end
    [ story.(['max_drift_' dirs_ran{i}]) ] = fn_drift_profile( node.(['disp_' dirs_ran{i} '_TH']), story, node );
    
    % Load Mode shape data and period
    if analysis.run_eigen
        periods = dlmread([output_dir filesep 'period.txt']);
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
end

%% Save Specific Data
save([output_dir filesep 'element_analysis.mat'],'element')
save([output_dir filesep 'node_analysis.mat'],'node')
save([output_dir filesep 'hinge_analysis.mat'],'hinge')
save([output_dir filesep 'story_analysis.mat'],'story')
if analysis.type == 1 % Dynamic Analysis
    save([output_dir filesep 'element_TH.mat'],'element_TH')
    save([output_dir filesep 'gm_data.mat'],'eq','dirs_ran','ground_motion')
end
%% Save All Data
% save([output_dir filesep 'post_process_data'])

end

