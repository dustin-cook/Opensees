function [ ] = main_post_process_opensees( analysis, model, story, node, element, hinge, ground_motion, opensees_dir )
% Main function that load raw opensees recorder data and transoforms it into something more readily usable 

%% Import Packages
import opensees.post_process.*

%% Load in Analysis data
if analysis.type == 1 % dynamic analysis
    % Ground mottion data
    dirs_ran = fieldnames(ground_motion);
    
    % Load element force data
    if analysis.write_xml
        [ element_force_recorders ] = fn_xml_read([opensees_dir filesep 'element_force.xml']);
    else
        element_force_recorders = dlmread([opensees_dir filesep 'element_force.txt'],' ');
    end
    
    % Define Time Step Vector from element force output
    time_step_vector = element_force_recorders(:,1)';
    
    % Omit Y direction if ran
    dirs_ran = dirs_ran(~strcmp(dirs_ran,'y'));

else % pushover analysis
    % Load element force data
    if analysis.write_xml
        [ element_force_recorders_x ] = fn_xml_read([opensees_dir filesep 'element_force_x.xml']);
        if strcmp(model.dimension,'3D')
            [ element_force_recorders_z ] = fn_xml_read([opensees_dir filesep 'element_force_z.xml']);
        end
    else
        element_force_recorders_x = dlmread([opensees_dir filesep 'element_forcee_x.txt'],' ');
        if strcmp(model.dimension,'3D')
            element_force_recorders_z = dlmread([opensees_dir filesep 'element_forcee_z.txt'],' ');
        end
    end
    
    % Define Direction Ran
    if strcmp(model.dimension,'3D')
        dirs_ran = {'x', 'z'};
    else
        dirs_ran = {'x'};
    end
end

%% Element Forces
% Force component IDs
comp_names = {'P_TH_1','V_TH_1','M_TH_1','M_TH_2'};
num_comps = 4;
comp_keys = [1,2,3,4];

%% Loop through elements and save data
for i = 1:length(element.id)
    % Force Time Histories
    if analysis.type == 1 % dynamic analysis
        ele_force_TH = element_force_recorders(:,((i-1)*num_comps+2):(i*num_comps+1));
    else % pushover analysis
        if strcmp(element.direction{i},'x')
            ele_force_TH = element_force_recorders_x(:,((i-1)*num_comps+2):(i*num_comps+1));
        elseif strcmp(element.direction{i},'z')
            ele_force_TH = element_force_recorders_z(:,((i-1)*num_comps+2):(i*num_comps+1));
        end
    end
    for j = 1:length(comp_names)
        element_TH.(['ele_' num2str(element.id(i))]).(comp_names{j}) = ele_force_TH(:,comp_keys(j))';
    end
    
    % Max Force for each element
    element.P_grav(i) = ele_force_TH(1,1);
    element.Pmax(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).P_TH_1));
    element.Pmin(i) = min(abs(element_TH.(['ele_' num2str(element.id(i))]).P_TH_1));
    element.Vmax(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).V_TH_1));
    element.Mmax(i) = max(abs([element_TH.(['ele_' num2str(element.id(i))]).M_TH_1,element_TH.(['ele_' num2str(element.id(i))]).M_TH_1]));
end
    
% clear raw opesees data
clear element_force_recorders

%% Load hinge moment and rotation
if analysis.nonlinear ~= 0
    if analysis.write_xml
        [ deformation_TH ] = fn_xml_read([opensees_dir filesep 'hinge_deformation_all.xml']);
        [ force_TH ] = fn_xml_read([opensees_dir filesep 'hinge_force_all.xml']);
    else
        deformation_TH = dlmread([opensees_dir filesep 'hinge_deformation_all.txt'],' ');
        force_TH = dlmread([opensees_dir filesep 'hinge_force_all.txt'],' ');
    end
    for i = 1:height(hinge)
        hinge.deformation_TH{i} = deformation_TH(:,2*i-1+1)';
        hinge.shear_TH{i} = -force_TH(:,2*i-1+1)';
        hinge.rotation_TH{i} = deformation_TH(:,2*i+1)';
        hinge.moment_TH{i} = -force_TH(:,2*i+1)'; % I think the forces here are coming in backward, but should triple check
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
        eq_analysis_timespace = time_step_vector;
        eq_analysis.(dirs_ran{i}) = interp1(eq_timespace,eq.(dirs_ran{i}),eq_analysis_timespace);
    end
    
   % EDP response history at each node
   if analysis.write_xml
       [ node_disp_raw ] = fn_xml_read([opensees_dir filesep 'nodal_disp_' dirs_ran{i} '.xml']);
       node_disp_raw = node_disp_raw'; % flip to be node per row
   else
       node_disp_raw = dlmread([opensees_dir filesep 'nodal_disp_' dirs_ran{i} '.txt'],' ')';
   end
   node.(['disp_' dirs_ran{i} '_TH']) = node_disp_raw(2:(height(node)+1),:);
   if analysis.type == 1 % Dynamic Analysis
       if analysis.write_xml
           [ node_accel_raw ] = fn_xml_read([opensees_dir filesep 'nodal_accel_' dirs_ran{i} '.xml']);
           node_accel_raw = node_accel_raw'; % flip to be node per row
       else
           node_accel_raw = dlmread([opensees_dir filesep 'nodal_accel_' dirs_ran{i} '.txt'],' ')';
       end
       node.(['accel_' dirs_ran{i} '_rel_TH']) = node_accel_raw(2:(height(node)+1),:)/386; % Convert to G
       node.(['accel_' dirs_ran{i} '_abs_TH']) = node_accel_raw(2:(height(node)+1),:)/386 + ones(height(node),1)*eq_analysis.(dirs_ran{i});
   elseif analysis.type == 2 % Pushover Analysis
       if analysis.write_xml
           [ node_reac_raw ] = fn_xml_read([opensees_dir filesep 'nodal_reaction_' dirs_ran{i} '.xml']);
           node_reac_raw = node_reac_raw'; % flip to be node per row
       else
           node_reac_raw = dlmread([opensees_dir filesep 'nodal_reaction_' dirs_ran{i} '.txt'],' ')';
       end
       node.(['reaction_' dirs_ran{i} '_TH']) = node_reac_raw(2:end,:);
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
        periods = dlmread([opensees_dir filesep 'period.txt']);
        if strcmp(dirs_ran{i},'x')
            % Save periods
            model.(['T1_' dirs_ran{i}]) = periods(1);
            % Save mode shapes
            if analysis.write_xml
                [ mode_shape_raw ] = fn_xml_read([opensees_dir filesep 'mode_shape_1.xml']);
            else
                mode_shape_raw = dlmread([opensees_dir filesep 'mode_shape_1.txt']);
            end
            mode_shape_norm = mode_shape_raw(1:2:end)/mode_shape_raw(end-1); % Extract odd rows and normalize by roof
            story.(['mode_shape_x']) = mode_shape_norm';
        elseif strcmp(dirs_ran{i},'z')
            % Save periods
            model.(['T1_' dirs_ran{i}]) = periods(2);
            % Save mode shapes
            if analysis.write_xml
                [ mode_shape_raw ] = fn_xml_read([opensees_dir filesep 'mode_shape_2.xml']);
            else
                mode_shape_raw = dlmread([opensees_dir filesep 'mode_shape_2.txt']);
            end
            mode_shape_norm = mode_shape_raw(1:2:end)/mode_shape_raw(end-1); % Extract odd rows and normalize by roof
            story.(['mode_shape_z']) = mode_shape_norm';
        end
    end
end

%% Save Specific Data
save([opensees_dir filesep 'model_analysis.mat'],'model')
save([opensees_dir filesep 'element_analysis.mat'],'element')
save([opensees_dir filesep 'node_analysis.mat'],'node')
save([opensees_dir filesep 'hinge_analysis.mat'],'hinge')
save([opensees_dir filesep 'story_analysis.mat'],'story')
save([opensees_dir filesep 'element_TH.mat'],'element_TH')
if analysis.type == 1 % Dynamic Analysis
    save([opensees_dir filesep 'gm_data.mat'],'eq','dirs_ran','ground_motion','eq_analysis_timespace','eq_analysis')
elseif analysis.type == 2 % Pushover Analysis
    pushover_dir = ['outputs' filesep model.name{1} filesep analysis.proceedure filesep 'pushover'];
    fn_make_directory( pushover_dir )
    save([pushover_dir filesep 'node_analysis.mat'],'node') 
    save([pushover_dir filesep 'element_TH.mat'],'element_TH')
    save([pushover_dir filesep 'hinge_analysis.mat'],'hinge')
    save([pushover_dir filesep 'analysis_options.mat'],'analysis')
end

end

