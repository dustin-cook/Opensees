function [ ] = main_post_process_opensees( analysis, model, story, node, element, joint, hinge, ground_motion, opensees_dir )
% Main function that load raw opensees recorder data and transoforms it into something more readily usable 

%% Initial Setup
% Import Packages
import opensees.post_process.*

% Make Ouput Directories 
if analysis.type == 2 % Pushover Analysis
    pushover_dir = [analysis.out_dir filesep 'pushover'];
    fn_make_directory( pushover_dir )
elseif analysis.type == 3 % Cyclic Analysis
    cyclic_dir = [analysis.out_dir filesep 'cyclic'];
    fn_make_directory( cyclic_dir )
end

% Define how much of the end to clip
clip = 5;

% Load timespace of analysis
converge_tol_data = load([opensees_dir filesep 'converge_tol_file.txt'],'-ascii');
eq_analysis_timespace = converge_tol_data(:,1);

%% Element Forces
if ~analysis.simple_recorders   
    % Force component IDs
    if strcmp(model.dimension,'3D')
        comp_names = {'P_TH','V_TH_1','V_TH_oop_1','M_TH_oop_1','M_TH_1','V_TH_2','V_TH_oop_2','M_TH_oop_2','M_TH_2'};
        comp_keys = [2,3,4,6,7,9,10,12,13];
    elseif strcmp(model.dimension,'2D')
        comp_names = {'P_TH','V_TH_1','M_TH_1','V_TH_2','M_TH_2'};
        comp_keys = [2,3,4,6,7];
    end

    % Loop through elements and save data
    for i = 1:length(element.id)
        % Force Time Histories
        if analysis.type == 1 % dynamic analysis
            ele_force_TH = fn_xml_read([opensees_dir filesep 'element_force_' num2str(i) '.xml']);
        else % pushover analysis
            if strcmp(element.direction{i},'x')
                ele_force_TH_pos = fn_xml_read([opensees_dir filesep 'element_force_x_' num2str(i) '.xml']);
                ele_force_TH_neg = fn_xml_read([opensees_dir filesep 'element_force_-x_' num2str(i) '.xml']);
                if strcmp(model.dimension,'3D')
                    ele_force_TH_oop_pos = fn_xml_read([opensees_dir filesep 'element_force_z_' num2str(i) '.xml']);
                    ele_force_TH_oop_neg = fn_xml_read([opensees_dir filesep 'element_force_-z_' num2str(i) '.xml']);
                end
            elseif strcmp(element.direction{i},'z')
                ele_force_TH_pos = fn_xml_read([opensees_dir filesep 'element_force_z_' num2str(i) '.xml']);
                ele_force_TH_oop_pos = fn_xml_read([opensees_dir filesep 'element_force_x_' num2str(i) '.xml']);
                ele_force_TH_neg = fn_xml_read([opensees_dir filesep 'element_force_-z_' num2str(i) '.xml']);
                ele_force_TH_oop_neg = fn_xml_read([opensees_dir filesep 'element_force_-x_' num2str(i) '.xml']);
            end
            min_push_length = min(length(ele_force_TH_pos(:,1)),length(ele_force_TH_neg(:,1)));
            ele_force_TH = max(abs(ele_force_TH_pos(1:min_push_length,:)),abs(ele_force_TH_neg(1:min_push_length,:)));
            if strcmp(model.dimension,'3D')
                min_push_length_oop = min(length(ele_force_TH_oop_pos(:,1)),length(ele_force_TH_oop_neg(:,1)));
                ele_force_TH_oop = max(abs(ele_force_TH_oop_pos(1:min_push_length_oop,:)),abs(ele_force_TH_oop_neg(1:min_push_length_oop,:)));
            end
        end
        for j = 1:length(comp_names)
            if strcmp(model.dimension,'3D') && contains(comp_names{j},'oop') && analysis.type == 2 % Pushover out of plane
                element_TH.(['ele_' num2str(element.id(i))]).(comp_names{j}) = ele_force_TH_oop(1:(end-clip),comp_keys(j))';
            else
                element_TH.(['ele_' num2str(element.id(i))]).(comp_names{j}) = ele_force_TH(1:(end-clip),comp_keys(j))';
            end
        end

        % Max Force for each element
        element.P_grav(i) = abs(element_TH.(['ele_' num2str(element.id(i))]).P_TH(1));
        element.Pmax(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).P_TH));
        element.Pmin(i) = min(abs(element_TH.(['ele_' num2str(element.id(i))]).P_TH));
        element.Vmax_1(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).V_TH_1));
        element.Vmax_2(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).V_TH_2));
        element.Mmax_1(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).M_TH_1));
        element.Mmax_2(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).M_TH_2));
        element.Mgrav_1(i) = abs(element_TH.(['ele_' num2str(element.id(i))]).M_TH_1(1));
        element.Mgrav_2(i) = abs(element_TH.(['ele_' num2str(element.id(i))]).M_TH_2(1));

        if strcmp(model.dimension,'3D')
            element.Vmax_oop_1(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).V_TH_oop_1));
            element.Vmax_oop_2(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).V_TH_oop_2));
            element.Mmax_oop_1(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).M_TH_oop_1));
            element.Mmax_oop_2(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).M_TH_oop_2));
        end
    end

    % Save Element Time History
    for i = 1:height(element)
        ele_TH = element_TH.(['ele_' num2str(element.id(i))]);
        save([opensees_dir filesep 'element_TH_' num2str(element.id(i)) '.mat'],'ele_TH')
        if analysis.type == 2 % Pushover Analysis
            save([pushover_dir filesep 'element_TH_' num2str(element.id(i)) '.mat'],'ele_TH')
        end
    end

    % clear raw opesees data
    clear ele_force_TH
    clear element_TH
end

%% Load hinge reactions and deformations
% Load Hinge Data
if analysis.nonlinear ~= 0 && ~isempty(hinge)
    max_story = max(hinge.story);
    for s = 1:max_story
        if exist([opensees_dir filesep 'S' num2str(s) '_hinge_rotation_x' '.xml'],'file')
            [ hinge_deformation_x ] = fn_xml_read([opensees_dir filesep 'S' num2str(s) '_hinge_rotation_x' '.xml']);
            [ hinge_force_x ] = fn_xml_read([opensees_dir filesep 'S' num2str(s) '_hinge_moment_x' '.xml']);
            hinge_ids = hinge(strcmp(hinge.ele_direction,'x') & strcmp(hinge.direction,'primary') & strcmp(hinge.type,'rotational') & hinge.story == s,:);
            for h = 1:height(hinge_ids)
                hinge_TH.(['hinge_' num2str(hinge_ids.id(h))]).deformation_TH = hinge_deformation_x(1:(end-clip),1+h);
                hinge_TH.(['hinge_' num2str(hinge_ids.id(h))]).force_TH = hinge_force_x(1:(end-clip),1+h);
            end
        end
        if strcmp(model.dimension,'3D')
            if exist([opensees_dir filesep 'S' num2str(s) '_hinge_deformation_z' '.xml'],'file')
                [ hinge_deformation_z ] = fn_xml_read([opensees_dir filesep 'S' num2str(s) '_hinge_deformation_z' '.xml']);
                [ hinge_force_z ] = fn_xml_read([opensees_dir filesep 'S' num2str(s) '_hinge_shear_z' '.xml']);
                hinge_ids = hinge(strcmp(hinge.ele_direction,'z') & strcmp(hinge.direction,'primary') & strcmp(hinge.type,'shear') & hinge.story == s,:);
                for h = 1:height(hinge_ids)
                    hinge_TH.(['hinge_' num2str(hinge_ids.id(h))]).deformation_TH = hinge_deformation_z(1:(end-clip),1+h);
                    hinge_TH.(['hinge_' num2str(hinge_ids.id(h))]).force_TH = hinge_force_z(1:(end-clip),1+h);
                end
            end

            if exist([opensees_dir filesep 'S' num2str(s) '_hinge_rotation_z_oop' '.xml'],'file')
                [ hinge_deformation_z_oop ] = fn_xml_read([opensees_dir filesep 'S' num2str(s) '_hinge_rotation_z_oop' '.xml']);
                [ hinge_force_z_oop ] = fn_xml_read([opensees_dir filesep 'S' num2str(s) '_hinge_moment_z_oop' '.xml']);
                hinge_ids = hinge(strcmp(hinge.ele_direction,'z') & strcmp(hinge.direction,'oop') & strcmp(hinge.type,'rotational') & hinge.story == s,:);
                for h = 1:height(hinge_ids)
                    hinge_TH.(['hinge_' num2str(hinge_ids.id(h))]).deformation_TH = hinge_deformation_z_oop(1:(end-clip),1+h);
                    hinge_TH.(['hinge_' num2str(hinge_ids.id(h))]).force_TH = hinge_force_z_oop(1:(end-clip),1+h);
                end
            end

            if exist([opensees_dir filesep 'S' num2str(s) '_hinge_rotation_x_oop' '.xml'],'file')
                [ hinge_deformation_x_oop ] = fn_xml_read([opensees_dir filesep 'S' num2str(s) '_hinge_rotation_x_oop' '.xml']);
                [ hinge_force_x_oop ] = fn_xml_read([opensees_dir filesep 'S' num2str(s) '_hinge_moment_x_oop' '.xml']);
                hinge_ids = hinge(strcmp(hinge.ele_direction,'x') & strcmp(hinge.direction,'oop') & strcmp(hinge.type,'rotational') & hinge.story == s,:);
                for h = 1:height(hinge_ids)
                    hinge_TH.(['hinge_' num2str(hinge_ids.id(h))]).deformation_TH = hinge_deformation_x_oop(1:(end-clip),1+h);
                    hinge_TH.(['hinge_' num2str(hinge_ids.id(h))]).force_TH = hinge_force_x_oop(1:(end-clip),1+h);
                end
            end
        end
    end
    
    % Calculate Dissapated Hinge Engery
    for h = 1:height(hinge)
        if analysis.simple_recorders && strcmp(hinge.direction(h),'oop')
            hinge.total_engergy_ft_lbs(h) = NaN;
        else
            hinge.total_engergy_ft_lbs(h) = (1/12)*trapz(hinge_TH.(['hinge_' num2str(hinge.id(h))]).deformation_TH, hinge_TH.(['hinge_' num2str(hinge.id(h))]).force_TH);
            hinge_TH.(['hinge_' num2str(hinge.id(h))]).energy_ft_lbs = (1/12)*cumtrapz(hinge_TH.(['hinge_' num2str(hinge.id(h))]).deformation_TH, hinge_TH.(['hinge_' num2str(hinge.id(h))]).force_TH);
        end
    end
end

%% Load Joint reactions and deformations
joint_TH = [];
if analysis.nonlinear ~= 0 && analysis.joint_explicit == 1 && ~analysis.simple_recorders % Nonlinear Joints
    % Rotational Hinges x direction - primary
    if exist([opensees_dir filesep 'joint_rotation' '.xml'],'file')
        [ joint_deformation ] = fn_xml_read([opensees_dir filesep 'joint_rotation' '.xml']);
        [ joint_force ] = fn_xml_read([opensees_dir filesep 'joint_moment' '.xml']);
        nonlin_joints = joint(joint.story <= analysis.stories_nonlinear,:);
        for j = 1:height(nonlin_joints)
            joint_TH.(['joint_' num2str(nonlin_joints.id(j))]).deformation_TH = joint_deformation(1:(end-clip),1+j);
            joint_TH.(['joint_' num2str(nonlin_joints.id(j))]).force_TH = joint_force(1:(end-clip),1+j);
        end
    end
end

%% Calculate Nodal Displacments, Accels, Reactions and Eigen values and vectors
if analysis.type == 1 % dynamic analysis
    % Ground motion data
    dirs_ran = fieldnames(ground_motion);
    % Omit Y direction if ran
    dirs_ran = dirs_ran(~strcmp(dirs_ran,'y'));
    % feild names are the directions ran (for compatibility with pushover)
    fld_names = dirs_ran;
elseif analysis.type == 2 || analysis.type == 3
    % Define Direction Ran
    if strcmp(model.dimension,'3D')
        dirs_ran = {'x', '-x', 'z', '-z'};
        fld_names = {'x', 'x_neg', 'z', 'z_neg'};
    else
        dirs_ran = {'x', '-x'};
        fld_names = {'x', 'x_neg'};
    end
end

for i = 1:length(dirs_ran)
   % EDP response history at each node
   if analysis.type == 1 % Dynamic Analysis       
       for n = 1:height(node)
           if node.record_disp(n)
               [ node_disp_raw ] = fn_xml_read([opensees_dir filesep 'nodal_disp_' num2str(node.id(n)) '.xml']);
               node_disp_raw = node_disp_raw'; % flip to be node per row

               if strcmp(dirs_ran{i},'x')
                   node_TH.(['node_' num2str(node.id(n)) '_TH']).(['disp_' dirs_ran{i} '_TH']) = node_disp_raw(2,:); 
                   node.(['max_disp_' dirs_ran{i}])(n) = max(abs(node_disp_raw(2,:)));
                   end_of_motion = 5/ground_motion.(dirs_ran{i}).eq_dt;
                   node.(['residual_disp_' dirs_ran{i}])(n) = abs(mean(node_disp_raw(2,(end-end_of_motion):end)));
               else
                   node_TH.(['node_' num2str(node.id(n)) '_TH']).(['disp_' dirs_ran{i} '_TH']) = node_disp_raw(3,:);  
                   node.(['max_disp_' dirs_ran{i}])(n) = max(abs(node_disp_raw(3,:)));
                   end_of_motion = 5/ground_motion.(dirs_ran{i}).eq_dt;
                   node.(['residual_disp_' dirs_ran{i}])(n) = abs(mean(node_disp_raw(3,(end-end_of_motion):end)));
               end   
           else
               node_TH.(['node_' num2str(node.id(n)) '_TH']).(['disp_' dirs_ran{i} '_TH']) = [];
               node.(['max_disp_' dirs_ran{i}])(n) = NaN;
               node.(['residual_disp_' dirs_ran{i}])(n) = NaN;
           end

           if ~analysis.simple_recorders && node.record_accel(n)
               [ node_accel_raw ] = fn_xml_read([opensees_dir filesep 'nodal_accel_' num2str(node.id(n)) '.xml']);
               node_accel_raw = node_accel_raw'; % flip to be node per row
               
              % Scale acceleration EQ (linear interpolation) based on
              % uniform time step points
              eq_analysis = [];
              if ~isfield(eq_analysis,dirs_ran{i})
                   eq.(dirs_ran{i}) = load([ground_motion.(dirs_ran{i}).eq_dir{1} filesep ground_motion.(dirs_ran{i}).eq_name{1}]);
                   eq_length = ground_motion.(dirs_ran{i}).eq_length;
                   eq_dt = ground_motion.(dirs_ran{i}).eq_dt;
                   eq_timespace = linspace(eq_dt,eq_length*eq_dt,eq_length);
                   if strcmp(dirs_ran{i},'x')
                        node_accel_interp = interp1(eq_analysis_timespace,node_accel_raw(2,:),eq_timespace);
                   else
                       node_accel_interp = interp1(eq_analysis_timespace,node_accel_raw(3,:),eq_timespace);
                   end
                   eq_analysis.(dirs_ran{i}) = eq_timespace*analysis.ground_motion_scale_factor;
              end
           
              % Filter out High Frequency Noise
%               low_freq = 0; % hardcode to no high pass filter
%               [ node_accel_filtered ] = fn_fft_accel_filter( node_accel_interp, eq_dt, eq_timespace, analysis.filter_high_freq, low_freq );
%               [ node_accel_filtered_temp ] = fn_accel_filter( node_accel_raw', 100, 5, 'low');
%               node_accel_filtered = node_accel_filtered_temp';
              
              % compile nodal accels into fields
              node_TH.(['node_' num2str(node.id(n)) '_TH']).(['accel_' dirs_ran{i} '_rel_TH']) = node_accel_interp/386; % Convert to G  
              node_TH.(['node_' num2str(node.id(n)) '_TH']).(['accel_' dirs_ran{i} '_abs_TH']) = node_accel_interp/386 + eq.(dirs_ran{i})';
              node.(['max_accel_' dirs_ran{i} '_rel'])(n) = max(abs(node_accel_interp/386));
              node.(['max_accel_' dirs_ran{i} '_abs'])(n) = max(abs(node_accel_interp/386 + eq.(dirs_ran{i})'));
           else
               node_TH.(['node_' num2str(node.id(n)) '_TH']).(['accel_' dirs_ran{i} '_rel_TH']) = [];
               node_TH.(['node_' num2str(node.id(n)) '_TH']).(['accel_' dirs_ran{i} '_abs_TH']) = [];
               node.(['max_accel_' dirs_ran{i} '_rel'])(n) = NaN;
               node.(['max_accel_' dirs_ran{i} '_abs'])(n) = NaN;
           end
       end
   elseif analysis.type == 2 || analysis.type == 3 % Pushover Analysis or Cyclic
        for n = 1:height(node)
           if node.record_disp(n)
               [ node_disp_raw ] = fn_xml_read([opensees_dir filesep 'nodal_disp_' dirs_ran{i} '_' num2str(node.id(n)) '.xml']);
               node_disp_raw = node_disp_raw'; % flip to be node per row
               node_TH.(['node_' num2str(node.id(n)) '_TH']).(['disp_' fld_names{i} '_TH']) = node_disp_raw(2,1:(end-clip)); 
               node.(['max_disp_' fld_names{i}])(n) = max(abs(node_disp_raw(2,:)));
               node.(['residual_disp_' fld_names{i}])(n) = NaN;
           else
               node_TH.(['node_' num2str(node.id(n)) '_TH']).(['disp_' fld_names{i} '_TH']) = [];
               node.(['max_disp_' fld_names{i}])(n) = NaN;
               node.(['residual_disp_' fld_names{i}])(n) = NaN;
           end       
        end
   end
   
    % EDP Profiles
    [ story.(['max_disp_' fld_names{i}]) ] = fn_calc_max_repsonse_profile( node.(['max_disp_' fld_names{i}]), story, node, 0 );
    story.(['max_disp_center_' fld_names{i}]) = node.(['max_disp_' fld_names{i}])(node.center == 1 & node.record_disp == 1 & node.on_slab == 1);
    [ story.(['ave_disp_' fld_names{i}]) ] = fn_calc_max_repsonse_profile( node.(['max_disp_' fld_names{i}]), story, node, 1 );
    [ story.(['residual_disp_' fld_names{i}]) ] = fn_calc_max_repsonse_profile( node.(['residual_disp_' fld_names{i}]), story, node, 1 );
    story.(['torsional_factor_' fld_names{i}]) = story.(['max_disp_' fld_names{i}]) ./ story.(['ave_disp_' fld_names{i}]);
    if ~analysis.simple_recorders && analysis.type == 1 % Dynamic Analysis
        [ story.(['max_accel_' fld_names{i}]) ] = fn_calc_max_repsonse_profile( node.(['max_accel_' fld_names{i} '_abs']), story, node, 0 );
        story.(['max_accel_center_' fld_names{i}]) = node.(['max_accel_' fld_names{i} '_abs'])(node.center == 1 & node.record_accel == 1 & node.on_slab == 1);
    end
    if analysis.simple_recorders
        [ story.(['max_drift_' fld_names{i}]) ] = fn_drift_profile( node_TH, story, node(node.primary_story == 1,:), fld_names{i} );
    else
        [ story.(['max_drift_' fld_names{i}]) ] = fn_drift_profile( node_TH, story, node, fld_names{i} );
        [ story.(['max_twist_' fld_names{i}]) ] = fn_twist_profile( node_TH, story, node, fld_names{i} );
    end

    % Base shear reactions
    [ base_node_reactions ] = fn_xml_read([opensees_dir filesep 'nodal_base_reaction_' dirs_ran{i} '.xml']);
    story_TH.(['base_shear_' fld_names{i} '_TH']) = sum(base_node_reactions(1:(end-clip),2:end),2)';
    story.(['max_reaction_' fld_names{i}])(1) = max(abs(story_TH.(['base_shear_' fld_names{i} '_TH'])));
    
    % Load Mode shape data and period
    if analysis.run_eigen
        periods = dlmread([opensees_dir filesep 'period.txt']);
        if strcmp(dirs_ran{i},'x')
            % Save periods
            model.(['T1_' dirs_ran{i}]) = periods(1);
            % Save mode shapes
            [ mode_shape_raw ] = fn_xml_read([opensees_dir filesep 'mode_shape_1.xml']);
            mode_shape_norm = mode_shape_raw(1:3:end)/mode_shape_raw(end-2); % Extract every three values and normalize by roof
            story.(['mode_shape_x'])(story.id > 0) = mode_shape_norm';
        elseif strcmp(dirs_ran{i},'z')
            % Save periods
            model.(['T1_' dirs_ran{i}]) = periods(3);
            % Save mode shapes
            [ mode_shape_raw ] = fn_xml_read([opensees_dir filesep 'mode_shape_3.xml']); % fixed as NS first mode = 3rd mode
            mode_shape_norm = mode_shape_raw(2+(1:3:end))/mode_shape_raw(end); % Extract z values and normalize by roof
            story.(['mode_shape_z'])(story.id > 0) = mode_shape_norm';
        end
    end
    
    %% Collect response info for each element
    if strcmp(dirs_ran{i},'x') || strcmp(dirs_ran{i},'z') % dont do for negative pushover directions
        for e = 1:height(element)
            element.(['disp_' dirs_ran{i}])(e,1) = node.(['max_disp_' dirs_ran{i}])(node.id == element.node_2(e)); % Taking the drift at the top of column or wall or the right side of beam
            if element.story(e) == 1
                element.(['drift_' dirs_ran{i}])(e,1) = element.(['disp_' dirs_ran{i}])(e)/story.story_ht(story.id == element.story(e));
            else
                nodes_at_story_below = node(node.story == (element.story(e)-1),:);
                [~, closest_node_idx] = min(sqrt((nodes_at_story_below.x-node.x(node.id == element.node_2(e))).^2 + (nodes_at_story_below.z-node.z(node.id == element.node_2(e))).^2)); % Min pathagorean distance to the closest point
                node_below = nodes_at_story_below(closest_node_idx,:);
                element.(['drift_' dirs_ran{i}])(e,1) = abs(element.(['disp_' dirs_ran{i}])(e)-node_below.max_disp_x)/story.story_ht(story.id == element.story(e));
            end
            if analysis.nonlinear ~= 0 && analysis.type == 1 % nonlinear dynamic analysis
                ele_hinges = hinge(hinge.element_id == element.id(e) & strcmp(hinge.direction,'primary'),:);
                if ~isempty(ele_hinges)
                    if strcmp(ele_hinges.type{1},'rotational') && height(ele_hinges) == 2
                        ele_hinge_TH_1 = hinge_TH.(['hinge_' num2str(ele_hinges.id(1))]);
                        ele_hinge_TH_2 = hinge_TH.(['hinge_' num2str(ele_hinges.id(2))]);
                        rot_1_TH = ele_hinge_TH_1.deformation_TH;
                        rot_2_TH = ele_hinge_TH_2.deformation_TH;
                        element.rot_1(e,1) = max(abs(rot_1_TH));
                        if max(rot_1_TH) > abs(min(rot_1_TH))
                            element.rot_1_dir{e,1} = 'pos';
                        else
                            element.rot_1_dir{e,1} = 'neg';
                        end
                        element.rot_2(e,1) = max(abs(rot_2_TH));
                        if max(rot_2_TH) > abs(min(rot_2_TH))
                            element.rot_2_dir{e,1} = 'pos';
                        else
                            element.rot_2_dir{e,1} = 'neg';
                        end
                    elseif strcmp(ele_hinges.type{1},'shear') && height(ele_hinges) == 1
                        ele_hinge_TH = hinge_TH.(['hinge_' num2str(ele_hinges.id(1))]);
                        element.shear_deform(e,1) = max(ele_hinge_TH.deformation_TH);
                    end
                end
            end
        end
    end
end

%% Save Specific Data
save([opensees_dir filesep 'model_analysis.mat'],'model')
save([opensees_dir filesep 'element_analysis.mat'],'element')
save([opensees_dir filesep 'joint_analysis.mat'],'joint')
save([opensees_dir filesep 'node_analysis.mat'],'node')
save([opensees_dir filesep 'hinge_analysis.mat'],'hinge')
save([opensees_dir filesep 'story_analysis.mat'],'story')
if ~analysis.simple_recorders && analysis.type == 1 % Dynamic Analysis
    save([opensees_dir filesep 'gm_data.mat'],'eq','dirs_ran','ground_motion','eq_analysis_timespace','eq_analysis')
elseif analysis.type == 2 % Pushover Analysis
    save([pushover_dir filesep 'element_analysis.mat'],'element')
    save([pushover_dir filesep 'node_analysis.mat'],'node') 
    save([pushover_dir filesep 'hinge_analysis.mat'],'hinge')
    save([pushover_dir filesep 'analysis_options.mat'],'analysis')
elseif analysis.type == 3 % Cyclic Analysis
    save([cyclic_dir filesep 'hinge_analysis.mat'],'hinge')
    save([cyclic_dir filesep 'analysis_options.mat'],'analysis')
end

% Save Time History Data
for i = 1:height(node)
    nd_TH = node_TH.(['node_' num2str(node.id(i)) '_TH']);
    save([opensees_dir filesep 'node_TH_' num2str(node.id(i)) '.mat'],'nd_TH')
    if analysis.type == 2 % Pushover Analysis
        save([pushover_dir filesep 'node_TH_' num2str(node.id(i)) '.mat'],'nd_TH')
    end
end

if analysis.nonlinear
    for i = 1:height(hinge)
        if isfield(hinge_TH,['hinge_' num2str(hinge.id(i))])
            hin_TH = hinge_TH.(['hinge_' num2str(hinge.id(i))]);
            save([opensees_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat'],'hin_TH')
            if analysis.type == 2 % Pushover Analysis
                save([pushover_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat'],'hin_TH')
            end
        end
    end
    if analysis.joint_explicit == 1
        for i = 1:height(joint)
            if isfield(joint_TH,['joint_' num2str(joint.id(i))])
                jnt_TH = joint_TH.(['joint_' num2str(joint.id(i))]);
                save([opensees_dir filesep 'joint_TH_' num2str(joint.id(i)) '.mat'],'jnt_TH')
                if analysis.type == 2 % Pushover Analysis
                    save([pushover_dir filesep 'joint_TH_' num2str(joint.id(i)) '.mat'],'jnt_TH')
                end
            end
        end
    end
end

if ~analysis.simple_recorders
    save([opensees_dir filesep 'story_TH.mat'],'story_TH')
    if analysis.type == 2 % Pushover Analysis
        save([pushover_dir filesep 'story_TH.mat'],'story_TH')
    end
end


end

