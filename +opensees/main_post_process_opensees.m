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

%% Load in Analysis data
if analysis.type == 1 % dynamic analysis
    % Load Hinge Data
    if analysis.nonlinear ~= 0
        hinge_groups = unique(hinge.group);
        for i = 1:length(hinge_groups)
            if analysis.write_xml
                [ hinge_group_deformation_TH ] = fn_xml_read([opensees_dir filesep 'hinge_deformation_group_' num2str(hinge_groups(i)) '.xml']);
                [ hinge_group_force_TH ] = fn_xml_read([opensees_dir filesep 'hinge_force_group_' num2str(hinge_groups(i)) '.xml']);
            else
                hinge_group_deformation_TH = dlmread([opensees_dir filesep 'hinge_deformation_group_' num2str(hinge_groups(i)) '.txt'],' ');
                hinge_group_force_TH = dlmread([opensees_dir filesep 'hinge_force_group_' num2str(hinge_groups(i)) '.txt'],' ');
            end
            hinge_group = hinge(hinge.group == hinge_groups(i),:);
            for h = 1:height(hinge_group)
                hinge_deformation_TH.(['hinge_' num2str(hinge_group.id(h))]) = hinge_group_deformation_TH(:,1+h);
                if strcmp(model.dimension,'3D')
                    hinge_force_TH.(['hinge_' num2str(hinge_group.id(h))]) = hinge_group_force_TH(:,(1+(h-1)*5+1):(1+h*5));
                elseif strcmp(model.dimension,'2D')
                    hinge_force_TH.(['hinge_' num2str(hinge_group.id(h))]) = hinge_group_force_TH(:,(1+(h-1)*2+1):(1+h*2));
                end
            end
        end
    end
else % pushover analysis or Cyclic
    % Hinges
    if analysis.nonlinear ~= 0
        hinge_groups = unique(hinge.group);
        for i = 1:length(hinge_groups)
            if analysis.write_xml
                [ hinge_group_deformation_TH_x ] = fn_xml_read([opensees_dir filesep 'hinge_deformation_x_group_' num2str(hinge_groups(i)) '.xml']);
                [ hinge_group_force_TH_x ] = fn_xml_read([opensees_dir filesep 'hinge_force_x_group_' num2str(hinge_groups(i)) '.xml']);
                if strcmp(model.dimension,'3D')
                    [ hinge_group_deformation_TH_z ] = fn_xml_read([opensees_dir filesep 'hinge_deformation_z_group_' num2str(hinge_groups(i)) '.xml']);
                    [ hinge_group_force_TH_z ] = fn_xml_read([opensees_dir filesep 'hinge_force_z_group_' num2str(hinge_groups(i)) '.xml']);
                end
            else
                hinge_group_deformation_TH_x = dlmread([opensees_dir filesep 'hinge_deformation_x_group_' num2str(hinge_groups(i)) '.txt'],' ');
                hinge_group_force_TH_x = dlmread([opensees_dir filesep 'hinge_force_x_group_' num2str(hinge_groups(i)) '.txt'],' ');
                if strcmp(model.dimension,'3D')
                    hinge_group_deformation_TH_z = dlmread([opensees_dir filesep 'hinge_deformation_z_group_' num2str(hinge_groups(i)) '.txt'],' ');
                    hinge_group_force_TH_z = dlmread([opensees_dir filesep 'hinge_force_z_group_' num2str(hinge_groups(i)) '.txt'],' ');
                end
            end
            hinge_group = hinge(hinge.group == hinge_groups(i),:);
            for h = 1:height(hinge_group)
                hinge_deformation_TH_x.(['hinge_' num2str(hinge_group.id(h))]) = hinge_group_deformation_TH_x(:,1+h);
                hinge_force_TH_x.(['hinge_' num2str(hinge_group.id(h))]) = hinge_group_force_TH_x(:,(1+(h-1)*3+1):(1+h*3));
                if strcmp(model.dimension,'3D')
                    hinge_deformation_TH_z.(['hinge_' num2str(hinge_group.id(h))]) = hinge_group_deformation_TH_z(:,1+h);
                    hinge_force_TH_z.(['hinge_' num2str(hinge_group.id(h))]) = hinge_group_force_TH_z(:,(1+(h-1)*3+1):(1+h*3));
                end
            end
        end
    end
end

% [ joint_force_recorders_all ] = fn_xml_read([opensees_dir filesep 'joint_force_all.xml']);

%% Element Forces
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
        if analysis.write_xml
            ele_force_TH = fn_xml_read([opensees_dir filesep 'element_force_' num2str(i) '.xml']);
        else
            ele_force_TH = dlmread([opensees_dir filesep 'element_force_' num2str(i) '.txt'],' ');
        end
    else % pushover analysis
        if analysis.write_xml
            if strcmp(element.direction{i},'x')
                ele_force_TH = fn_xml_read([opensees_dir filesep 'element_force_x_' num2str(i) '.xml']);
                ele_force_TH_oop = fn_xml_read([opensees_dir filesep 'element_force_z_' num2str(i) '.xml']);
            elseif strcmp(element.direction{i},'z')
                ele_force_TH = fn_xml_read([opensees_dir filesep 'element_force_z_' num2str(i) '.xml']);
                ele_force_TH_oop = fn_xml_read([opensees_dir filesep 'element_force_x_' num2str(i) '.xml']);
            end
        else
            if strcmp(element.direction{i},'x')
                ele_force_TH = dlmread([opensees_dir filesep 'element_force_x_' num2str(i) '.txt'],' ');
                ele_force_TH_oop = dlmread([opensees_dir filesep 'element_force_z_' num2str(i) '.txt'],' ');
            elseif strcmp(element.direction{i},'z')
                ele_force_TH = dlmread([opensees_dir filesep 'element_force_z_' num2str(i) '.txt'],' ');
                ele_force_TH_oop = dlmread([opensees_dir filesep 'element_force_x_' num2str(i) '.txt'],' ');
            end
        end
    end
    for j = 1:length(comp_names)
        if contains(comp_names{j},'oop') && analysis.type == 2 % Pushover out of plane
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
    
    if strcmp(model.dimension,'3D')
        element.Vmax_oop_1(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).V_TH_oop_1));
        element.Vmax_oop_2(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).V_TH_oop_2));
        element.Mmax_oop_1(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).M_TH_oop_1));
        element.Mmax_oop_2(i) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).M_TH_oop_2));
    end
end

%% Save Element Time History
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

%% Load hinge moment and rotation
% Force component IDs
if strcmp(model.dimension,'3D')
    comp_names = {'Vx','Vy','Vz','Mx','Mz'};
    comp_keys = [1,2,3,4,5];
elseif strcmp(model.dimension,'2D')
    comp_names = {'Vx','Vx','Mz'};
    comp_keys = [1,2,3];
end

if analysis.nonlinear ~= 0
    if analysis.type == 1 % dynamic analysis
        for i = 1:height(hinge)
            ele = element(element.id == hinge.element_id(i),:);
            hinge_TH.(['hinge_' num2str(hinge.id(i))]).deformation_TH = hinge_deformation_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip))';
            if strcmp(element.direction(element.id == hinge.element_id(i)),'x') && strcmp(hinge.direction(i),'oop')
                hinge_TH.(['hinge_' num2str(hinge.id(i))]).moment_TH = hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),4)';
                hinge_TH.(['hinge_' num2str(hinge.id(i))]).shear_TH = -hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),3)';
            elseif strcmp(element.direction(element.id == hinge.element_id(i)),'x')
                if strcmp(model.dimension,'3D')
                    if strcmp(ele.type,'beam')
                        hinge_TH.(['hinge_' num2str(hinge.id(i))]).moment_TH = -hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),5)'; % I think the forces here are coming in backward, but should triple check
                        hinge_TH.(['hinge_' num2str(hinge.id(i))]).shear_TH = -hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),2)';
                    else % walls and columns
                        hinge_TH.(['hinge_' num2str(hinge.id(i))]).moment_TH = -hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),5)'; % I think the forces here are coming in backward, but should triple check
                        hinge_TH.(['hinge_' num2str(hinge.id(i))]).shear_TH = hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),1)';
                    end
                elseif strcmp(model.dimension,'2D')
                    if strcmp(ele.type,'beam')
                        hinge_TH.(['hinge_' num2str(hinge.id(i))]).moment_TH = -hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),3)';
                        hinge_TH.(['hinge_' num2str(hinge.id(i))]).shear_TH = -hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),2)';
                    else % walls and columns
                        hinge_TH.(['hinge_' num2str(hinge.id(i))]).moment_TH = -hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),3)';
                        hinge_TH.(['hinge_' num2str(hinge.id(i))]).shear_TH = hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),1)';
                    end
                end
            elseif strcmp(element.direction(element.id == hinge.element_id(i)),'z') && strcmp(hinge.direction(i),'oop')
                hinge_TH.(['hinge_' num2str(hinge.id(i))]).moment_TH = hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),5)';
                hinge_TH.(['hinge_' num2str(hinge.id(i))]).shear_TH = -hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),1)';
            elseif strcmp(element.direction(element.id == hinge.element_id(i)),'z')
                hinge_TH.(['hinge_' num2str(hinge.id(i))]).moment_TH = -hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),4)';
                if strcmp(ele.type,'beam')
                    hinge_TH.(['hinge_' num2str(hinge.id(i))]).shear_TH = hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),3)';
                else
                    hinge_TH.(['hinge_' num2str(hinge.id(i))]).shear_TH = -hinge_force_TH.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),2)';
                end
            end
        end
        clear hinge_deformation_TH
        clear hinge_force_TH
    else % pushover analysis
        for i = 1:height(hinge) % assumes shear in z and bending in x
            ele = element(element.id == hinge.element_id(i),:);
            if (strcmp(ele.direction,'x') && strcmp(hinge.direction(i),'primary')) || (strcmp(ele.direction,'z') && strcmp(hinge.direction(i),'oop'))
                hinge_TH.(['hinge_' num2str(hinge.id(i))]).deformation_TH = hinge_deformation_TH_x.(['hinge_' num2str(hinge.id(i))])(1:(end-clip))';
                hinge_TH.(['hinge_' num2str(hinge.id(i))]).moment_TH = -hinge_force_TH_x.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),3)'; % I think the forces here are coming in backward, but should triple check
                if strcmp(ele.type,'beam')
                    hinge_TH.(['hinge_' num2str(hinge.id(i))]).shear_TH = -hinge_force_TH_x.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),2)';
                else
                    hinge_TH.(['hinge_' num2str(hinge.id(i))]).shear_TH = hinge_force_TH_x.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),1)';
                end
            elseif (strcmp(ele.direction,'z') && strcmp(hinge.direction(i),'primary')) || (strcmp(ele.direction,'x') && strcmp(hinge.direction(i),'oop'))
                hinge_TH.(['hinge_' num2str(hinge.id(i))]).deformation_TH = hinge_deformation_TH_z.(['hinge_' num2str(hinge.id(i))])(1:(end-clip))';
                if strcmp(ele.type,'wall')
                    hinge_TH.(['hinge_' num2str(hinge.id(i))]).moment_TH = -hinge_force_TH_z.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),3)'; % I think the forces here are coming in backward, but should triple check
                else
                    hinge_TH.(['hinge_' num2str(hinge.id(i))]).moment_TH = hinge_force_TH_z.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),3)'; % I think the forces here are coming in backward, but should triple check
                end
                if strcmp(ele.type,'beam')
                    hinge_TH.(['hinge_' num2str(hinge.id(i))]).shear_TH = -hinge_force_TH_z.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),1)';
                else
                    hinge_TH.(['hinge_' num2str(hinge.id(i))]).shear_TH = -hinge_force_TH_z.(['hinge_' num2str(hinge.id(i))])(1:(end-clip),2)';
                end
            end
        end
    end   
end

% clear raw opesees data
clear hinge_deformation_TH
clear hinge_deformation_TH_x
clear hinge_deformation_TH_z
clear hinge_force_TH
clear hinge_force_TH_z
clear hinge_force_TH_x


%% Calculate Nodal Displacments, Accels, Reactions and Eigen values and vectors
if analysis.type == 1 % dynamic analysis
    % Ground motion data
    dirs_ran = fieldnames(ground_motion);
    % Omit Y direction if ran
    dirs_ran = dirs_ran(~strcmp(dirs_ran,'y'));
else
    % Define Direction Ran
    if strcmp(model.dimension,'3D')
        dirs_ran = {'x', 'z'};
    else
        dirs_ran = {'x'};
    end
end

for i = 1:length(dirs_ran)
   % EDP response history at each node
   if analysis.type == 1 % Dynamic Analysis       
       for n = 1:height(node)
           if node.record_disp(n)
               if analysis.write_xml
                   [ node_disp_raw ] = fn_xml_read([opensees_dir filesep 'nodal_disp_' num2str(node.id(n)) '.xml']);
                   node_disp_raw = node_disp_raw'; % flip to be node per row
               else
                   node_disp_raw = dlmread([opensees_dir filesep 'nodal_disp_' num2str(node.id(n)) '.txt'],' ')';
               end
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

           if node.record_accel(n)
               if analysis.write_xml
                   [ node_accel_raw ] = fn_xml_read([opensees_dir filesep 'nodal_accel_' num2str(node.id(n)) '.xml']);
                   node_accel_raw = node_accel_raw'; % flip to be node per row
               else
                   node_accel_raw = dlmread([opensees_dir filesep 'nodal_accel_' num2str(node.id(n)) '.txt'],' ')';
               end
               
              % Scale acceleration EQ (linear interpolation) based on
              % uniform time step points
              eq_analysis = [];
              if ~isfield(eq_analysis,dirs_ran{i})
                   eq.(dirs_ran{i}) = load([ground_motion.(dirs_ran{i}).eq_dir{1} filesep ground_motion.(dirs_ran{i}).eq_name{1}]);
                   eq_length = ground_motion.(dirs_ran{i}).eq_length;
                   eq_dt = ground_motion.(dirs_ran{i}).eq_dt;
                   eq_timespace = linspace(eq_dt,eq_length*eq_dt,eq_length);
                   eq_analysis_timespace = node_accel_raw(1,:);
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
               if analysis.write_xml
                   [ node_disp_raw ] = fn_xml_read([opensees_dir filesep 'nodal_disp_' dirs_ran{i} '_' num2str(node.id(n)) '.xml']);
                   node_disp_raw = node_disp_raw'; % flip to be node per row
               else
                   node_disp_raw = dlmread([opensees_dir filesep 'nodal_disp_' dirs_ran{i} '_' num2str(node.id(n)) '.txt'],' ')';
               end
                   node_TH.(['node_' num2str(node.id(n)) '_TH']).(['disp_' dirs_ran{i} '_TH']) = node_disp_raw(2,1:(end-clip)); 
                   node.(['max_disp_' dirs_ran{i}])(n) = max(abs(node_disp_raw(2,:)));
                   node.(['residual_disp_' dirs_ran{i}])(n) = NaN;
           else
               node_TH.(['node_' num2str(node.id(n)) '_TH']).(['disp_' dirs_ran{i} '_TH']) = [];
               node.(['max_disp_' dirs_ran{i}])(n) = NaN;
               node.(['residual_disp_' dirs_ran{i}])(n) = NaN;
           end       
        end
        
%         % Base shear reactions 
%         [ base_node_reactions ] = fn_xml_read([opensees_dir filesep 'nodal_base_reaction_' dirs_ran{i} '.xml']);
%         analysis.(['base_shear_' dirs_ran{i}]) = abs(sum(base_node_reactions(:,2:end),2))';
   end
   
   % Base shear reactions
   if analysis.write_xml
       [ base_node_reactions ] = fn_xml_read([opensees_dir filesep 'nodal_base_reaction_' dirs_ran{i} '.xml']);
       story_TH.(['base_shear_' dirs_ran{i} '_TH']) = sum(base_node_reactions(1:(end-clip),2:end),2)';
   else
       base_node_reactions = dlmread([opensees_dir filesep 'nodal_base_reaction_' dirs_ran{i} '.txt'],' ')';
       story_TH.(['base_shear_' dirs_ran{i} '_TH']) = sum(base_node_reactions(2:end,1:(end-clip)),1);
   end
   story.(['max_reaction_' dirs_ran{i}])(1) = max(abs(story_TH.(['base_shear_' dirs_ran{i} '_TH'])));
   clear base_node_reactions
   
        
    % EDP Profiles
    [ story.(['max_disp_' dirs_ran{i}]) ] = fn_calc_max_repsonse_profile( node.(['max_disp_' dirs_ran{i}]), story, node, 0 );
    story.(['max_disp_center_' dirs_ran{i}]) = node.(['max_disp_' dirs_ran{i}])(node.center == 1 & node.record_disp == 1 & node.story > 0);
    [ story.(['ave_disp_' dirs_ran{i}]) ] = fn_calc_max_repsonse_profile( node.(['max_disp_' dirs_ran{i}]), story, node, 1 );
    [ story.(['residual_disp_' dirs_ran{i}]) ] = fn_calc_max_repsonse_profile( node.(['residual_disp_' dirs_ran{i}]), story, node, 1 );
    story.(['torsional_factor_' dirs_ran{i}]) = story.(['max_disp_' dirs_ran{i}]) ./ story.(['ave_disp_' dirs_ran{i}]);
    if analysis.type == 1 % Dynamic Analysis
        [ story.(['max_accel_' dirs_ran{i}]) ] = fn_calc_max_repsonse_profile( node.(['max_accel_' dirs_ran{i} '_abs']), story, node, 0 );
        story.(['max_accel_center_' dirs_ran{i}]) = node.(['max_accel_' dirs_ran{i} '_abs'])(node.center == 1 & node.record_accel == 1 & node.story > 0);
    end
    [ story.(['max_drift_' dirs_ran{i}]) ] = fn_drift_profile( node_TH, story, node, dirs_ran{i} );
    
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
    
    %% Collect response info for each element
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

%% Save Specific Data
save([opensees_dir filesep 'model_analysis.mat'],'model')
save([opensees_dir filesep 'element_analysis.mat'],'element')
save([opensees_dir filesep 'joint_analysis.mat'],'joint')
save([opensees_dir filesep 'node_analysis.mat'],'node')
save([opensees_dir filesep 'hinge_analysis.mat'],'hinge')
save([opensees_dir filesep 'story_analysis.mat'],'story')
if analysis.type == 1 % Dynamic Analysis
    save([opensees_dir filesep 'gm_data.mat'],'eq','dirs_ran','ground_motion','eq_analysis_timespace','eq_analysis')
elseif analysis.type == 2 % Pushover Analysis
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

for i = 1:height(hinge)
    hin_TH = hinge_TH.(['hinge_' num2str(hinge.id(i))]);
    save([opensees_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat'],'hin_TH')
    if analysis.type == 2 % Pushover Analysis
        save([pushover_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat'],'hin_TH')
    end
end

save([opensees_dir filesep 'story_TH.mat'],'story_TH')
if analysis.type == 2 % Pushover Analysis
    save([pushover_dir filesep 'story_TH.mat'],'story_TH')
end


end

