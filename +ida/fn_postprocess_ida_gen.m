function [ ] = fn_postprocess_ida_gen( analysis, model, node, element, ground_motion, opensees_dir, summary_dir )
% Main function that load raw opensees recorder data and transoforms it into something more readily usable 

%% Initial Setup
% Import Packages
% import opensees.post_process.*

% Define how much of the end to clip
clip = 5;

% Load timespace of analysis
converge_tol_data = load([opensees_dir filesep 'converge_tol_file.txt'],'-ascii');
eq_analysis_timespace = converge_tol_data(:,1);

% Load summary data
load([summary_dir filesep 'summary_results.mat'],'summary')

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
for i = 1:length(element)
    % Force Time Histories
    ele_force_TH = fn_xml_read([opensees_dir filesep 'element_force_' num2str(element(i)) '.xml']);
    for j = 1:length(comp_names)
        element_TH.(['ele_' num2str(element(i))]).(comp_names{j}) = ele_force_TH(1:(end-clip),comp_keys(j))';
    end

    % Max Force for each element
    summary.element.id(i,1) = element(i);
    summary.element.P_grav(i,1) = abs(element_TH.(['ele_' num2str(element(i))]).P_TH(1));
    summary.element.Pmax(i,1) = max(abs(element_TH.(['ele_' num2str(element(i))]).P_TH));
    summary.element.Pmin(i,1) = min(abs(element_TH.(['ele_' num2str(element(i))]).P_TH));
    summary.element.Vmax_1(i,1) = max(abs(element_TH.(['ele_' num2str(element(i))]).V_TH_1));
    summary.element.Vmax_2(i,1) = max(abs(element_TH.(['ele_' num2str(element(i))]).V_TH_2));
    summary.element.Mmax_1(i,1) = max(abs(element_TH.(['ele_' num2str(element(i))]).M_TH_1));
    summary.element.Mmax_2(i,1) = max(abs(element_TH.(['ele_' num2str(element(i))]).M_TH_2));
    summary.element.Mgrav_1(i,1) = abs(element_TH.(['ele_' num2str(element(i))]).M_TH_1(1));
    summary.element.Mgrav_2(i,1) = abs(element_TH.(['ele_' num2str(element(i))]).M_TH_2(1));

    if strcmp(model.dimension,'3D')
        summary.element.Vmax_oop_1(i,1) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).V_TH_oop_1));
        summary.element.Vmax_oop_2(i,1) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).V_TH_oop_2));
        summary.element.Mmax_oop_1(i,1) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).M_TH_oop_1));
        summary.element.Mmax_oop_2(i,1) = max(abs(element_TH.(['ele_' num2str(element.id(i))]).M_TH_oop_2));
    end
end

summary.element = struct2table(summary.element);

% clear raw opesees data
clear ele_force_TH
clear element_TH

%% Calculate Nodal Displacments, Accels, Reactions and Eigen values and vectors
% Ground motion data
dirs_ran = fieldnames(ground_motion);
dir_ids = [2,3]; % define output columns to pull from
dirs_ran = dirs_ran(~strcmp(dirs_ran,'y')); % Omit Y direction if ran

for i = 1:length(dirs_ran)
   % Load earthquake information
   eq.(dirs_ran{i}) = load([ground_motion.(dirs_ran{i}).eq_dir{1} filesep ground_motion.(dirs_ran{i}).eq_name{1}]);
   eq_length = ground_motion.(dirs_ran{i}).eq_length;
   eq_dt = ground_motion.(dirs_ran{i}).eq_dt;
   eq_timespace = linspace(eq_dt,eq_length*eq_dt,eq_length);
   
   % EDP response history at each node  
   for n = 1:length(node.id)
       [ node_disp_raw ] = fn_xml_read([opensees_dir filesep 'nodal_disp_' num2str(node.id(n)) '.xml']);
       [ node_accel_raw ] = fn_xml_read([opensees_dir filesep 'nodal_accel_' num2str(node.id(n)) '.xml']);
       [ node_reactions ] = fn_xml_read([opensees_dir filesep 'nodal_reaction_' num2str(node.id(n)) '.xml']);

       % Change interpolate accel to original time step so it can be
       % combined with ground motion data
       node_accel_interp = interp1(eq_analysis_timespace,node_accel_raw(:,dir_ids(i)),eq_timespace);
    
       % Define number of points in the last 5 seconds of the ground motion
       end_of_motion = 5/ground_motion.(dirs_ran{i}).eq_dt;
       
       % Find peak nodal values
       summary.node.id(n,1) = node.id(n);
       summary.node.(['max_disp_' dirs_ran{i}])(n,1) = max(abs(node_disp_raw(:,dir_ids(i))));
       summary.node.(['max_reaction_' dirs_ran{i}])(n,1) = max(abs(node_reactions(:,dir_ids(i))));
       try
       summary.node.(['residual_disp_' dirs_ran{i}])(n,1) = abs(mean(node_disp_raw((end-end_of_motion):end,dir_ids(i))));
       catch
           test = 5;
       end
       summary.node.(['max_accel_' dirs_ran{i} '_rel'])(n,1) = max(abs(node_accel_interp/analysis.g_unit));
       summary.node.(['max_accel_' dirs_ran{i} '_abs'])(n,1) = max(abs(node_accel_interp'/analysis.g_unit + analysis.ground_motion_scale_factor*eq.(dirs_ran{i})));
   end
end

summary.node = struct2table(summary.node);

%% Save Summary Data
save([summary_dir filesep 'summary_results.mat'],'summary')

%% Delete post processed middle man opensees data
rmdir(opensees_dir, 's')
end

