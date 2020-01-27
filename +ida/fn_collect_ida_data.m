function [ ] = fn_collect_ida_data(analysis, model, gm_set_table, ida_results, main_dir, write_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import packages
import plotting_tools.*

% Defined fixed parames
% params = {'b','e','io','ls','cp','euro_th_NC','euro_th_SD','euro_th_DL'};
% params = {'b','io','ls','cp','euro_th_NC','euro_th_SD','euro_th_DL'};
params = {'b','io','ls','cp'};

% Load model data
model_dir = [main_dir '/' 'opensees_data'];
asce41_dir = [main_dir '/' 'asce_41_data'];
load([asce41_dir filesep 'element_analysis.mat']);

% Collect IDA data
id = 0;
id_missing = 0;
for gm = 1:height(gm_set_table)
    gm_dir = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_table.set_id(gm)) '_' num2str(gm_set_table.pair(gm))];
    scale_folders = dir([gm_dir filesep 'Scale_*']);
    for s = 1:length(scale_folders)
        % Load data
        outputs_dir = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_table.set_id(gm)) '_' num2str(gm_set_table.pair(gm)) '/' scale_folders(s).name];
        outputs_file = [outputs_dir filesep 'summary_results.mat'];
        hinge_file = [outputs_dir filesep 'hinge_analysis.mat'];
        story_file = [outputs_dir filesep 'story_analysis.mat'];
        if exist(outputs_file,'file') && exist(hinge_file,'file')
            load(outputs_file)
            load(hinge_file)
            load(story_file)
            if isfield(summary,'max_drift_x')
            id = id + 1;
            ida.id(id,1) = id;
            ida.eq_name{id,1} = gm_set_table.eq_name{gm};
            ida.scale(id,1) = str2double(regexp(scale_folders(s).name,'(?<=_).+$','match'));
            
            if analysis.run_z_motion
                ida.max_drift(id,1) = max(summary.max_drift_x,summary.max_drift_z);
            else
                ida.max_drift(id,1) = summary.max_drift_x;
            end
            
            % X direction
            ida.sa_x(id,1) = summary.sa_x;
            ida.mce_ratio_x(id,1) = ida.sa_x(id,1)/ida_results.mce(1);
            ida.drift_x(id,1) = summary.max_drift_x;
            
            % z direction 
            if analysis.run_z_motion
                ida.sa_z(id,1) = summary.sa_z;
                ida.mce_ratio_z(id,1) = ida.sa_z(id,1)/ida_results.mce(2);
                ida.drift_z(id,1) = summary.max_drift_z;
            end
            
            % Collapse metrics
            ida.collapse(id,1) = summary.collapse;
            if summary.collapse > 0
                ida.collapse_direction{id,1} = summary.collapse_direction;
                ida.collapse_mech{id,1} = summary.collaspe_mech;
                load(hinge_file)
                if contains(summary.collaspe_mech,'column')
                    ida.collapse_comps(id,1) = sum(hinge.b_ratio(strcmp(hinge.ele_type,'column')) >= 1);
                else
                    ida.collapse_comps(id,1) = sum(hinge.b_ratio >= 1);
                end
            else
                ida.collapse_direction{id,1} = 'NA';
                ida.collapse_mech{id,1} = 'NA';
                ida.collapse_comps(id,1) = NaN;
            end
            
%            % Dissapated Energy
%            col_hinges = hinge(hinge.story == 1 & strcmp(hinge.direction,'primary') & strcmp(hinge.ele_type,'column'),:);
%            in_plane_col_base = hinge(hinge.story == 1 & hinge.ele_side == 1 & strcmp(hinge.direction,'primary') & strcmp(hinge.ele_type,'column'),:);
%            
%            ew_pushover_energy = 0;
%            for e = 1:height(col_hinges)
%                pushover_TH = load([pushover_dir filesep 'hinge_TH_' num2str(col_hinges.id(e)) '.mat']);
%                ew_pushover_energy = max(pushover_TH.hin_TH.energy_ft_lbs) + ew_pushover_energy;
%            end
% %            ns_pushover_energy = 0;
% %            for e = 1:height(wall_hinges)
% %                pushover_TH = load([pushover_dir filesep 'hinge_TH_' num2str(wall_hinges.id(e)) '.mat']);
% %                ns_pushover_energy = max(pushover_TH.hin_TH.energy_ft_lbs) + ns_pushover_energy;
% %            end
%            ida.total_energy_ft_lbs(id,1) = sum(cols_walls_1_hinges.total_engergy_ft_lbs);
%            ida.total_energy_ew(id,1) = sum(col_hinges.total_engergy_ft_lbs);
% %            ida.total_energy_ns(id,1) = sum(wall_hinges.total_engergy_ft_lbs);
% %            ida.norm_energy_tot(id,1) = ida.total_energy_ft_lbs(id,1) / (ew_pushover_energy + ns_pushover_energy);
%            ida.norm_energy_tot(id,1) = ida.total_energy_ft_lbs(id,1) / (ew_pushover_energy);
%            ida.norm_energy_ew(id,1) = ida.total_energy_ew(id,1) / ew_pushover_energy;
% %            ida.norm_energy_ns(id,1) = ida.total_energy_ns(id,1) / ns_pushover_energy;
% %            ida.norm_energy_max(id,1) = max([ida.norm_energy_ew(id,1),ida.norm_energy_ns(id,1)]);
%            ida.norm_energy_max(id,1) = ida.norm_energy_ew(id,1);
           
            % Gravity and Lateral Capacity Remaining
            for i = 1:height(story)
                gravity_load(i) = sum(story.story_dead_load(i:end)) + sum(story.story_live_load(i:end));
                col_hinges_1 = hinge(hinge.story == i & strcmp(hinge.ele_type,'column') & hinge.ele_side == 1,:);
                col_hinges_2 = hinge(hinge.story == i & strcmp(hinge.ele_type,'column') & hinge.ele_side == 2,:);
                if ~isempty(col_hinges_1)
                    grav_cap_1 = sum(col_hinges_1.P_capacity);
                    grav_cap_2 = sum(col_hinges_2.P_capacity);
                    gravity_capacity(i) = min([grav_cap_1,grav_cap_2]);
                else
                    gravity_capacity(i) = inf;
                end
            end
            axial_dcr = gravity_load ./ gravity_capacity;
            ida.gravity_dcr(id,1) = max(axial_dcr);
            
            % Lateral Capacity Remaining (first story)
            % full strength
            all_cols = element(element.story == 1 & strcmp(element.type,'column'),:);
            lat_cap_tot = sum(all_cols.Mn_pos_1) + sum(all_cols.Mn_pos_2);
            
            % remaining hinges and columns
            failed_col_any = hinge.element_id(hinge.story == 1 & strcmp(hinge.ele_type,'column') & hinge.a_ratio > 1);
            failed_col_1 = hinge.element_id(hinge.story == 1 & strcmp(hinge.ele_type,'column') & hinge.ele_side == 1 & hinge.a_ratio > 1);
            failed_col_2 = hinge.element_id(hinge.story == 1 & strcmp(hinge.ele_type,'column') & hinge.ele_side == 2 & hinge.a_ratio > 1);
            remianing_cols = element(element.story == 1 & strcmp(element.type,'column') & ~ismember(element.id,unique(failed_col_any)),:);
            remianing_cols_1 = element(element.story == 1 & strcmp(element.type,'column') & ~ismember(element.id,failed_col_1),:);
            remianing_cols_2 = element(element.story == 1 & strcmp(element.type,'column') & ~ismember(element.id,failed_col_2),:);
            
            % remaining columns with no damage
            lat_cap_remain = sum(remianing_cols.Mn_pos_1) + sum(remianing_cols.Mn_pos_2);
            ida.lat_cap_ratio_any(id,1) = lat_cap_remain / lat_cap_tot;
            
            % remaining columns with only 1 side damaged
            lat_cap_remain = sum(remianing_cols_1.Mn_pos_1) + sum(remianing_cols_2.Mn_pos_2);
            ida.lat_cap_ratio_both(id,1) = lat_cap_remain / lat_cap_tot;
            
            % max side of columns remaining
            lat_cap_remain = max([sum(remianing_cols_1.Mn_pos_1),sum(remianing_cols_2.Mn_pos_2)]);
            ida.lat_cap_ratio_max(id,1) = lat_cap_remain / (lat_cap_tot/2); % only works when columns are symmetric  
            
            % how many full column failures are there
            ida.num_full_col_fails(id,1) = sum(~ismember(all_cols.id,remianing_cols_1.id) & ~ismember(all_cols.id,remianing_cols_2.id));
            
            else
                id_missing = id_missing + 1;
                missing_ida.eq_name{id_missing,1} = gm_set_table.eq_name{gm};
                missing_ida.scale(id_missing,1) = str2double(regexp(scale_folders(s).name,'(?<=_).+$','match'));
            end
        else
            id_missing = id_missing + 1;
            missing_ida.eq_name{id_missing,1} = gm_set_table.eq_name{gm};
            missing_ida.scale(id_missing,1) = str2double(regexp(scale_folders(s).name,'(?<=_).+$','match'));
        end
    end
end

ida_table = struct2table(ida);

% Remove all cases that failed to converge yet did not get far enough
failed_convergence = ida_table(ida_table.collapse == 5,:);
ida_table(ida_table.collapse == 5,:) = []; % filter out failed models

% Go through and define collapse props for each ground motion
for gm = 1:height(gm_set_table)
    filt_collapse = ida_table.collapse > 0 & strcmp(ida_table.eq_name,gm_set_table.eq_name{gm});
    
    % Remove all GM's that do not collapse
    if sum(filt_collapse) == 0
        ida_table(strcmp(ida_table.eq_name,gm_set_table.eq_name{gm}),:) = [];
    else
        collapse_idx = find(filt_collapse,1,'first');
        gm_set_table.scale_collapse(gm) = ida_table.scale(collapse_idx);
        gm_set_table.scale_jbc(gm) = ida_table.scale(collapse_idx - 1);
        gm_set_table.collapse_dir(gm) = ida_table.collapse_direction(collapse_idx);
        gm_set_table.collapse_mech(gm) = ida_table.collapse_mech(collapse_idx);
        gm_set_table.collapse_comps(gm) = ida_table.collapse_comps(collapse_idx);
    end
end
gm_set_table(gm_set_table.scale_collapse == 0,:) = [];

% Go through IDA table and get hinge properties based on collapse mechanism
for i = 1:height(ida_table)
    outputs_dir = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_table.set_id(strcmp(gm_set_table.eq_name,ida_table.eq_name{i}))) '_' num2str(gm_set_table.pair(strcmp(gm_set_table.eq_name,ida_table.eq_name{i}))) '/Scale_' num2str(ida_table.scale(i))];
    hinge_file = [outputs_dir filesep 'hinge_analysis.mat'];
    load(hinge_file)
    
    % Identify collapse case for this GM
    gm_stripes = ida_table(strcmp(ida_table.eq_name,ida_table.eq_name{i}),:);
    collapse_idx = find(gm_stripes.collapse > 0,1,'first');
    collapse_hinge_file = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_table.set_id(strcmp(gm_set_table.eq_name,ida_table.eq_name{i}))) '_' num2str(gm_set_table.pair(strcmp(gm_set_table.eq_name,ida_table.eq_name{i}))) '/Scale_' num2str(gm_stripes.scale(collapse_idx)) filesep 'hinge_analysis.mat'];
    collapse_hinge = load(collapse_hinge_file);
    
    % Get element group filters from collapse case
%     if contains(gm_stripes.collapse_mech{collapse_idx},'column')
%         % collapse mechanism is all column hinges that have failed in the first stripe that collapse
%         mech_filter = strcmp(collapse_hinge.hinge.direction,'primary') & collapse_hinge.hinge.b_ratio >= 1 & strcmp(collapse_hinge.hinge.ele_type,'column'); 
%     else
%         % collapse mechanism is all hinges that have failed in the first stripe that collapse
%         mech_filter = strcmp(collapse_hinge.hinge.direction,'primary') & collapse_hinge.hinge.b_ratio >= 1; 
%     end
    
    mech_filter = strcmp(collapse_hinge.hinge.direction,'primary') & collapse_hinge.hinge.b_ratio >= 1;
    
    ida_table.num_comps(i) = sum(mech_filter);
    mech_hinges = hinge(mech_filter,:);

    % For Each accetance criteria listed above
    for p = 1:length(params)
        [ num_eles, percent_eles, num_eles_15, max_ele, min_ele, mean_ele, range_ele, std_ele, cov_ele ] = fn_collect_component_data(params{p}, ida_table.collapse(i), ida_table.collapse_direction{i}, mech_hinges);
        ida_table.(['num_' params{p}])(i) = num_eles;
        ida_table.(['percent_' params{p}])(i) = percent_eles;
        ida_table.(['num_' params{p} '_15'])(i) = num_eles_15;
        ida_table.(['max_' params{p}])(i) = max_ele;
        ida_table.(['min_' params{p}])(i) = min_ele;
        ida_table.(['mean_' params{p}])(i) = mean_ele;
        ida_table.(['range_' params{p}])(i) = range_ele;
        ida_table.(['std_' params{p}])(i) = std_ele;
        ida_table.(['cov_' params{p}])(i) = cov_ele;
    end

%     % Unacceptable Response
%     if ida_table.collapse(i) == 1 || ida_table.collapse(i) == 3 || ida_table.num_b_15(i) > 0
%         ida_table.UR(i) = 1;
%     else
%         ida_table.UR(i) = 0;
%     end
% 
%     % Gravity Load Lost
%     first_story_elements = element(ismember(element.id, mech_hinges.element_id),:);
%     hinges_lost_grav = mech_hinges(mech_hinges.b_ratio > 1,:);
%     elements_lost_grav = element(ismember(element.id, hinges_lost_grav.element_id),:);
%     grav_load_lost = sum(elements_lost_grav.P_grav);
%     total_grav_load = sum(first_story_elements.P_grav);
%     ida_table.gravity_load_lost_ratio(i) = grav_load_lost / total_grav_load;
end

% Save Tabular Results as CSVs
writetable(gm_set_table,[write_dir filesep 'gm_table.csv'])
writetable(ida_table,[write_dir filesep 'ida_table.csv'])
if exist('missing_ida','var')
    writetable(struct2table(missing_ida),[write_dir filesep 'idas_missing.csv'])
end
writetable(failed_convergence,[write_dir filesep 'idas_failed_convergence.csv'])

end

function [ num_eles, percent_eles, num_eles_15, max_ele, min_ele, mean_ele, range_ele, std_ele, cov_ele ] = fn_collect_component_data(var_name, collapse_flag, collaspe_dir, ele_hinges)

ele_ratios = ele_hinges.([var_name '_ratio']);

num_eles = 0;
num_eles_15 = 0;
for e = 1:length(ele_ratios)
    if (collapse_flag == 3 || collapse_flag == 1)
        if strcmp(ele_hinges.ele_direction(e),collaspe_dir)
            num_eles = num_eles + 1; % if collapse this gm, in this direction, set this element to 1
            num_eles_15 = num_eles_15 + 1;
        end
    elseif ele_ratios(e) >= 1.5
        num_eles_15 = num_eles_15 + 1;
        num_eles = num_eles + 1;
    elseif ele_ratios(e) >= 1
        num_eles = num_eles + 1;
    end
end
percent_eles = num_eles / length(ele_ratios);
max_ele = max(ele_ratios);
min_ele = min(ele_ratios);
mean_ele = mean(ele_ratios);
range_ele = max_ele - min_ele;
std_ele = std(ele_ratios);
cov_ele = std_ele/mean_ele;
end