function [ ] = fn_collect_ida_data(analysis, model, gm_set_table, ida_results, write_dir, pushover_dir, model_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import packages
import plotting_tools.*

% Defined fixed parames
params = {'b','e','b_e','io','ls','cp','euro_th_NC','euro_th_SD','euro_th_DL'};
mechs = { 'cols_1', 'walls_1', 'cols_walls_1'};

% Load model data
if ~exist('pushover_dir','var')
    model_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'opensees_data'];
end
if ~exist('pushover_dir','var')
    pushover_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'pushover'];
end
load([model_dir filesep 'element_analysis.mat']);
load([model_dir filesep 'node_analysis.mat']);

% Collect IDA data
id = 0;
id_missing = 0;
for gm = 1:height(gm_set_table)
    gm_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_table.set_id(gm)) '_' num2str(gm_set_table.pair(gm))];
    scale_folders = dir([gm_dir filesep 'Scale_*']);
    for s = 1:length(scale_folders)
        % Load data
        outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_table.set_id(gm)) '_' num2str(gm_set_table.pair(gm)) '/' scale_folders(s).name];
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
            ida.max_drift(id,1) = max(summary.max_drift_x,summary.max_drift_z);
            
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
            if isfield(summary,'collapse_direction')
                ida.collapse_direction{id,1} = summary.collapse_direction;
            else
                summary.collapse_direction = 'NA';
            end
            if isfield(summary,'collaspe_mech')
                ida.collapse_mech{id,1} = summary.collaspe_mech;
            end
            
            % Get element group filters
            first_story_col_filter = hinge.story == 1 & strcmp(hinge.direction,'primary') & strcmp(hinge.ele_type,'column');
            first_story_wall_filter = hinge.story == 1 & strcmp(hinge.direction,'primary')  & strcmp(hinge.ele_type,'wall');
            ida.num_comps_cols_1(id,1) = sum(first_story_col_filter);
            ida.num_comps_walls_1(id,1) = sum(first_story_wall_filter);
            ida.num_comps_cols_walls_1(id,1) = sum(first_story_col_filter | first_story_wall_filter);
            mech_hinges{1}.prime = hinge(first_story_col_filter,:);
            mech_hinges{1}.second = [];
            mech_hinges{2}.prime = hinge(first_story_wall_filter,:);
            mech_hinges{2}.second = [];
            mech_hinges{3}.prime = hinge(first_story_col_filter,:);
            mech_hinges{3}.second = hinge(first_story_wall_filter,:);

            % For Each accetance criteria listed above
            for m = 1:length(mechs)
                for p = 1:length(params)
                    [ num_eles, percent_eles, num_eles_15, max_ele, min_ele, mean_ele, range_ele, std_ele, cov_ele ] = fn_collect_component_data(params{p}, summary.collapse, summary.collapse_direction, mech_hinges{m}.prime, mech_hinges{m}.second);
                    ida.([mechs{m} '_num_' params{p}])(id,1) = num_eles;
                    ida.([mechs{m} '_percent_' params{p}])(id,1) = percent_eles;
                    ida.([mechs{m} '_num_' params{p} '_15'])(id,1) = num_eles_15;
                    ida.([mechs{m} '_max_' params{p}])(id,1) = max_ele;
                    ida.([mechs{m} '_min_' params{p}])(id,1) = min_ele;
                    ida.([mechs{m} '_mean_' params{p}])(id,1) = mean_ele;
                    ida.([mechs{m} '_range_' params{p}])(id,1) = range_ele;
                    ida.([mechs{m} '_std_' params{p}])(id,1) = std_ele;
                    ida.([mechs{m} '_cov_' params{p}])(id,1) = cov_ele;
                end
            end
            
            % Unacceptable Response
            if summary.collapse == 1 || summary.collapse == 3 || ida.cols_walls_1_num_b_e_15(id,1) > 0
                ida.UR(id,1) = 1;
            else
                ida.UR(id,1) = 0;
            end
            
            % Test hinge values
            wall_hinges = hinge(first_story_wall_filter,:);
            col_hinges = hinge(hinge.story == 1 & strcmp(hinge.direction,'primary') & strcmp(hinge.ele_type,'column'),:);
            in_plane_col_base = hinge(hinge.story == 1 & hinge.ele_side == 1 & strcmp(hinge.direction,'primary') & strcmp(hinge.ele_type,'column'),:);
            in_plane_col_top = hinge(hinge.story == 1 & hinge.ele_side == 2 & strcmp(hinge.direction,'primary') & strcmp(hinge.ele_type,'column'),:);
            oop_col_base = hinge(hinge.story == 1 & hinge.ele_side == 1 & strcmp(hinge.direction,'oop') & strcmp(hinge.ele_direction,'x'),:);
            oop_col_top = hinge(hinge.story == 1 & hinge.ele_side == 2 & strcmp(hinge.direction,'oop'),:);
            
            % Gravity Load Lost
            cols_walls_1_hinges = hinge(first_story_col_filter | first_story_wall_filter,:);
            first_story_elements = element(ismember(element.id, cols_walls_1_hinges.element_id),:);
            hinges_lost_grav = cols_walls_1_hinges(cols_walls_1_hinges.b_ratio > 1 | cols_walls_1_hinges.e_ratio > 1,:);
            wall_fails = cols_walls_1_hinges(cols_walls_1_hinges.e_ratio > 1,:);
            side_1_fails = cols_walls_1_hinges(cols_walls_1_hinges.b_ratio > 1 & cols_walls_1_hinges.ele_side == 1,:);
            side_2_fails = cols_walls_1_hinges(cols_walls_1_hinges.b_ratio > 1 & cols_walls_1_hinges.ele_side == 2,:);
            elements_lost_grav = element(ismember(element.id, hinges_lost_grav.element_id),:);
            elements_lost_grav_2 = element(ismember(element.id, wall_fails.element_id) | (ismember(element.id, side_1_fails.element_id) & ismember(element.id, side_2_fails.element_id)),:);
            grav_load_lost = sum(elements_lost_grav.P_grav);
            grav_load_lost_2 = sum(elements_lost_grav_2.P_grav);
            total_grav_load = sum(first_story_elements.P_grav);
%             total_grav_load = sum(story.story_dead_load + story.story_live_load); % take the whole build wt since I am assessing the first story
            ida.gravity_load_lost_ratio(id,1) = grav_load_lost / total_grav_load;
            ida.gravity_load_lost_ratio_alt(id,1) = grav_load_lost_2 / total_grav_load;

            % Adjacent components (focus on just columns for now)
            cols_walls_1_hinge_nodes = node(ismember(node.id,cols_walls_1_hinges.node_1),:);
            ida.adjacent_failure_any(id,1) = 0;
            ida.adjacent_failure_any_frame(id,1) = 0;
            ida.adjacent_failure_all(id,1) = 0;
            col_hinges_fail = hinges_lost_grav(strcmp(hinges_lost_grav.ele_type,'column'),:);
            ida.num_adjacent_failure_any(id,1) = 0;
            ida.num_adjacent_failure_any_frame(id,1) = 0;
            ida.num_adjacent_failure_all(id,1) = 0;
            for h = 1:height(col_hinges_fail)
                hin = col_hinges_fail(h,:);
                hin_x = node.x(node.id == hin.node_1); 
                hin_z = node.z(node.id == hin.node_1);
                node_east = cols_walls_1_hinge_nodes(cols_walls_1_hinge_nodes.x == hin_x + 300 & cols_walls_1_hinge_nodes.z == hin_z,:);
                node_west = cols_walls_1_hinge_nodes(cols_walls_1_hinge_nodes.x == hin_x - 300 & cols_walls_1_hinge_nodes.z == hin_z,:);
                node_north = cols_walls_1_hinge_nodes(cols_walls_1_hinge_nodes.x == hin_x & cols_walls_1_hinge_nodes.z == hin_z + 300,:);
                node_south = cols_walls_1_hinge_nodes(cols_walls_1_hinge_nodes.x == hin_x & cols_walls_1_hinge_nodes.z == hin_z - 300,:);
                hin_east = cols_walls_1_hinges(ismember(cols_walls_1_hinges.node_1,node_east.id) | ismember(cols_walls_1_hinges.node_2,node_east.id),:);
                hin_west = cols_walls_1_hinges(ismember(cols_walls_1_hinges.node_1,node_west.id) | ismember(cols_walls_1_hinges.node_2,node_west.id),:);
                hin_north = cols_walls_1_hinges(ismember(cols_walls_1_hinges.node_1,node_north.id) | ismember(cols_walls_1_hinges.node_2,node_north.id),:);
                hin_south = cols_walls_1_hinges(ismember(cols_walls_1_hinges.node_1,node_south.id) | ismember(cols_walls_1_hinges.node_2,node_south.id),:);
                
                % At least 1 adjacent component
                if any(ismember(col_hinges_fail.id,[hin_east.id; hin_west.id; hin_north.id; hin_south.id]))
                    ida.adjacent_failure_any(id,1) = 1;
                    ida.num_adjacent_failure_any(id,1) = 1 + ida.num_adjacent_failure_any(id,1);
                end
                
                % At least 1 adjacent component in frame line
                if any(ismember(col_hinges_fail.id,[hin_east.id; hin_west.id]))
                    ida.adjacent_failure_any_frame(id,1) = 1;
                    ida.num_adjacent_failure_any_frame(id,1) = 1 + ida.num_adjacent_failure_any_frame(id,1);
                end

                % All adjacent components 
               if any(ismember(col_hinges_fail.id,hin_east.id)) && any(ismember(col_hinges_fail.id,hin_west.id)) &&  any(ismember(col_hinges_fail.id,hin_north.id)) && any(ismember(col_hinges_fail.id,hin_south.id))
                    ida.adjacent_failure_all(id,1) = 1;
                    ida.num_adjacent_failure_all(id,1) = 1 + ida.num_adjacent_failure_all(id,1);
               end
            end
            
           % Dissapated Energy
           ew_pushover_energy = 0;
           for e = 1:height(col_hinges)
               pushover_TH = load([pushover_dir filesep 'hinge_TH_' num2str(col_hinges.id(e)) '.mat']);
               ew_pushover_energy = max(pushover_TH.hin_TH.energy_ft_lbs) + ew_pushover_energy;
           end
           ns_pushover_energy = 0;
           for e = 1:height(wall_hinges)
               pushover_TH = load([pushover_dir filesep 'hinge_TH_' num2str(wall_hinges.id(e)) '.mat']);
               ns_pushover_energy = max(pushover_TH.hin_TH.energy_ft_lbs) + ns_pushover_energy;
           end
           ida.total_energy_ft_lbs(id,1) = sum(cols_walls_1_hinges.total_engergy_ft_lbs);
           ida.total_energy_ew(id,1) = sum(col_hinges.total_engergy_ft_lbs);
           ida.total_energy_ns(id,1) = sum(wall_hinges.total_engergy_ft_lbs);
           ida.norm_energy_tot(id,1) = ida.total_energy_ft_lbs(id,1) / (ew_pushover_energy + ns_pushover_energy);
           ida.norm_energy_ew(id,1) = ida.total_energy_ew(id,1) / ew_pushover_energy;
           ida.norm_energy_ns(id,1) = ida.total_energy_ns(id,1) / ns_pushover_energy;
           ida.norm_energy_max(id,1) = max([ida.norm_energy_ew(id,1),ida.norm_energy_ns(id,1)]);
           
           % Gravity Capacity Remaining
           gravity_capacity = 24*24*7.5*(24 - sum(in_plane_col_base.b_ratio >= 1));
           ida.gravity_dcr(id,1) = 11841/gravity_capacity;
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

% filter non_collapse 
% Remove all cases that failed to converge yet did not get far enough
ida_table = struct2table(ida);
failed_convergence = ida_table(ida_table.collapse == 5,:);
ida_table(ida_table.collapse == 5,:) = []; % filter out failed models

% Save Tabular Results as CSVs
writetable(ida_table,[write_dir filesep 'ida_table.csv'])
if exist('missing_ida','var')
    writetable(struct2table(missing_ida),[write_dir filesep 'idas_missing.csv'])
end
writetable(failed_convergence,[write_dir filesep 'idas_failed_convergence.csv'])

end

function [ num_eles, percent_eles, num_eles_15, max_ele, min_ele, mean_ele, range_ele, std_ele, cov_ele ] = fn_collect_component_data(var_name, collapse_flag, collaspe_dir, ele_hinges, ele_hinges_alt)

ele_ratios_alt = [];

if strcmp(var_name,'b_e')
    ele_ratios = ele_hinges.b_ratio; % combo of both b and e values
    if ~isempty(ele_hinges_alt)
        ele_ratios_alt = ele_hinges_alt.e_ratio; % combo of both b and e values
    end
else
    ele_ratios = ele_hinges.([var_name '_ratio']);
    if ~isempty(ele_hinges_alt)
        ele_ratios_alt = ele_hinges_alt.([var_name '_ratio']);
    end
end

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

% Additional criteria (ie e ratio)
num_eles_alt = 0;
num_eles_15_alt = 0;
for e = 1:length(ele_ratios_alt)
    if (collapse_flag == 3 || collapse_flag == 1)
        if strcmp(ele_hinges_alt.ele_direction(e),collaspe_dir)
            num_eles_alt = num_eles_alt + 1; % if collapse this gm, in this direction, set this element to 1
            num_eles_15_alt = num_eles_15_alt + 1;
        end
    elseif ele_ratios_alt(e) >= 1.5
        num_eles_alt = num_eles_alt + 1;
        num_eles_15_alt = num_eles_15_alt + 1;
    elseif ele_ratios_alt(e) >= 1
        num_eles_alt = num_eles_alt + 1;
    end
end
num_eles = num_eles + num_eles_alt;
num_eles_15 = num_eles_15 + num_eles_15_alt;
percent_eles_alt = num_eles_alt / length(ele_ratios_alt);
percent_eles = max(percent_eles,percent_eles_alt);

max_ele = max([ele_ratios;ele_ratios_alt]);
min_ele = min([ele_ratios;ele_ratios_alt]);
mean_ele = mean([ele_ratios;ele_ratios_alt]);
range_ele = max_ele - min_ele;
std_ele = std([ele_ratios;ele_ratios_alt]);
cov_ele = std_ele/mean_ele;
end