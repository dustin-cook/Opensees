function [ ] = fn_collect_and_plot_ida_simple(analysis, gm_set_table, node, story, main_dir, write_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import packages
import plotting_tools.*

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
        if exist(outputs_file,'file')
            % Load summary data
            load(outputs_file)
            
            if isfield(summary,'node')

            % Define EQ and scale data in table
            id = id + 1;
            ida.id(id,1) = id;
            ida.eq_name{id,1} = gm_set_table.eq_name{gm};
            ida.scale(id,1) = str2double(regexp(scale_folders(s).name,'(?<=_).+$','match'));
           
            % Define shaking
            ida.sa_x(id,1) = summary.sa_x;
            ida.sa_geo(id,1) = summary.sa_x;
            if analysis.run_z_motion
                ida.sa_z(id,1) = summary.sa_z;
                ida.sa_geo(id,1) = geomean([summary.sa_x,summary.sa_z]);
            end
            
            % Summarize story response (Currently only 
            summary.story = table;
            for str = 1:length(story.story_ht)
                summary.story.id(str) = str;
                summary.story.story_ht(str) = story.story_ht(str);
                node_above = node.primary_nodes(str);
                
                % Disp above
                node_above_disp_x = summary.node.max_disp_x(summary.node.id == node_above);
                if analysis.run_z_motion
                    node_above_disp_z = summary.node.max_disp_z(summary.node.id == node_above);
                end
                
                % disp below
                if str == 1
                    node_below_disp_x = 0;
                    node_below_disp_z = 0;
                else
                    node_below = node.primary_nodes(str-1);
                    node_below_disp_x = summary.node.max_disp_x(summary.node.id == node_below);
                    if analysis.run_z_motion
                        node_below_disp_z = summary.node.max_disp_z(summary.node.id == node_below);
                    end
                end
                
                summary.story.peak_drift_x(str) = (node_above_disp_x - node_below_disp_x)/story.story_ht(str);
                if analysis.run_z_motion
                    summary.story.peak_drift_z(str) = (node_above_disp_z - node_below_disp_z)/story.story_ht(str);
                end
            end
            
            % Save summary data
            save(outputs_file, 'summary');
            
            % Define Drift
            if analysis.run_z_motion
                ida.drift_x(id,1) = max(summary.story.peak_drift_x);
                ida.drift_z(id,1) = max(summary.story.peak_drift_z);
                ida.max_drift(id,1) = max([summary.story.peak_drift_x; summary.story.peak_drift_z]);
                ida.max_drift_story(id,1) = summary.story.id(max(summary.story.peak_drift_x, summary.story.peak_drift_x) == ida.max_drift(id,1));
                if ida.drift_x(id,1) < ida.drift_z(id,1)
                    ida.max_drift_dir{id,1} = 'z';
                else
                    ida.max_drift_dir{id,1} = 'x';
                end
            else
                ida.drift_x(id,1) = max(summary.story.peak_drift_x);
                ida.max_drift(id,1) = max(summary.story.peak_drift_x);
                ida.max_drift_story(id,1) = summary.story.id(summary.story.peak_drift_x == ida.max_drift(id,1));
            end

            % Collapse metrics
            ida.collapse(id,1) = summary.collapse;
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
        gm_set_table.sa_collapse(gm) = ida_table.sa_x(collapse_idx);
        gm_set_table.scale_jbc(gm) = ida_table.scale(collapse_idx - 1);
        gm_set_table.sa_jbc(gm) = ida_table.sa_x(collapse_idx - 1);
    end
end
gm_set_table(gm_set_table.scale_collapse == 0,:) = [];

% Save Tabular Results as CSVs
writetable(gm_set_table,[write_dir filesep 'gm_table.csv'])
writetable(ida_table,[write_dir filesep 'ida_table.csv'])
if exist('missing_ida','var')
    writetable(struct2table(missing_ida),[write_dir filesep 'idas_missing.csv'])
end
writetable(failed_convergence,[write_dir filesep 'idas_failed_convergence.csv'])

% Plot IDA curves
hold on
for gms = 1:height(gm_set_table)
    ida_plt = plot(ida_table.max_drift(strcmp(ida_table.eq_name,gm_set_table.eq_name{gms})),ida_table.sa_x(strcmp(ida_table.eq_name,gm_set_table.eq_name{gms})),'-o','color',[0.75 0.75 0.75],'HandleVisibility','off');
end
xlim([0 0.1])
xlabel('Max Drift')
ylabel('Sa(T_1) (g)')
fn_format_and_save_plot( write_dir,'IDA Plot', 4 )

% Create Collapse Fragility Curve
sa_dist = gm_set_table.sa_collapse;
sa_dist(isnan(sa_dist)) = [];
sa_dist(sa_dist <= 0) = 1e-6;
if ~isempty(sa_dist)
    [pHat, ~] = lognfit(sa_dist);
    collapse.theta = exp(pHat(1));
    collapse.beta = pHat(2);
else
    collapse.theta = NaN;
    collapse.beta = NaN;
end

% Plot Collapse Fragility 
x_points = 0.01:0.01:4;
hold on
rank_sa = sort(sa_dist);
rank_val = 1:length(rank_sa);
sct = scatter(rank_sa,rank_val./length(rank_sa),'k','filled','HandleVisibility','off');
cdf = logncdf(x_points,log(collapse.theta),collapse.beta);
plot(x_points,cdf,'k','DisplayName','Collapse')
xlabel('Sa(T_1) (g)')
ylabel('P[Collapse]')
xlim([0,1])
fn_format_and_save_plot( write_dir, 'Collapse Fragility', 4 )

end
