function [ ] = fn_create_fragilities(analysis, model, IDA_scale_factors, gm_set_table, max_dir_spectra, ida_results)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import packages
import plotting_tools.*

% Defined fixed parames
params = {'b','e','b_e','io','ls','cp','euro_th_NC','euro_th_SD','euro_th_DL'};
mechs = { 'cols_1', 'walls_1', 'cols_walls_1'};
frag_probs = [10 25 50 75 100];

% Load model data
model_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'opensees_data'];
load([model_dir filesep 'element_analysis.mat']);
load([model_dir filesep 'node_analysis.mat']);

% Collect IDA data
write_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Fragility Data'];
if ~exist(write_dir,'dir')
    mkdir(write_dir)
end
id = 0;
id_missing = 0;
for i = 1:length(IDA_scale_factors)
    for gms = 1:height(gm_set_table)   
        % Load data
        outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id '/' 'IDA' '/' 'Summary Data' '/' 'Scale_' num2str(IDA_scale_factors(i)) '/' 'GM_' num2str(gm_set_table.set_id(gms)) '_' num2str(gm_set_table.pair(gms))];
        outputs_file = [outputs_dir filesep 'summary_results.mat'];
        hinge_file = [outputs_dir filesep 'hinge_analysis.mat'];
        story_file = [outputs_dir filesep 'story_analysis.mat'];
        if exist(outputs_file,'file') && exist(hinge_file,'file')
            id = id + 1;
            load(outputs_file)
            load(hinge_file)
            load(story_file)
            ida.id(id,1) = id;
            ida.scale(id,1) = IDA_scale_factors(i);
            ida.gm_name{id,1} = gm_set_table.eq_name{gms};
            
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
            if i > 1
                gm_prev_stripe_collapse = ida.collapse(strcmp(ida.gm_name(1:end-1),gm_set_table.eq_name{gms}));
                if summary.collapse > 0 && any(gm_prev_stripe_collapse)
                    ida.collapse_prev_stripe(id,1) = 1;
                else
                    ida.collapse_prev_stripe(id,1) = 0;
                end
            end
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
%             first_story_beam_filter = hinge.story == 1 & strcmp(hinge.direction,'primary')  & strcmp(hinge.ele_type,'beam');
            first_story_wall_filter = hinge.story == 1 & strcmp(hinge.direction,'primary')  & strcmp(hinge.ele_type,'wall');
            num_comps_cols_1 = sum(first_story_col_filter);
%             num_comps_bms_1 = sum(first_story_beam_filter);
            num_comps_walls_1 = sum(first_story_wall_filter);
            num_comps_cols_walls_1 = sum(first_story_col_filter | first_story_wall_filter);
            mech_hinges{1}.prime = hinge(first_story_col_filter,:);
            mech_hinges{1}.second = [];
%             mech_hinges{2}.prime = hinge(first_story_beam_filter,:);
%             mech_hinges{2}.second = [];
            mech_hinges{2}.prime = hinge(first_story_wall_filter,:);
            mech_hinges{2}.second = [];
            mech_hinges{3}.prime = hinge(first_story_col_filter,:);
            mech_hinges{3}.second = hinge(first_story_wall_filter,:);

            % For Each accetance criteria listed above
            for m = 1:length(mechs)
                for p = 1:length(params)
                    [ num_eles, percent_eles, num_eles_15, max_ele_ratio, mean_ele_ratio ] = fn_collect_ida_data(params{p}, summary.collapse, summary.collapse_direction, mech_hinges{m}.prime, mech_hinges{m}.second);
                    ida.([mechs{m} '_num_' params{p}])(id,1) = num_eles;
                    ida.([mechs{m} '_percent_' params{p}])(id,1) = percent_eles;
                    ida.([mechs{m} '_num_' params{p} '_15'])(id,1) = num_eles_15;
                    ida.([mechs{m} '_max_' params{p}])(id,1) = max_ele_ratio;
                    ida.([mechs{m} '_mean_' params{p}])(id,1) = mean_ele_ratio;
                end
            end
            
            % Gravity Load Lost
            cols_walls_1_hinges = hinge(first_story_col_filter | first_story_wall_filter,:);
            hinges_lost_grav = cols_walls_1_hinges(cols_walls_1_hinges.b_ratio > 1 | cols_walls_1_hinges.e_ratio > 1,:);
            elements_lost_grav = element(ismember(element.id, hinges_lost_grav.element_id),:);
            grav_load_lost = sum(elements_lost_grav.P_grav);
            total_grav_load = sum(story.story_dead_load + story.story_live_load); % take the whole build wt since I am assessing the first story
            ida.gravity_load_lost_ratio(id,1) = grav_load_lost / total_grav_load;
%             cols_walls_1_elements = element(ismember(element.id,mech_hinges{3}.prime.element_id) | ismember(element.id,mech_hinges{3}.second.element_id),:);

            % Adjacent components (focus on just columns for now)
            cols_walls_1_hinge_nodes = node(ismember(node.id,cols_walls_1_hinges.node_1),:);
            ida.adjacent_failure_any(id,1) = 0;
            ida.adjacent_failure_any_frame(id,1) = 0;
            ida.adjacent_failure_all(id,1) = 0;
            col_hinges_fail = hinges_lost_grav(strcmp(hinges_lost_grav.ele_type,'column'),:);
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
                end
                
                % At least 1 adjacent component in frame line
                if any(ismember(col_hinges_fail.id,[hin_east.id; hin_west.id]))
                    ida.adjacent_failure_any_frame(id,1) = 1;
                end

                % All adjacent components 
               if any(ismember(col_hinges_fail.id,hin_east.id)) && any(ismember(col_hinges_fail.id,hin_west.id)) &&  any(ismember(col_hinges_fail.id,hin_north.id)) && any(ismember(col_hinges_fail.id,hin_south.id))
                    ida.adjacent_failure_all(id,1) = 1;
                end
            end
        else
            id_missing = id_missing + 1;
            missing_ida.scale(id_missing,1) = IDA_scale_factors(i);
            missing_ida.gm_set_id(id_missing,1) = gm_set_table.set_id(gms);
            missing_ida.gm_set_pair_id(id_missing,1) = gm_set_table.pair(gms);
        end
    end
end

% filter non_collapse 
% Remove all cases that failed to converge yet did not get far enough
ida_table = struct2table(ida);
failed_convergence = ida_table(ida_table.collapse == 5,:);
ida_table(ida_table.collapse == 5,:) = []; % filter out failed models
ida_table_full = ida_table;
ida_table(ida_table.collapse_prev_stripe == 1,:) = []; % fitler out models where it previously collapsed

% Save Tabular Results as CSVs
writetable(ida_table,[write_dir filesep 'ida_table.csv'])
if exist('missing_ida','var')
    writetable(struct2table(missing_ida),[write_dir filesep 'idas_missing.csv'])
end
writetable(failed_convergence,[write_dir filesep 'idas_failed_convergence.csv'])
writetable(ida_table_full,[write_dir filesep 'ida_table_full.csv'])

% %% Collect Stripe Info
% for i = 1:length(IDA_scale_factors)
%     this_stripe = ida_table(ida_table.scale == IDA_scale_factors(i),:);
%     
%     %% Regular Collapse
%     stripe.num_gms(i,1) = height(this_stripe);
%     stripe.num_gms_drift(i,1) = sum(this_stripe.collapse ~= 3);
%     stripe.num_gms_convergence(i,1) = sum(this_stripe.collapse ~= 1);
%     stripe.num_gms_ew(i,1) = sum(~strcmp(this_stripe.collapse_direction,'z'));
%     stripe.num_gms_ns(i,1) = sum(~strcmp(this_stripe.collapse_direction,'x'));
%     stripe.num_collapse_drift(i,1) = sum(this_stripe.collapse == 1);
%     stripe.num_collapse_convergence(i,1) = sum(this_stripe.collapse == 3);
%     stripe.num_collapse_tot(i,1) = stripe.num_collapse_drift(i,1) + stripe.num_collapse_convergence(i,1);
%     stripe.num_collapse_ew(i,1) =  sum(strcmp(this_stripe.collapse_direction,'x'));
%     stripe.num_collapse_ns(i,1) =  sum(strcmp(this_stripe.collapse_direction,'z'));
%     stripe.num_accept_15(i,1) = sum(this_stripe.cols_walls_1_num_b_e_15 > 0);
%     stripe.num_unacceptable_response(i,1) = sum(this_stripe.collapse == 1 | this_stripe.collapse == 3 | this_stripe.cols_walls_1_num_b_e_15 > 0);
% 
%     %% Drift checks
%     for d = 1:5
%         stripe.(['num_drift_' num2str(d) '_percent'])(i,1) = sum(this_stripe.drift_x >= d/100 | this_stripe.drift_z >= d/100);
%     end
%     
%     %% Gravity Load lost
%     for f = 1:length(frag_probs)
%         stripe.(['num_grav_lost_' num2str(frag_probs(f)) '_percent'])(i,1) = sum(this_stripe.gravity_load_lost_ratio >= frag_probs(f)/100);
%     end
%     
%     %% Adjacent Components
%     stripe.adjacent_failure_any(i,1) = sum(this_stripe.adjacent_failure_any == 1);
%     stripe.adjacent_failure_any_frame(i,1) = sum(this_stripe.adjacent_failure_any_frame == 1);
%     stripe.adjacent_failure_all(i,1) = sum(this_stripe.adjacent_failure_all == 1);
%     
%     %% X Direction
%     % Lognormal Mean Spectral Accelration
%     stripe.p695_sa_ew(i,1) = exp(mean(log(this_stripe.sa_x)));
%     stripe.p_15_sa_ew(i,1) = prctile(this_stripe.sa_x,15);
%     stripe.p_85_sa_ew(i,1) = prctile(this_stripe.sa_x,85);
% 
%     % Max Direction Spectral Acceleration
%     [~,idx] = min(abs(max_dir_spectra.period - ida_results.period(1)));
%     stripe.asce41_sa_ew(i,1) = max_dir_spectra.sa(idx)*IDA_scale_factors(i);
% 
%     % Drift
%     stripe.mean_idr_ew(i,1) = mean(this_stripe.drift_x);
%     
%     %% Z Direction
%     if analysis.run_z_motion
%         % Lognormal Mean Spectral Accelration
%         stripe.p695_sa_ns(i,1) = exp(mean(log(this_stripe.sa_z)));
%         stripe.p_15_sa_ns(i,1) = prctile(this_stripe.sa_z,15);
%         stripe.p_85_sa_ns(i,1) = prctile(this_stripe.sa_z,85);
% 
%         % Max Direction Spectral Acceleration
%         [~,idx_ns] = min(abs(max_dir_spectra.period - ida_results.period(2)));
%         stripe.asce41_sa_ns(i,1) = max_dir_spectra.sa(idx_ns)*IDA_scale_factors(i);
% 
%         % Drift
%         stripe.mean_idr_ns(i,1) = exp(mean(log(this_stripe.drift_z)));
%     end
%     
%     %% Component mechanisms
%     for m = 1:length(mechs)
%         for p = 1:length(params)
%             stripe.([mechs{m} '_' 'num_first_' params{p}])(i,1) = sum(this_stripe.([mechs{m} '_num_' params{p}]) >= 1);
%             for pr = 1:length(frag_probs)
%                 stripe.([mechs{m} '_' 'num_' num2str(frag_probs(pr)) '_' params{p}])(i,1) = sum(this_stripe.([mechs{m} '_percent_' params{p}]) >= frag_probs(pr)/100);
%             end
%         end
%     end
% end
% 
% % Save frag data
% stripe = struct2table(stripe);
% save([write_dir filesep 'stripe_table.mat'],'stripe')

%% Create Fragility Curves based on Baker MLE
% sa_points = stripe.asce41_sa_ew; % fix use EW period based spectral acceleration points

% collape and unnacceptable response
[frag_curves.collapse.theta, frag_curves.collapse.beta] = fn_fit_fragility_MOM(ida_table.sa_x(ida_table.collapse > 0));
[frag_curves.collapse_drift.theta, frag_curves.collapse_drift.beta] = fn_fit_fragility_MOM(ida_table.sa_x(ida_table.collapse == 1));
[frag_curves.collapse_convergence.theta, frag_curves.collapse_convergence.beta] = fn_fit_fragility_MOM(ida_table.sa_x(ida_table.collapse == 3));
[frag_curves.UR_accept_15.theta, frag_curves.UR_accept_15.beta] = fn_fit_fragility_MOM(ida_table.sa_x(ida_table.cols_walls_1_num_b_e_15 > 0));
[frag_curves.UR.theta, frag_curves.UR.beta] = fn_fit_fragility_MOM(ida_table.sa_x(ida_table.collapse == 1 | ida_table.collapse == 3 | ida_table.cols_walls_1_num_b_e_15 > 0));

% non directional component fragilities 
for p = 1:length(params)
    [frag_curves.cols_walls_1.(params{p})] = fn_multi_frag_curves(ida_table, 'cols_walls_1', params{p}, frag_probs, num_comps_cols_walls_1);
end

% 1% to 5% drift fragilities
for d = 1:5 
    [frag_curves.drift.(['idr_' num2str(d)]).theta, frag_curves.drift.(['idr_' num2str(d)]).beta] = fn_fit_fragility_MOM(ida_table.sa_x(ida_table.drift_x >= d/100 | ida_table.drift_z >= d/100));
end

% grav load lost fragilities
for f = 1:length(frag_probs) 
    [frag_curves.gravity.(['percent_lost_' num2str(frag_probs(f))]).theta, frag_curves.gravity.(['percent_lost_' num2str(frag_probs(f))]).beta] = fn_fit_fragility_MOM(ida_table.sa_x(ida_table.gravity_load_lost_ratio >= frag_probs(f)/100));
end

% adjacent components
[frag_curves.adjacent_comp.any.theta, frag_curves.adjacent_comp.any.beta] = fn_fit_fragility_MOM(ida_table.sa_x(ida_table.adjacent_failure_any == 1));
[frag_curves.adjacent_comp.any_frame.theta, frag_curves.adjacent_comp.any_frame.beta] = fn_fit_fragility_MOM(ida_table.sa_x(ida_table.adjacent_failure_any_frame == 1));
[frag_curves.adjacent_comp.all.theta, frag_curves.adjacent_comp.all.beta] = fn_fit_fragility_MOM(ida_table.sa_x(ida_table.adjacent_failure_all == 1));

% X direction Curves
[frag_curves.ew_collapse.theta, frag_curves.ew_collapse.beta] = fn_fit_fragility_MOM(ida_table.sa_x(strcmp(ida_table.collapse_direction,'x')));
for p = 1:length(params)
    [frag_curves.cols_1.(params{p})] = fn_multi_frag_curves(ida_table, 'cols_1', params{p}, frag_probs, num_comps_cols_1);
end

% Z direction curves
if analysis.run_z_motion
    [frag_curves.ns_collapse.theta, frag_curves.ns_collapse.beta] = fn_fit_fragility_MOM(ida_table.sa_x(strcmp(ida_table.collapse_direction,'z')));
    for p = 1:length(params)
        [frag_curves.walls_1.(params{p})] = fn_multi_frag_curves(ida_table, 'walls_1', params{p}, frag_probs, num_comps_walls_1);
    end
end

% Save Frag Curve Data
save([write_dir filesep 'frag_curves.mat'],'frag_curves')

end

function [theta, beta] = fn_fit_fragility_MOM(limit_state_dist)
[pHat, ~] = lognfit(limit_state_dist);
theta = exp(pHat(1));
beta = pHat(2);
end

function [frag_curves] = fn_multi_frag_curves(ida_table, mech, param, frag_probs, num_comp_mech)
frag_curves = table;
frag_curves.num_comp(1) = 1;
frag_curves.prct_mech(1) = round(1/num_comp_mech,3);
[frag_curves.theta(1), frag_curves.beta(1)] = fn_fit_fragility_MOM(ida_table.sa_x(ida_table.([mech '_num_' param]) > 0));
for pr = 1:length(frag_probs)
    frag_curves.num_comp(pr+1) = ceil(num_comp_mech*frag_probs(pr)/100);
    frag_curves.prct_mech(pr+1) = frag_probs(pr)/100;
    [frag_curves.theta(pr+1), frag_curves.beta(pr+1)] = fn_fit_fragility_MOM(ida_table.sa_x(ida_table.([mech '_percent_' param]) >= frag_probs(pr)/100));
end
end

function [ num_eles, percent_eles, num_eles_15, max_ele_ratio, mean_ele_ratio ] = fn_collect_ida_data(var_name, collapse_flag, collaspe_dir, ele_hinges, ele_hinges_alt)

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


max_ele_ratio = max([ele_ratios;ele_ratios_alt]);
mean_ele_ratio = mean([ele_ratios;ele_ratios_alt]);
end