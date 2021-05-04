function [ ] = fn_create_fragilities(analysis, write_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import packages
import plotting_tools.*

% Defined fixed parames
% params = {'b','e','b_e','io','ls','cp','euro_th_NC','euro_th_SD','euro_th_DL'};
% mechs = { 'cols_1', 'walls_1', 'cols_walls_1'};
if analysis.run_z_motion 
%     attrs = {'sa_geo', 'max_drift', 'drift_x', 'drift_z', 'gravity_dcr', 'lat_cap_ratio_any', 'lat_cap_ratio_both', 'lat_cap_ratio_max', ...
%             'num_cp', 'percent_cp','max_cp','min_cp','mean_cp','range_cp','std_cp','cov_cp',... 
%             'num_b','percent_b','min_b','max_b','mean_b', 'num_full_col_fails',...
%             'num_e','percent_e','min_e','max_e','mean_e'};
    params = {'be','io','ls','cp','ls_43','be_15', 'be_std', 'be_half_std'};
%     params_first = {};
else
%     attrs = {'sa_geo', 'max_drift', 'drift_x', 'gravity_dcr', 'lat_cap_ratio_any', 'lat_cap_ratio_both', 'lat_cap_ratio_max', ...
%         'num_cp', 'percent_cp','max_cp','min_cp','mean_cp','range_cp','std_cp','cov_cp',... 
%         'num_b','percent_b','min_b','max_b','mean_b', 'num_full_col_fails'};
    params = {'b','io','ls','cp','ls_43','b_15', 'b_std', 'b_half_std'};
%     params_first = {'b_15', 'b_std', 'b_std_half'};
end
frag_probs = [10 25 50 75 100];

% col_vars = {'collapse', 'collapse_drift', 'collapse_grav', 'collapse_grav_or_lat', 'collapse_nonconvergence'};
col_vars = {'collapse_drift', 'collapse_grav', 'collapse_grav_or_lat'};

%% Load Data and pre-allocate variables
ida_table = readtable([write_dir filesep 'ida_table.csv']);
gm_table = readtable([write_dir filesep 'gm_table.csv']);

%% Collect info for each ground motion (new fragilities)
% gm_data.collapse = table;
% gm_data.collapse_2 = table;
% gm_data.gravity = table;
% gm_data.collapse_x = table;
% gm_data.collapse_2_x = table;
% if analysis.run_z_motion 
%     gm_data.collapse_z = table;
%     gm_data.collapse_2_z = table;
% end

% for gm = 1:height(gm_table)
%     gm_response = ida_table(strcmp(ida_table.eq_name,gm_table.eq_name(gm)),:);
%     if ~isempty(gm_response)
%         gm_response_no_col = ida_table(strcmp(ida_table.eq_name,gm_table.eq_name(gm)) & ida_table.collapse == 0,:);
%         gm_response_no_col.collapse(gm_response_no_col.sa_geo == max(gm_response_no_col.sa_geo)) = 1; % Flag just before collapse point
%         jbc = gm_response_no_col;
% 
%         % Just Before Collapse
%         gm_data.collapse.eq_name{gm} = gm_table.eq_name{gm};
%         gm_data.collapse_x.eq_name{gm} = gm_table.eq_name{gm};
%         if analysis.run_z_motion 
%             gm_data.collapse_z.eq_name{gm} = gm_table.eq_name{gm};
%         end
%         for a = 1:length(attrs)
%             gm_data.collapse.(attrs{a})(gm) = min(jbc.(attrs{a})(jbc.collapse == 1));
%             if strcmp(gm_table.collapse_dir{gm},'x')
%                 gm_data.collapse_x.(attrs{a})(gm) = min(jbc.(attrs{a})(jbc.collapse == 1));
%                 if analysis.run_z_motion 
%                     gm_data.collapse_z.(attrs{a})(gm) = NaN;
%                 end
%             elseif strcmp(gm_table.collapse_dir{gm},'z')
%                 gm_data.collapse_x.(attrs{a})(gm) = NaN;
%                 gm_data.collapse_z.(attrs{a})(gm) = min(jbc.(attrs{a})(jbc.collapse == 1));
%             end
%         end
%         
%         % Post Process Collapse
%         pp_collapse = jbc;
%         pp_collapse.collapse(pp_collapse.max_drift >= 0.06) = 1;
%         pp_collapse.collapse(pp_collapse.gravity_dcr >= 1) = 1;
%         gm_data.collapse_2.eq_name{gm} = gm_table.eq_name{gm};
%         gm_data.collapse_2_x.eq_name{gm} = gm_table.eq_name{gm};
%         if analysis.run_z_motion 
%             gm_data.collapse_2_z.eq_name{gm} = gm_table.eq_name{gm};
%         end
%         for a = 1:length(attrs)
%             gm_data.collapse_2.(attrs{a})(gm) = min(pp_collapse.(attrs{a})(pp_collapse.collapse == 1));
%             if strcmp(gm_table.collapse_dir{gm},'x')
%                 gm_data.collapse_2_x.(attrs{a})(gm) = min(pp_collapse.(attrs{a})(pp_collapse.collapse == 1));
%                 if analysis.run_z_motion 
%                     gm_data.collapse_2_z.(attrs{a})(gm) = NaN;
%                 end
%             elseif strcmp(gm_table.collapse_dir{gm},'z')
%                 gm_data.collapse_2_x.(attrs{a})(gm) = NaN;
%                 gm_data.collapse_2_z.(attrs{a})(gm) = min(pp_collapse.(attrs{a})(pp_collapse.collapse == 1));
%             end
%         end
% 
%     %     % Full Gravity Exceedance
%     %     gm_data.full_gravity.eq_name{gm} = gm_set_table.eq_name{gm};
%     %     for a = 1:length(attrs)
%     %         gm_data.full_gravity.(attrs{a})(gm) = min(jbc.(attrs{a})(jbc.gravity_dcr >= 1));
%     %     end
%     end
% end
% 
% % Save Tabular Results as CSVs
% save([write_dir filesep 'gm_data.mat'],'gm_data')
% 
% %% Create Fragility Curves based on Baker MLE (New Fragilities)
% % collape and unnacceptable response
% field_names = fieldnames(gm_data);
% for i = 1:length(field_names)
%     fld = field_names{i};
%     for j = 2:width(gm_data.(fld))
%         var_name = gm_data.(fld).Properties.VariableNames{j};
%         [new_frag_curves.(fld).(var_name)] = fn_fit_fragility_MOM(gm_data.(fld).(var_name));
%     end
% end
% 
% % Save Frag Curve Data
% save([write_dir filesep 'new_frag_curves.mat'],'new_frag_curves')

%% Collect info for each ground motion (traditional fragilities)
for gm = 1:height(gm_table)
    gm_response = ida_table(strcmp(ida_table.eq_name,gm_table.eq_name(gm)),:);
    
    % Global limit states
    gm_table.sa_collapse(gm) = max([min(gm_response.sa_geo(gm_response.collapse > 0)),NaN]);
    gm_table.sa_collapse_drift(gm) = max([min(gm_response.sa_geo(gm_response.collapse == 1)),NaN]);
%     gm_table.sa_collapse_drift_6(gm) = max([min(gm_response.sa_geo(gm_response.max_drift >= 0.06)),NaN]);
    gm_table.sa_collapse_grav(gm) = max([min(gm_response.sa_geo(gm_response.gravity_dcr >= 1)),NaN]);
    
%     gm_table.sa_collapse_2(gm) = max([min(gm_response.sa_geo(gm_response.collapse > 0 | gm_response.max_drift >= 0.06 | gm_response.gravity_dcr > 1)),NaN]);
    gm_table.sa_collapse_grav_or_lat(gm) = max([min(gm_response.sa_geo(gm_response.collapse == 1 | gm_response.gravity_dcr > 1)),NaN]);
    
%     gm_table.sa_collapse_x(gm) = max([min(gm_response.sa_geo(gm_response.collapse_x > 0)),NaN]);
%     gm_table.sa_collapse_2_x(gm) = max([min(gm_response.sa_geo(gm_response.collapse_x > 0 | gm_response.drift_x >= 0.06 | gm_response.gravity_dcr > 1)),NaN]);
%     if analysis.run_z_motion 
%         gm_table.sa_collapse_z(gm) = max([min(gm_response.sa_geo(gm_response.collapse_z > 0)),NaN]);
%         gm_table.sa_collapse_2_z(gm) = max([min(gm_response.sa_geo(gm_response.collapse_z > 0 | gm_response.drift_z >= 0.06 | gm_response.gravity_dcr > 1)),NaN]);
%     end

    gm_table.sa_collapse_nonconvergence(gm) = max([min(gm_response.sa_geo(gm_response.collapse == 2 | gm_response.collapse == 3)),NaN]);
    
%     gm_table.sa_UR(gm) = max([min(gm_response.sa_geo(gm_response.collapse == 1 | gm_response.collapse == 3 | gm_response.num_b_15 > 0)),NaN]);

%     % Local Limit States
%     % 1.5*b
%     gm_table.sa_b_15(gm) = max([min(gm_response.sa_geo(gm_response.num_b_15 > 0)),NaN]);
%     
%     % b +stdev
%     gm_table.sa_b_std(gm) = max([min(gm_response.sa_geo(gm_response.num_b_std > 0)),NaN]);
%     
%     % b +stdev/2
%     gm_table.sa_b_half_std(gm) = max([min(gm_response.sa_geo(gm_response.num_b_half_std > 0)),NaN]);
    
    % 1% to 10% drift fragilities
    for d = 1:10 
        ls_lab = ['drift_' num2str(d)];
        gm_table.(['sa_' ls_lab])(gm) = max([min(gm_response.sa_geo(gm_response.max_drift >= d/100)),NaN]);
        for c = 1:length(col_vars)
            gm_table.(['CIR3_' ls_lab '_' col_vars{c}])(gm) = ...
                gm_table.(['sa_' col_vars{c}])(gm) / gm_table.(['sa_' ls_lab])(gm);
        end
%         if analysis.run_z_motion
%             gm_table.(['sa_drift_x_' num2str(d)])(gm) = max([min(gm_response.sa_geo(gm_response.drift_x >= d/100)),NaN]);
%             gm_table.(['sa_drift_z_' num2str(d)])(gm) = max([min(gm_response.sa_geo(gm_response.drift_z >= d/100)),NaN]);
%         end
    end
    
    % Drift normalized by pushover deformation capacity
    for r = 1:10
        rat = r*20; % ratio of peak dynamic drift to pushover capacity (in %) 
        ls_lab = ['push_cap_' num2str(rat)];
        gm_table.(['sa_' ls_lab])(gm) = max([min(gm_response.sa_geo(gm_response.push_deform_ratio_max_x >= rat/100)),NaN]);
        if analysis.run_z_motion
            max_push_cap = max(gm_response.push_deform_ratio_max_x,gm_response.push_deform_ratio_max_z);
            gm_table.(['sa_' ls_lab])(gm) = max([min(gm_response.sa_geo(max_push_cap >= rat/100)),NaN]);
        end
        for c = 1:length(col_vars)
            gm_table.(['CIR3_' ls_lab '_' col_vars{c}])(gm) = ...
                gm_table.(['sa_' col_vars{c}])(gm) / gm_table.(['sa_' ls_lab])(gm);
        end
    end
    
    % FEMA P-2018 Story Rating
    for r = 1:9
        ls_lab = ['p2018_story_rat_' num2str(r)];
        gm_table.(['sa_' ls_lab])(gm) = max([min(gm_response.sa_geo(gm_response.story_rating >= r/10)),NaN]);
%         if analysis.run_z_motion
%             gm_table.(['p2018_story_rat_z_' num2str(rat)])(gm) = max([min(gm_response.sa_geo(gm_response.story_rating >= r/10)),NaN]);
%         end
        for c = 1:length(col_vars)
            gm_table.(['CIR3_' ls_lab '_' col_vars{c}])(gm) = ...
                gm_table.(['sa_' col_vars{c}])(gm) / gm_table.(['sa_' ls_lab])(gm);
        end
    end
    
%     % grav load dcr
%     for f = 1:10
%         gm_table.(['sa_gravity_dcr_' num2str(f*10)])(gm) = max([min(gm_response.sa_geo(gm_response.gravity_dcr >= f/10)),NaN]);
%     end
%     
%     % lateral cap remaining
%     for f = 1:9
%         gm_table.(['sa_lat_cap_any_' num2str(f*10)])(gm) = max([min(gm_response.sa_geo(gm_response.lat_cap_ratio_any < f/10)),NaN]);
%         gm_table.(['sa_lat_cap_both_' num2str(f*10)])(gm) = max([min(gm_response.sa_geo(gm_response.lat_cap_ratio_both < f/10)),NaN]);
%         gm_table.(['sa_lat_cap_max_' num2str(f*10)])(gm) = max([min(gm_response.sa_geo(gm_response.lat_cap_ratio_max < f/10)),NaN]);
%     end
    
%     % adjacent components
%     gm_set_table.sa_adjacent_component_any(gm) = max([min(gm_response.sa_geo(gm_response.adjacent_failure_any == 1)),NaN]);
%     gm_set_table.sa_adjacent_component_any_frame(gm) = max([min(gm_response.sa_geo(gm_response.adjacent_failure_any_frame == 1)),NaN]);
%     gm_set_table.sa_adjacent_component_all(gm) = max([min(gm_response.sa_geo(gm_response.adjacent_failure_all == 1)),NaN]);
%     
%     % Energy
%     gm_table.collapse_energy(gm) = max([min(gm_response.total_energy_ft_lbs(gm_response.collapse == 1)),NaN]);
%     for f = 1:length(frag_probs)
%         gm_table.(['sa_collapse_energy_percent_' num2str(frag_probs(f))])(gm) = max([min(gm_response.sa_geo(gm_response.total_energy_ft_lbs/gm_table.collapse_energy(gm) >= frag_probs(f)/100)),NaN]);
%     end
    
    % non directional component fragilities 
    for p = 1:length(params)
        ls_lab = ['first_' params{p}];
        gm_table.(['sa_' ls_lab])(gm) = max([min(gm_response.sa_geo(gm_response.(['num_' params{p}]) > 0)),NaN]);
        for c = 1:length(col_vars)
            gm_table.(['CIR3_' ls_lab '_' col_vars{c}])(gm) = ...
                gm_table.(['sa_' col_vars{c}])(gm) / gm_table.(['sa_' ls_lab])(gm);
        end
%         if any(strcmp(gm_response.collapse_direction,'x'))
%             gm_table.(['sa_first_' params_full{p} '_x'])(gm) = max([min(gm_response.sa_geo(gm_response.(['num_' params_full{p}]) > 0)),NaN]);
%             if analysis.run_z_motion 
%                 gm_table.(['sa_first_' params_full{p} '_z'])(gm) = NaN;
%             end
%         elseif any(strcmp(gm_response.collapse_direction,'z'))
%             gm_table.(['sa_first_' params_full{p} '_x'])(gm) = NaN;
%             gm_table.(['sa_first_' params_full{p} '_z'])(gm) = max([min(gm_response.sa_geo(gm_response.(['num_' params_full{p}]) > 0)),NaN]);
%         end
        for pr = 1:length(frag_probs)
            ls_lab = [num2str(frag_probs(pr)) '_percent_' params{p}];
            gm_table.(['sa_' ls_lab])(gm) = ...
                max([min(gm_response.sa_geo(gm_response.(['percent_' params{p}]) >= frag_probs(pr)/100)),NaN]);
                for c = 1:length(col_vars)
                    gm_table.(['CIR3_' ls_lab '_' col_vars{c}])(gm) = ...
                    gm_table.(['sa_' col_vars{c}])(gm) / gm_table.(['sa_' ls_lab])(gm);
                end
%             if any(strcmp(gm_response.collapse_direction,'x'))
%                 gm_table.(['sa_' num2str(frag_probs(pr)) '_percent_' params_full{p} '_x'])(gm) = max([min(gm_response.sa_geo(gm_response.(['percent_' params_full{p}]) >= frag_probs(pr)/100)),NaN]);
%                 if analysis.run_z_motion 
%                     gm_table.(['sa_' num2str(frag_probs(pr)) '_percent_' params_full{p} '_z'])(gm) = NaN;
%                 end
%             elseif any(strcmp(gm_response.collapse_direction,'z'))
%                 gm_table.(['sa_' num2str(frag_probs(pr)) '_percent_' params_full{p} '_x'])(gm) = NaN;
%                 gm_table.(['sa_' num2str(frag_probs(pr)) '_percent_' params_full{p} '_z'])(gm) = max([min(gm_response.sa_geo(gm_response.(['percent_' params_full{p}]) >= frag_probs(pr)/100)),NaN]);
%             end
        end
    end
    
%     % X direction Curves
%     gm_set_table.sa_collapse_x(gm) = max([min(gm_response.sa_geo(strcmp(gm_response.collapse_direction,'x'))),NaN]);
%     for p = 1:length(params)
%         if any(strcmp(gm_response.collapse_direction,'x'))
%             gm_set_table.(['sa_cols_1_first_' params{p}])(gm) = max([min(gm_response.sa_geo(gm_response.(['cols_1_num_' params{p}]) > 0)),NaN]);
%         else
%             gm_set_table.(['sa_cols_1_first_' params{p}])(gm) = NaN;
%         end
%         for pr = 1:length(frag_probs)
%             if any(strcmp(gm_response.collapse_direction,'x'))
%                 gm_set_table.(['sa_cols_1_' num2str(frag_probs(pr)) '_percent_' params{p}])(gm) = max([min(gm_response.sa_geo(gm_response.(['cols_1_percent_' params{p}]) >= frag_probs(pr)/100)),NaN]);
%             else
%                 gm_set_table.(['sa_cols_1_' num2str(frag_probs(pr)) '_percent_' params{p}])(gm) = NaN;
%             end
%         end
%     end
%     
% 
%     % Z direction curves
%     if analysis.run_z_motion
%         gm_set_table.sa_collapse_z(gm) = max([min(gm_response.sa_geo(strcmp(gm_response.collapse_direction,'z'))),NaN]);
%         for p = 1:length(params)
%             gm_set_table.(['sa_walls_1_first_' params{p}])(gm) = max([min(gm_response.sa_geo(gm_response.(['walls_1_num_' params{p}]) > 0)),NaN]);
%             for pr = 1:length(frag_probs)
%                 gm_set_table.(['sa_walls_1_' num2str(frag_probs(pr)) '_percent_' params{p}])(gm) = max([min(gm_response.sa_geo(gm_response.(['walls_1_percent_' params{p}]) >= frag_probs(pr)/100)),NaN]);
%             end
%         end
%     end
end

%% Create Fragility Curves based on Baker MLE (traditional Fragilities)
% Global limit states
for v = 1:length(col_vars)
    [stats] = fn_fit_fragility_MOM(gm_table.(['sa_' col_vars{v}]));
    frag_curves.global.limit_state{v,1} = col_vars{v};
    frag_curves.global.median(v,1) = stats.theta;
    frag_curves.global.beta(v,1) = stats.beta;
end

% Local Limit states
% c = 1;
% frag_curves.local.limit_state_param{c,1} = '1.5*b';
% frag_curves.local.limit_state_degree{c,1} = 'first';
% [stats] = fn_fit_fragility_MOM(gm_table.sa_b_15);
% frag_curves.local.median(c,1) = stats.theta;
% frag_curves.local.beta(c,1) = stats.beta;
% 
% c = c + 1;
% frag_curves.local.limit_state_param{c,1} = 'b+stdev';
% frag_curves.local.limit_state_degree{c,1} = 'first';
% [stats] = fn_fit_fragility_MOM(gm_table.sa_b_std);
% frag_curves.local.median(c,1) = stats.theta;
% frag_curves.local.beta(c,1) = stats.beta;
% 
% c = c + 1;
% frag_curves.local.limit_state_param{c,1} = 'b+stdev/2';
% frag_curves.local.limit_state_degree{c,1} = 'first';
% [stats] = fn_fit_fragility_MOM(gm_table.sa_b_half_std);
% frag_curves.local.median(c,1) = stats.theta;
% frag_curves.local.beta(c,1) = stats.beta;

discrete_sa_vals = 0.2:0.05:3;
pdf_x_points = 0.01:0.01:4;

% non directional component fragilities 
c = 0;
for p = 1:length(params)
    c = c + 1;
    frag_curves.local.limit_state_param{c,1} = params{p};
    frag_curves.local.limit_state_degree{c,1} = 'first';
    ls_lab = ['first_' params{p}];
    frag_curves.local.limit_state_lab{c,1} = ls_lab;
    [stats] = fn_fit_fragility_MOM(gm_table.(['sa_' ls_lab]));
    frag_curves.local.median(c,1) = stats.theta;
    frag_curves.local.beta(c,1) = stats.beta;
    
    for v = 1:length(col_vars)
        CIR2b = [];
        for sa = 1:length(discrete_sa_vals)
            num_col = sum(gm_table.(['sa_' col_vars{v}]) <= discrete_sa_vals(sa));
            num_ls = sum(gm_table.(['sa_' 'first_' params{p}]) <= discrete_sa_vals(sa));
            CIR2b(sa) = num_col/num_ls;
        end
        frag_curves.local.sa_vals{c,1} = discrete_sa_vals;
        frag_curves.local.(['CIR2b_' col_vars{v}]){c,1} = CIR2b;
        
        % Reliability Index
        loc_pdf = lognpdf(pdf_x_points,log(stats.theta),stats.beta);
        filt = strcmp(frag_curves.global.limit_state,col_vars{v});
        glo_pdf = lognpdf(pdf_x_points,log(frag_curves.global.median(filt)),frag_curves.global.beta(filt));
        if frag_curves.global.median(filt) <= stats.theta
            r_index = 1;
        else
            combo_pdf = min(loc_pdf,glo_pdf);
            r_index = trapz(pdf_x_points,combo_pdf);
        end
        frag_curves.local.(['RI_' col_vars{v}])(c,1) = r_index;
            
        % P[Col|LS] dist
        filt = strcmp(frag_curves.global.limit_state,col_vars{v});
        col_med = frag_curves.global.median(filt);
        col_beta = frag_curves.global.beta(filt);
        p_col_giv_ls = logncdf(gm_table.(['sa_' ls_lab]),log(col_med),col_beta);
        gm_table.(['CIR2c_' ls_lab '_' col_vars{v}]) = p_col_giv_ls;
    end

    for pr = 1:length(frag_probs)
        c = c + 1;
        frag_curves.local.limit_state_param{c,1} = params{p};
        frag_curves.local.limit_state_degree{c,1} = [num2str(frag_probs(pr)) '%'];
        ls_lab = [num2str(frag_probs(pr)) '_percent_' params{p}];
        frag_curves.local.limit_state_lab{c,1} = ls_lab;
        [stats] = fn_fit_fragility_MOM(gm_table.(['sa_' ls_lab]));
        frag_curves.local.median(c,1) = stats.theta;
        frag_curves.local.beta(c,1) = stats.beta;
        
        for v = 1:length(col_vars)
            CIR2b = [];
            for sa = 1:length(discrete_sa_vals)
                num_col = sum(gm_table.(['sa_' col_vars{v}]) <= discrete_sa_vals(sa));
                num_ls = sum(gm_table.(['sa_' num2str(frag_probs(pr)) '_percent_' params{p}]) <= discrete_sa_vals(sa));
                CIR2b(sa) = num_col/num_ls;
            end
            frag_curves.local.sa_vals{c,1} = discrete_sa_vals;
            frag_curves.local.(['CIR2b_' col_vars{v}]){c,1} = CIR2b;
            
            % Reliability Index
            loc_pdf = lognpdf(pdf_x_points,log(stats.theta),stats.beta);
            filt = strcmp(frag_curves.global.limit_state,col_vars{v});
            glo_pdf = lognpdf(pdf_x_points,log(frag_curves.global.median(filt)),frag_curves.global.beta(filt));
            if frag_curves.global.median(filt) <= stats.theta
                r_index = 1;
            else
                combo_pdf = min(loc_pdf,glo_pdf);
                r_index = trapz(pdf_x_points,combo_pdf);
            end
            frag_curves.local.(['RI_' col_vars{v}])(c,1) = r_index;
            
            % P[Col|LS] dist
            filt = strcmp(frag_curves.global.limit_state,col_vars{v});
            col_med = frag_curves.global.median(filt);
            col_beta = frag_curves.global.beta(filt);
            p_col_giv_ls = logncdf(gm_table.(['sa_' ls_lab]),log(col_med),col_beta);
            gm_table.(['CIR2c_' ls_lab '_' col_vars{v}]) = p_col_giv_ls;
        end
    
    end
%     if analysis.run_z_motion 
%         [frag_curves.([params{p} '_x'])] = fn_multi_frag_curves(gm_table, params{p}, '_x', frag_probs, ida_table.num_comps(1));
%         [frag_curves.([params{p} '_z'])] = fn_multi_frag_curves(gm_table, params{p}, '_z', frag_probs, ida_table.num_comps(1));
%     end
end

% 1% to 10% drift fragilities
for d = 1:10 
    c = c + 1;
    frag_curves.local.limit_state_param{c,1} = 'peak drift';
    frag_curves.local.limit_state_degree{c,1} = [num2str(d) '%'];
    ls_lab = ['drift_' num2str(d)];
    frag_curves.local.limit_state_lab{c,1} = ls_lab;
    [stats] = fn_fit_fragility_MOM(gm_table.(['sa_' ls_lab]));
    frag_curves.local.median(c,1) = stats.theta;
    frag_curves.local.beta(c,1) = stats.beta;
        
    for v = 1:length(col_vars)
        CIR2b = [];
        for sa = 1:length(discrete_sa_vals)
            num_col = sum(gm_table.(['sa_' col_vars{v}]) <= discrete_sa_vals(sa));
            num_ls = sum(gm_table.(['sa_drift_' num2str(d)]) <= discrete_sa_vals(sa));
            CIR2b(sa) = num_col/num_ls;
        end
        frag_curves.local.sa_vals{c,1} = discrete_sa_vals;
        frag_curves.local.(['CIR2b_' col_vars{v}]){c,1} = CIR2b;
        
        % Reliability Index
        loc_pdf = lognpdf(pdf_x_points,log(stats.theta),stats.beta);
        filt = strcmp(frag_curves.global.limit_state,col_vars{v});
        glo_pdf = lognpdf(pdf_x_points,log(frag_curves.global.median(filt)),frag_curves.global.beta(filt));
        if frag_curves.global.median(filt) <= stats.theta
            r_index = 1;
        else
            combo_pdf = min(loc_pdf,glo_pdf);
            r_index = trapz(pdf_x_points,combo_pdf);
        end
        frag_curves.local.(['RI_' col_vars{v}])(c,1) = r_index;
        
        % P[Col|LS] dist
        filt = strcmp(frag_curves.global.limit_state,col_vars{v});
        col_med = frag_curves.global.median(filt);
        col_beta = frag_curves.global.beta(filt);
        p_col_giv_ls = logncdf(gm_table.(['sa_' ls_lab]),log(col_med),col_beta);
        gm_table.(['CIR2c_' ls_lab '_' col_vars{v}]) = p_col_giv_ls;
    end
    
%     if analysis.run_z_motion 
%         [frag_curves.drift.(['idr_x_' num2str(d)])] = fn_fit_fragility_MOM(gm_table.(['sa_drift_x_' num2str(d)]));
%         [frag_curves.drift.(['idr_z_' num2str(d)])] = fn_fit_fragility_MOM(gm_table.(['sa_drift_z_' num2str(d)]));
%     end
end

% pushover deform capacity fragilities
for r = 1:10
    rat = r*20;
    c = c + 1;
    frag_curves.local.limit_state_param{c,1} = 'deform_cap';
    frag_curves.local.limit_state_degree{c,1} = [num2str(rat) '%'];
    ls_lab = ['push_cap_' num2str(rat)];
    frag_curves.local.limit_state_lab{c,1} = ls_lab;
    [stats] = fn_fit_fragility_MOM(gm_table.(['sa_' ls_lab]));
    frag_curves.local.median(c,1) = stats.theta;
    frag_curves.local.beta(c,1) = stats.beta;
        
    for v = 1:length(col_vars)
        CIR2b = [];
        for sa = 1:length(discrete_sa_vals)
            num_col = sum(gm_table.(['sa_' col_vars{v}]) <= discrete_sa_vals(sa));
            num_ls = sum(gm_table.(['sa_push_cap_' num2str(rat)]) <= discrete_sa_vals(sa));
            CIR2b(sa) = num_col/num_ls;
        end
        frag_curves.local.sa_vals{c,1} = discrete_sa_vals;
        frag_curves.local.(['CIR2b_' col_vars{v}]){c,1} = CIR2b;
        
        % Reliability Index
        loc_pdf = lognpdf(pdf_x_points,log(stats.theta),stats.beta);
        filt = strcmp(frag_curves.global.limit_state,col_vars{v});
        glo_pdf = lognpdf(pdf_x_points,log(frag_curves.global.median(filt)),frag_curves.global.beta(filt));
        if frag_curves.global.median(filt) <= stats.theta
            r_index = 1;
        else
            combo_pdf = min(loc_pdf,glo_pdf);
            r_index = trapz(pdf_x_points,combo_pdf);
        end
        frag_curves.local.(['RI_' col_vars{v}])(c,1) = r_index;
        
        % P[Col|LS] dist
        filt = strcmp(frag_curves.global.limit_state,col_vars{v});
        col_med = frag_curves.global.median(filt);
        col_beta = frag_curves.global.beta(filt);
        p_col_giv_ls = logncdf(gm_table.(['sa_' ls_lab]),log(col_med),col_beta);
        gm_table.(['CIR2c_' ls_lab '_' col_vars{v}]) = p_col_giv_ls;
    end
    
%     if analysis.run_z_motion 
%         [frag_curves.drift.(['idr_x_' num2str(d)])] = fn_fit_fragility_MOM(gm_table.(['sa_drift_x_' num2str(d)]));
%         [frag_curves.drift.(['idr_z_' num2str(d)])] = fn_fit_fragility_MOM(gm_table.(['sa_drift_z_' num2str(d)]));
%     end
end

% FEMA P-2018 Story rating
for r = 1:9
    c = c + 1;
    frag_curves.local.limit_state_param{c,1} = 'p2018_rating';
    frag_curves.local.limit_state_degree{c,1} = [num2str(r/10)];
    ls_lab = ['p2018_story_rat_' num2str(r)];
    frag_curves.local.limit_state_lab{c,1} = ls_lab;
    [stats] = fn_fit_fragility_MOM(gm_table.(['sa_' ls_lab]));
    frag_curves.local.median(c,1) = stats.theta;
    frag_curves.local.beta(c,1) = stats.beta;
    
    for v = 1:length(col_vars)
        CIR2b = [];
        for sa = 1:length(discrete_sa_vals)
            num_col = sum(gm_table.(['sa_' col_vars{v}]) <= discrete_sa_vals(sa));
            num_ls = sum(gm_table.(['sa_p2018_story_rat_' num2str(r)]) <= discrete_sa_vals(sa));
            CIR2b(sa) = num_col/num_ls;
        end
        frag_curves.local.sa_vals{c,1} = discrete_sa_vals;
        frag_curves.local.(['CIR2b_' col_vars{v}]){c,1} = CIR2b;
        
        % Reliability Index
        loc_pdf = lognpdf(pdf_x_points,log(stats.theta),stats.beta);
        filt = strcmp(frag_curves.global.limit_state,col_vars{v});
        glo_pdf = lognpdf(pdf_x_points,log(frag_curves.global.median(filt)),frag_curves.global.beta(filt));
        if frag_curves.global.median(filt) <= stats.theta
            r_index = 1;
        else
            combo_pdf = min(loc_pdf,glo_pdf);
            r_index = trapz(pdf_x_points,combo_pdf);
        end
        frag_curves.local.(['RI_' col_vars{v}])(c,1) = r_index;
        
        % P[Col|LS] dist
        filt = strcmp(frag_curves.global.limit_state,col_vars{v});
        col_med = frag_curves.global.median(filt);
        col_beta = frag_curves.global.beta(filt);
        p_col_giv_ls = logncdf(gm_table.(['sa_' ls_lab]),log(col_med),col_beta);
        gm_table.(['CIR2c_' ls_lab '_' col_vars{v}]) = p_col_giv_ls;
    end
end
% hold on
% sorted_val = sort(gm_table.(['p2018_story_rat_' num2str(r)]));
% rank = (1:length(sorted_val)) / length(sorted_val);
% scatter(sorted_val,rank)
% x_vals = 0.01:0.01:1;
% cdf = logncdf(x_vals,log(stats.theta),stats.beta);
% plot(x_vals,cdf)
% % grav dcr
% for f = 1:10
%     [frag_curves.gravity.(['dcr_' num2str(f*10)])] = fn_fit_fragility_MOM(gm_table.(['sa_gravity_dcr_' num2str(f*10)]));
% end

% % lat cap
% for f = 1:9
%     [frag_curves.lateral.(['cap_any_' num2str(f*10)])] = fn_fit_fragility_MOM(gm_table.(['sa_lat_cap_any_' num2str(f*10)]));
%     [frag_curves.lateral.(['cap_both_' num2str(f*10)])] = fn_fit_fragility_MOM(gm_table.(['sa_lat_cap_both_' num2str(f*10)]));
%     [frag_curves.lateral.(['cap_max_' num2str(f*10)])] = fn_fit_fragility_MOM(gm_table.(['sa_lat_cap_max_' num2str(f*10)]));
% end

% % adjacent components
% [frag_curves.adjacent_comp.any] = fn_fit_fragility_MOM(gm_set_table.sa_adjacent_component_any);
% [frag_curves.adjacent_comp.any_frame] = fn_fit_fragility_MOM(gm_set_table.sa_adjacent_component_any_frame);
% [frag_curves.adjacent_comp.all] = fn_fit_fragility_MOM(gm_set_table.sa_adjacent_component_all);

% % Collapse Energy
% for f = 1:length(frag_probs) 
%     [frag_curves.energy.(['percent_collapse_' num2str(frag_probs(f))])] = fn_fit_fragility_MOM(gm_table.(['sa_collapse_energy_percent_' num2str(frag_probs(f))]));
% end

% % X direction Curves
% [frag_curves.ew_collapse] = fn_fit_fragility_MOM(gm_set_table.sa_collapse_x);
% for p = 1:length(params)
%     [frag_curves.cols_1.(params{p})] = fn_multi_frag_curves(gm_set_table, 'cols_1', params{p}, frag_probs, ida_table.num_comps_cols_1(1));
% end
% 
% % Z direction curves
% if analysis.run_z_motion
%     [frag_curves.ns_collapse] = fn_fit_fragility_MOM(gm_set_table.sa_collapse_z);
%     for p = 1:length(params)
%         [frag_curves.walls_1.(params{p})] = fn_multi_frag_curves(gm_set_table, 'walls_1', params{p}, frag_probs, ida_table.num_comps_walls_1(1));
%     end
% end

%% Save Data 
% Save Tabular GM Results as CSVs
writetable(gm_table,[write_dir filesep 'gm_table.csv'])

% Save Frag Curve Data
% writetable(struct2table(frag_curves.global),[write_dir filesep 'global_limit_states.csv'])
% writetable(struct2table(frag_curves.local),[write_dir filesep 'local_limit_states.csv'])
frag_curves.global = struct2table(frag_curves.global);
frag_curves.local = struct2table(frag_curves.local);
save([write_dir filesep 'frag_curves.mat'],'frag_curves')

% % Write outputs summary file
% outputs.med_sa_collapse = frag_curves.collapse.theta;
% outputs.med_sa_collapse_grav = frag_curves.collapse_2.theta;
% outputs.med_sa_cp = frag_curves.cp.theta(1);
% outputs.med_sa_cp_10 = frag_curves.cp.theta(frag_curves.cp.prct_mech == 0.1);
% outputs.med_sa_cp_25 = frag_curves.cp.theta(find(frag_curves.cp.prct_mech == 0.25,1,'first'));
% outputs.med_sa_cp_50 = frag_curves.cp.theta(find(frag_curves.cp.prct_mech == 0.5,1,'first'));
% outputs.med_sa_cp_75 = frag_curves.cp.theta(find(frag_curves.cp.prct_mech == 0.75,1,'first'));
% outputs.med_sa_b = frag_curves.b.theta(1);
% outputs.collapse_margin = outputs.med_sa_collapse/outputs.med_sa_cp;
% outputs.collapse_margin_grav = outputs.med_sa_collapse_grav/outputs.med_sa_cp;
% outputs.collapse_margin_grav_b = outputs.med_sa_collapse_grav/outputs.med_sa_b;
% % outputs.num_comps = median(gm_data.collapse_2.num_cp);
% % outputs.percent_comps = median(gm_data.collapse_2.percent_cp);
% % outputs.num_full_col_fails = median(gm_data.collapse_2.num_full_col_fails);
% % outputs.drift = new_frag_curves.collapse_2.drift_x.theta;
% outputs.med_sa_drift_2 = frag_curves.drift.idr_2.theta;
% outputs.med_sa_drift_3 = frag_curves.drift.idr_3.theta;
% outputs.med_sa_drift_4 = frag_curves.drift.idr_4.theta;
% % outputs.gravity_dcr = median(gm_data.collapse_2.gravity_dcr);
% % outputs.lat_cap_any = median(gm_data.collapse_2.lat_cap_ratio_any);
% % outputs.lat_cap_both = median(gm_data.collapse_2.lat_cap_ratio_both);
% % outputs.lat_cap_max = median(gm_data.collapse_2.lat_cap_ratio_max);
% % outputs.med_sa_lat_50 = frag_curves.lateral.cap_both_50.theta;
% % outputs.med_sa_lat_60 = frag_curves.lateral.cap_both_60.theta;
% % outputs.med_sa_lat_70 = frag_curves.lateral.cap_both_70.theta;
% % outputs.med_sa_lat_80 = frag_curves.lateral.cap_both_80.theta;
% % outputs.med_sa_lat_90 = frag_curves.lateral.cap_both_90.theta;
% % outputs.med_sa_grav_30 = frag_curves.gravity.dcr_30.theta;
% % outputs.med_sa_grav_40 = frag_curves.gravity.dcr_40.theta;
% % outputs.med_sa_grav_50 = frag_curves.gravity.dcr_50.theta;
% % outputs.med_sa_grav_60 = frag_curves.gravity.dcr_60.theta;
% % outputs.med_sa_grav_70 = frag_curves.gravity.dcr_70.theta;
% % outputs.med_sa_grav_80 = frag_curves.gravity.dcr_80.theta;
% % if analysis.run_z_motion
% %     outputs.med_sa_collapse_x = frag_curves.collapse_x.theta;
% %     outputs.med_sa_collapse_z = frag_curves.collapse_z.theta;
% %     outputs.med_sa_collapse_grav_x = frag_curves.collapse_2_x.theta;
% %     outputs.med_sa_collapse_grav_z = frag_curves.collapse_2_z.theta;
% %     outputs.med_sa_cp_x = frag_curves.cp_x.theta(1);
% %     outputs.med_sa_cp_z = frag_curves.cp_z.theta(1);
% % end
% writetable(struct2table(outputs),[write_dir filesep 'summary_outputs.csv'])
end

function [params] = fn_fit_fragility_MOM(limit_state_dist)
    limit_state_dist(isnan(limit_state_dist)) = [];
    limit_state_dist(limit_state_dist <= 0) = 1e-6;
    if ~isempty(limit_state_dist)
        [pHat, ~] = lognfit(limit_state_dist);
        params.theta = exp(pHat(1));
        params.beta = pHat(2);
    else
        params.theta = NaN;
        params.beta = NaN;
    end
end

function [frag_curves] = fn_multi_frag_curves(gm_set_table, param, dir, frag_probs, num_comp_mech)
    frag_curves = table;
    frag_curves.num_comp(1) = 1;
    frag_curves.prct_mech(1) = round(1/num_comp_mech,3);
    [fits] = fn_fit_fragility_MOM(gm_set_table.(['sa_' 'first_' param dir]));
    frag_curves.theta(1) = fits.theta;
    frag_curves.beta(1) = fits.beta;
    for pr = 1:length(frag_probs)
        frag_curves.num_comp(pr+1,1) = ceil(num_comp_mech*frag_probs(pr)/100);
        frag_curves.prct_mech(pr+1,1) = frag_probs(pr)/100;
        [fits] = fn_fit_fragility_MOM(gm_set_table.(['sa_' num2str(frag_probs(pr)) '_percent_' param dir]));
        frag_curves.theta(pr+1,1) = fits.theta;
        frag_curves.beta(pr+1,1) = fits.beta;
    end
end