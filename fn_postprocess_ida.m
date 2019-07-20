function [ ] = fn_postprocess_ida(analysis, model, story, element, node, hinge, gm_set_table, gm_idx, scale_factor)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% import packages
import opensees.main_post_process_opensees
import opensees.post_process.*

% Load element properties table
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

% Defin gms for this run
ground_motion.x = gm_set_table(gm_idx,:);
ground_motion.x.eq_dir = {['ground_motions' '/' analysis.gm_set '/' ground_motion.x.eq_name{1}]};
ground_motion.x.eq_name = {[ground_motion.x.eq_name{1} '.tcl']};
ground_motion.z = gm_set_table(gm_set_table.set_id == ground_motion.x.set_id & gm_set_table.pair ~= ground_motion.x.pair,:);
ground_motion.z.eq_dir = {['ground_motions' '/' analysis.gm_set '/' ground_motion.z.eq_name{1}]};
ground_motion.z.eq_name = {[ground_motion.z.eq_name{1} '.tcl']};

% Create Directories
opensees_outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'IDA' '/' 'Scale_' num2str(scale_factor) '/' 'GM_' num2str(ground_motion.x.set_id) '_' num2str(ground_motion.x.pair)];
ida_outputs_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' num2str(analysis.id) '/' 'IDA' '/' 'Summary Data' '/' 'Scale_' num2str(scale_factor) '/' 'GM_' num2str(ground_motion.x.set_id) '_' num2str(ground_motion.x.pair)];

%% Post process IDA using opensees post processor
main_post_process_opensees( analysis, model, story, node, element, [], hinge, ground_motion, opensees_outputs_dir )

%% Create IDA Summary DATA
% Load existing summary data
load([ida_outputs_dir filesep 'summary_results.mat'])

% Max Story Drifts
load([opensees_outputs_dir filesep 'story_analysis.mat']);
summary.max_drift_x = max(story.max_drift_x);
summary.max_drift_z = max(story.max_drift_z);

% Hinge Ratios and Acceptance Criteria
load([opensees_outputs_dir filesep 'hinge_analysis.mat'])
for i = 1:height(hinge)
    ele = element(element.id == hinge.element_id(i),:);
    ele_prop = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
    if exist([opensees_outputs_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat'],'file')
        load([opensees_outputs_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat']);
        hinge.hinge_deform(i) = max(abs(hin_TH.deformation_TH));
        hinge.ele_type(i) = ele_prop.type;

        if strcmp(hinge.ele_direction{i},'z') % in plane Walls
            hinge.tot_deform(i) = hinge.hinge_deform(i);
            hinge.d_value_tot(i) = ele.(['d_hinge_' num2str(hinge.ele_side(i))]);
            hinge.d_ratio(i) = hinge.tot_deform(i) / hinge.d_value_tot(i);
            hinge.e_value_tot(i) = ele.(['e_hinge_' num2str(hinge.ele_side(i))]);
            hinge.e_ratio(i) = hinge.tot_deform(i) / hinge.e_value_tot(i);
        else % Everything Else
            % Shear force Properties
            hinge.shear_demand(i) = ele.(['Vmax_' num2str(hinge.ele_side(i))]);
            hinge.asce41_shear_capacity(i) = ele.(['Vn_' num2str(hinge.ele_side(i))]);

            % Rotation properties
            K_elastic = 6*ele_prop.e*ele_prop.iz/ele.length;
            theta_yeild_total = ele.Mn_pos_1/K_elastic; % only for column bases currently

            % Modify hinge rotation to be element rotation
            if hinge.hinge_deform(i) >= (1/11)*theta_yeild_total
                hinge.elastic_deform(i) = theta_yeild_total;
                hinge.tot_deform(i) = hinge.hinge_deform(i) + (10/11)*theta_yeild_total;
            else
                percent_yield = hinge.hinge_deform(i)/((1/11)*theta_yeild_total);
                hinge.elastic_deform(i) = percent_yield*theta_yeild_total;
                hinge.tot_deform(i) = hinge.elastic_deform(i);
            end

            hinge.plastic_deform(i) = hinge.tot_deform(i) - hinge.elastic_deform(i);

            hinge.b_value_tot(i) = ele.(['b_hinge_' num2str(hinge.ele_side(i))]) + theta_yeild_total;
            hinge.b_ratio(i) = hinge.tot_deform(i) / hinge.b_value_tot(i);
            hinge.io_value_tot(i) = ele.(['io_' num2str(hinge.ele_side(i))]) + theta_yeild_total;
            hinge.io_ratio(i) = hinge.tot_deform(i) / hinge.io_value_tot(i);
            hinge.ls_value_tot(i) = ele.(['ls_' num2str(hinge.ele_side(i))]) + theta_yeild_total;
            hinge.ls_ratio(i) = hinge.tot_deform(i) / hinge.ls_value_tot(i);
            hinge.cp_value_tot(i) = ele.(['cp_' num2str(hinge.ele_side(i))]) + theta_yeild_total;
            hinge.cp_ratio(i) = hinge.tot_deform(i) / hinge.cp_value_tot(i);

            % EuroCode
            h = ele_prop.h;           % in
            b = ele_prop.w;           % in
            d = ele_prop.d_eff;       % in
            d_prm = ele_prop.h - ele_prop.d_eff; % in
            s = ele_prop.(['S_' num2str(hinge.ele_side(i))]);      % in
            Av = ele_prop.(['Av_' num2str(hinge.ele_side(i))]);    % in2
            As = sum(str2double(strsplit(strrep(strrep(ele_prop.As{1},']',''),'[',''),',')));  % in2
            db = mean(str2double(strsplit(strrep(strrep(ele_prop.d_b{1},']',''),'[',''),',')));  % in2
            fc = ele_prop.fc_e;    % psi
            fy = ele_prop.fy_e;   % psi
            P = ele.Pmax;    % lbs
            M = ele.(['Mmax_' num2str(hinge.ele_side(i))]);    % lbs-in
            V = ele.(['Vmax_' num2str(hinge.ele_side(i))]);    % lbs
            deform_pl = hinge.plastic_deform(i);
            cov = 1.5 + 0.25; % Concrete cover plus half of the tie (assuming #4)
            [hinge.euro_V_NC(i)] = fn_eurocode_column_shear_acceptance(h,b,d,d_prm,s,As,Av,db,fc,fy,P,M,V,deform_pl); % assumes colums are shear controlled
        %                 [~, hinge.euro_th_NC_value(i)] = parametric_study(L,h,b,s,Av,cov,fc,fy,P); 
            [hinge.euro_th_NC_value(i), hinge.euro_th_SD_value(i), hinge.euro_th_DL_value(i)] = fn_eurocode_rotation_acceptance(h,b,d,d_prm,s,As,Av,db,cov,fc,fy,P,M,V);
            hinge.euro_th_NC_ratio(i) = hinge.tot_deform(i) / hinge.euro_th_NC_value(i);
            hinge.euro_th_SD_ratio(i) = hinge.tot_deform(i) / hinge.euro_th_SD_value(i);
            hinge.euro_th_DL_ratio(i) = hinge.tot_deform(i) / hinge.euro_th_DL_value(i);
        end
    end
end

% Convergence Collapse Check
if summary.collapse == 3 && max(max(hinge.b_ratio,hinge.e_ratio)) < 1.5
    summary.collapse = 5; % if no element goes past 1.5b, dont treat convergence as collapse
end

% Collapse Direction and mechanism
if summary.collapse > 0
    if summary.max_drift_x > summary.max_drift_z
        summary.collapse_direction = 'x';
        if min(hinge.b_ratio(strcmp(hinge.ele_type,'beam') & strcmp(hinge.ele_direction,'x'))) > 1
            summary.collaspe_mech = 'beams';
        elseif min(hinge.b_ratio(strcmp(hinge.ele_type,'column') & strcmp(hinge.ele_direction,'x') & hinge.story == 1)) > 1
            summary.collaspe_mech = 'story 1';
        elseif min(hinge.b_ratio(strcmp(hinge.ele_type,'column') & strcmp(hinge.ele_direction,'x') & hinge.story == 2)) > 1
            summary.collaspe_mech = 'story 2';
        elseif min(hinge.b_ratio(strcmp(hinge.ele_type,'column') & strcmp(hinge.ele_direction,'x') & hinge.story == 3)) > 1
            summary.collaspe_mech = 'story 3';
        elseif min(hinge.b_ratio(strcmp(hinge.ele_type,'column') & strcmp(hinge.ele_direction,'x') & hinge.story == 4)) > 1
            summary.collaspe_mech = 'story 4';
        elseif min(hinge.b_ratio(strcmp(hinge.ele_type,'column') & strcmp(hinge.ele_direction,'x') & hinge.story == 5)) > 1
            summary.collaspe_mech = 'story 5';
        elseif min(hinge.b_ratio(strcmp(hinge.ele_type,'column') & strcmp(hinge.ele_direction,'x') & hinge.story == 6)) > 1
            summary.collaspe_mech = 'story 6';
        else
            summary.collaspe_mech = 'other';
        end
    else
        summary.collapse_direction = 'z';
        if min(hinge.e_ratio(strcmp(hinge.ele_type,'wall') & strcmp(hinge.ele_direction,'z') & hinge.story == 1)) > 1
            summary.collaspe_mech = 'story 1 walls';
        else
            summary.collaspe_mech = 'other';
        end
    end
end

% Save data for this run
fprintf('Writing IDA Results to Directory: %s \n',ida_outputs_dir)
save([ida_outputs_dir filesep 'hinge_analysis.mat'],'hinge')
save([ida_outputs_dir filesep 'story_analysis.mat'],'story')
save([ida_outputs_dir filesep 'summary_results.mat'],'summary')

end %function

