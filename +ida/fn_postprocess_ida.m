function [ ] = fn_postprocess_ida(analysis, model, story, element, node, hinge, joint, ground_motion, ida_opensees_dir, ida_summary_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% import packages
import opensees.main_post_process_opensees
import opensees.post_process.*
import ida.*
import plotting_tools.fn_curt_plot_2D

% Load element properties table
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

%% Post process IDA using opensees post processor
main_post_process_opensees( analysis, model, story, node, element, joint, hinge, ground_motion, ida_opensees_dir )

%% Create IDA Summary DATA
% Load existing summary data
load([ida_summary_dir filesep 'summary_results.mat'])

% Max Story Drifts
load([ida_opensees_dir filesep 'story_analysis.mat']);
summary.max_drift_x = max(story.max_drift_x);
if analysis.run_z_motion
    summary.max_drift_z = max(story.max_drift_z);
end

% Hinge Ratios and Acceptance Criteria
if analysis.nonlinear ~= 0 && ~isempty(hinge)
    load([ida_opensees_dir filesep 'hinge_analysis.mat'])
    for i = 1:height(hinge)
        ele = element(element.id == hinge.element_id(i),:);
        if analysis.model_type == 3
            ele_prop = ele;
        else
            ele_prop = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
        end
        if exist([ida_opensees_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat'],'file')
            load([ida_opensees_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat']);
            hinge.hinge_deform(i) = max(abs(hin_TH.deformation_TH));
            hinge.ele_type(i) = ele_prop.type;

            if strcmp(hinge.ele_direction{i},'z') % in plane Walls
                hinge.tot_deform(i) = hinge.hinge_deform(i);
                hinge.d_value_tot(i) = ele.(['d_hinge_' num2str(hinge.ele_side(i))]);
                hinge.d_ratio(i) = hinge.tot_deform(i) / hinge.d_value_tot(i);
                hinge.e_value_tot(i) = ele.(['e_hinge_' num2str(hinge.ele_side(i))]);
                hinge.e_ratio(i) = hinge.tot_deform(i) / hinge.e_value_tot(i);
                
                if analysis.model_type ~= 3 % not archetype models
                    hinge.io_value_tot(i) = ele.(['io_' num2str(hinge.ele_side(i))]);
                    hinge.io_ratio(i) = hinge.tot_deform(i) / hinge.io_value_tot(i);
                    hinge.ls_value_tot(i) = ele.(['ls_' num2str(hinge.ele_side(i))]);
                    hinge.ls_ratio(i) = hinge.tot_deform(i) / hinge.ls_value_tot(i);
                    hinge.cp_value_tot(i) = ele.(['cp_' num2str(hinge.ele_side(i))]);
                    hinge.cp_ratio(i) = hinge.tot_deform(i) / hinge.cp_value_tot(i);
                end
            else % Everything Else
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
                hinge.a_value_tot(i) = ele.(['a_hinge_' num2str(hinge.ele_side(i))]) + theta_yeild_total;
                hinge.a_ratio(i) = hinge.tot_deform(i) / hinge.a_value_tot(i);
                hinge.b_value_tot(i) = ele.(['b_hinge_' num2str(hinge.ele_side(i))]) + theta_yeild_total;
                hinge.b_ratio(i) = hinge.tot_deform(i) / hinge.b_value_tot(i);

                if analysis.model_type ~= 3 % not archetype models
                    hinge.io_value_tot(i) = ele.(['io_' num2str(hinge.ele_side(i))]) + theta_yeild_total;
                    hinge.io_ratio(i) = hinge.tot_deform(i) / hinge.io_value_tot(i);
                    hinge.ls_value_tot(i) = ele.(['ls_' num2str(hinge.ele_side(i))]) + theta_yeild_total;
                    hinge.ls_ratio(i) = hinge.tot_deform(i) / hinge.ls_value_tot(i);
                    hinge.cp_value_tot(i) = ele.(['cp_' num2str(hinge.ele_side(i))]) + theta_yeild_total;
                    hinge.cp_ratio(i) = hinge.tot_deform(i) / hinge.cp_value_tot(i);

                    % Shear force Properties
                    hinge.shear_demand(i) = ele.(['Vmax_' num2str(hinge.ele_side(i))]);
                    hinge.asce41_shear_capacity(i) = ele.(['Vn_' num2str(hinge.ele_side(i))]);

                    if strcmp(ele.type,'column')
                        % Calculate the column axial capacity based on Elwood 2004
                        % limit state (only for flexure-shear columns)
                        as = ele_prop.(['Av_' num2str(hinge.ele_side(i))]);
                        fyt = ele_prop.fy_e;
                        dc = ele_prop.w/2 - ele_prop.clear_cover;
                        theta = 65*pi/180; % assumed 65 deg from Elwood 2004 
                        s = ele_prop.(['S_' num2str(hinge.ele_side(i))]);
                        hinge.P_demand(i) = ele.Pmax;
                        hinge.P_capacity(i) = min([(as*fyt*dc*tan(theta)/s)*((4*(1+(tan(theta)^2))/(100*summary.max_drift_x)) - tan(theta)),ele.Pn_c]);
                        hinge.P_dcr(i) = hinge.P_demand(i) / hinge.P_capacity(i);
                    else % for beams, ignore axial stuff
                        hinge.P_demand(i) = 0;
                        hinge.P_capacity(i) = 0;
                        hinge.P_dcr(i) = 0;
                    end
                end

    %             % EuroCode
    %             h = ele_prop.h;           % in
    %             b = ele_prop.w;           % in
    %             d = ele_prop.d_eff;       % in
    %             d_prm = ele_prop.h - ele_prop.d_eff; % in
    %             s = ele_prop.(['S_' num2str(hinge.ele_side(i))]);      % in
    %             Av = ele_prop.(['Av_' num2str(hinge.ele_side(i))]);    % in2
    %             As = sum(str2double(strsplit(strrep(strrep(ele_prop.As{1},']',''),'[',''),',')));  % in2
    %             As_comp = 0.5*sum(str2double(strsplit(strrep(strrep(ele_prop.As{1},']',''),'[',''),',')));  % in2
    %             As_ten = 0.5*sum(str2double(strsplit(strrep(strrep(ele_prop.As{1},']',''),'[',''),',')));  % in2
    %             db = mean(str2double(strsplit(strrep(strrep(ele_prop.d_b{1},']',''),'[',''),',')));  % in2
    %             bi = [10,10,10,10,10,10,10,10];
    %             fc = ele_prop.fc_e;    % psi
    %             fy = ele_prop.fy_e;   % psi
    %             P = ele.Pmax;    % lbs
    %             M = ele.(['Mmax_' num2str(hinge.ele_side(i))]);    % lbs-in
    %             V = ele.(['Vmax_' num2str(hinge.ele_side(i))]);    % lbs
    %             av = 1; % av = 1 if shear cracking expected before flexure, av = 0 othewise 
    %             d_NA = h/2; % assume nuetral axis is half of the member depth
    %             phi_y = 0.003 / d_NA;   % Yeild Curvature of the end section (in/in^2) (assume is same as concrete crushing over half of the member depth)
    %             deform_pl = hinge.plastic_deform(i);
    %             cov = 1.5 + 0.25; % Concrete cover plus half of the tie (assuming #4)
    %             [hinge.euro_V_NC(i)] = fn_eurocode_column_shear_acceptance(h,b,d,d_prm,s,As,Av,db,fc,fy,P,M,V,deform_pl,phi_y,d_NA); % assumes colums are shear controlled
    %         %                 [~, hinge.euro_th_NC_value(i)] = parametric_study(L,h,b,s,Av,cov,fc,fy,P); 
    %             [hinge.euro_th_NC_value(i), hinge.euro_th_SD_value(i), hinge.euro_th_DL_value(i)] = fn_eurocode_rotation_acceptance(h,b,d,d_prm,cov,s,As_comp,As_ten,Av,db,bi,fc,fy,P,M,V,av,phi_y);
    %             hinge.euro_th_NC_ratio(i) = hinge.tot_deform(i) / hinge.euro_th_NC_value(i);
    %             hinge.euro_th_SD_ratio(i) = hinge.tot_deform(i) / hinge.euro_th_SD_value(i);
    %             hinge.euro_th_DL_ratio(i) = hinge.tot_deform(i) / hinge.euro_th_DL_value(i);
            end
        end
    end


    % Convergence Collapse Check
    if isfield(hinge,'e_ratio')
        max_ratio = max(max(hinge.b_ratio,hinge.e_ratio));
    else
        max_ratio = max(hinge.b_ratio);
    end
    if summary.collapse == 3 && max_ratio < 1.5
        summary.collapse = 5; % if no element goes past 1.5b, dont treat convergence as collapse
    end

    % Collapse Direction and mechanism
    if summary.collapse > 0
        if ~analysis.run_z_motion || (summary.max_drift_x > summary.max_drift_z)
            summary.collapse_direction = 'x';
            collapse_story = story.id(story.max_drift_x > 0.05);
            if length(collapse_story) == 1 % Assume single story means column mechanism
                summary.collaspe_mech = ['story ' num2str(collapse_story) ' columns'];
            elseif length(collapse_story) > 1 % Assume Multi-Story means beam mechanism
                summary.collaspe_mech = ['story ' num2str(collapse_story(1:end-1)') ' beams'];
            else
                summary.collaspe_mech = 'other';
            end
        else
            summary.collapse_direction = 'z';
            if min(hinge.e_ratio(strcmp(hinge.ele_type,'wall') & strcmp(hinge.ele_direction,'z') & strcmp(hinge.direction,'primary') & hinge.story == 1)) > 1
                summary.collaspe_mech = 'story 1 walls';
            else
                summary.collaspe_mech = 'other';
            end
        end
    end

    % Plot Collapse Mechanism 
%     fn_curt_plot_2D( hinge, element, node, story, 'Collapse Mechanism', ida_summary_dir )
end

% Save data for this run
fprintf('Writing IDA Results to Directory: %s \n',ida_summary_dir)
if ~analysis.simple_recorders 
    save([ida_summary_dir filesep 'hinge_analysis.mat'],'hinge')
end
save([ida_summary_dir filesep 'story_analysis.mat'],'story')
save([ida_summary_dir filesep 'summary_results.mat'],'summary')

% Delete post processed middle man opensees data
try
    rmdir(ida_opensees_dir, 's')
catch
    warning('Could not remove directory')
end
% delete([ida_opensees_dir filesep 'story_analysis.mat'])
% delete([ida_opensees_dir filesep 'node_analysis.mat'])
% delete([ida_opensees_dir filesep 'model_analysis.mat'])
% delete([ida_opensees_dir filesep 'joint_analysis.mat'])
% delete([ida_opensees_dir filesep 'hinge_analysis.mat'])
% delete([ida_opensees_dir filesep 'element_analysis.mat'])
% delete([ida_opensees_dir filesep 'node_TH_*.mat'])
% delete([ida_opensees_dir filesep 'hinge_TH_*.mat'])

end %function

