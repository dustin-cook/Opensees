function [ ] = fn_postprocess_ida(analysis, model, story, element, node, hinge, gm_set_table, gm_idx, scale_factor)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Post Process Data
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
ida_outputs_file = [ida_outputs_dir filesep 'summary_results.mat'];

% Load summary data
if exist(ida_outputs_file,'file')
    load([ida_outputs_dir filesep 'summary_results.mat'])

    % Nodal displacements
    for n = 1:height(node)
       node_id = node.id(n);
       if node.record_disp(n)
           [ node_disp_raw ] = fn_xml_read([opensees_outputs_dir filesep 'nodal_disp_' num2str(node.id(n)) '.xml']);
           node_disp_raw = node_disp_raw'; % flip to be node per row
           disp_TH.(['node_' num2str(node_id) '_TH']).(['disp_x_TH']) = node_disp_raw(2,:);
           disp_TH.(['node_' num2str(node_id) '_TH']).(['disp_z_TH']) = node_disp_raw(3,:);  % Currently Hard coded to three dimensions
       end
    end

    % Hinge Deformations
    if analysis.nonlinear ~= 0 && ~isempty(hinge)
        for i = 1:height(hinge)
            hinge_y = node.y(node.id == hinge.node_1(i));
            if hinge_y == 0 && strcmp(hinge.direction{i},'primary') && strcmp(hinge.ele_direction{i},'x')
                hinge_id = element.id(end) + hinge.id(i);
                ele = element(element.id == hinge.element_id(i),:);
                ele_prop = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
                
                K_elastic = 6*ele_prop.e*ele_prop.iz/ele.length;
                theta_yeild_total = ele.Mn_pos_1/K_elastic; % only for column bases currently
                
                [ hinge_deformation_TH ] = fn_xml_read([opensees_outputs_dir filesep 'hinge_deformation_' num2str(hinge_id) '.xml']);
                
                hinge.hinge_deform(i) = max(abs(hinge_deformation_TH(:,2)));
                
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
                L = ele.length;
                deform_pl = hinge.plastic_deform(i);
                [hinge.euro_V_NC(i)] = fn_eurocode_column_shear_acceptance(h,b,d,d_prm,s,As,Av,db,fc,fy,P,M,V,deform_pl); % assumes colums are shear controlled
%                 [~, hinge.euro_th_NC_value(i)] = parametric_study(L,h,b,s,Av,cov,fc,fy,P); 
                [hinge.euro_th_NC_value(i)] = fn_eurocode_rotation_acceptance(L,h,b,s,Av,cov,fc,fy,P);
            end
        end
    end

    % Check bad convergence models to see if they are close enough to
    % collapse
    if summary.collapse == 3 && max(hinge.b_ratio(hinge.b_ratio ~= 0 & hinge.b_ratio ~= inf)) < 1.5
        summary.collapse = 5;
    end

    % Calc Story Drift
    [ story.max_drift_x ] = fn_drift_profile( disp_TH, story, node, 'x' );
    [ story.max_drift_z ] = fn_drift_profile( disp_TH, story, node, 'z' );

    summary.max_drift_x = max(story.max_drift_x);
    summary.max_drift_z = max(story.max_drift_z);
    if summary.collapse > 0
        if summary.max_drift_x > summary.max_drift_z
            summary.collapse_direction = 'x';
        else
            summary.collapse_direction = 'z';
        end
    end

    % Save data for this run
    fprintf('Writing IDA Results to Directory: %s \n',ida_outputs_dir)
    save([ida_outputs_dir filesep 'hinge_analysis.mat'],'hinge')
    save([ida_outputs_dir filesep 'story_analysis.mat'],'story')
    save([ida_outputs_dir filesep 'summary_results.mat'],'summary')
end
end %function

