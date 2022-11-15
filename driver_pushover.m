%% Clear the Workspace
clear
close
clc
fclose('all');

%% Description: Method to build an Opensees model and run a ASCE 41-17 teir 3 seismic assessment.

% Created By: Dustin Cook
% Date Created: 1/2/2019

% Inputs:

% Outputs:

% Assumptions:

%% User Inputs (Think about changing this to a file read and command line execution)
% analysis.model_id = 2;
analysis.model_type = 3; % 1 = SDOF, 2 = MDOF (default), 3 = Archetype model
analysis.proceedure = 'Pushover'; % LDP or NDP or test
analysis.id = 'Dissertation_Study_LP'; % ID of the analysis for it to create its own directory
analysis.summit = 0; % Write tcl files to be run on summit and change location of opensees call
analysis.skip_2_outputs = 0; % Skip all the way to the plotters

% dynamic analysis inputs
% analysis.gm_seq_id = 16; % Maybe also make this part ot the defaults or model?

% Archetpye Inputs
analysis.additional_elements = 1; % this is the leaning column
analysis.eq_lat_load_factor = 1;
analysis.nonlinear_type = 'lumped'; % lumped or fiber

%% Initial Setup
% Import packages
import asce_41.main_ASCE_41
import plotting_tools.*

% Define remote directory
remote_dir = ['G:\My Drive\Dissertation Archetype Study\Archetypes RCMF\Archetype Model Responses'];

%% Secondary Inputs
[ analysis ] = fn_analysis_options( analysis );

%% Pull Model Data
% Define models to run
model_data = readtable(['inputs' filesep 'archetype_models.csv'],'ReadVariableNames',true);
model_data = model_data(model_data.num_stories > 4,:);
model_data = model_data(model_data.ie == 1,:);
model_data = model_data(model_data.design_drift == 0.02,:);
num_models = height(model_data);

%% Initiate Analysis
data = table;
for m = 1:num_models % run for each model    
    % Load basic model data
    model = model_data(m,:);
    analysis.model_id = model.id;
    fprintf('Running Model %i of %i: %s\n', m, num_models, model.name{1})
    
    %% Run Pushover
    tic
    main_ASCE_41( analysis )
    toc
    
    %% Collect and save data
    % Define read and write directories
    main_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id];
    opensees_dir = [main_dir '/' 'opensees_data'];
    push_dir = [main_dir '/' 'pushover'];
    write_dir = [remote_dir filesep model.name{1}];
    if ~exist(write_dir,'dir')
        mkdir(write_dir)
    end
    
    % Load pushover parameters
    model_analysis = load([opensees_dir filesep 'model_analysis.mat']);
    load([opensees_dir filesep 'story_analysis.mat'])
    load([opensees_dir filesep 'node_analysis.mat']);
    load([push_dir filesep 'story_TH.mat']);
    load([push_dir filesep 'element_analysis.mat']);
    load([push_dir filesep 'hinge_analysis.mat']);

    % Grab pushover THs
    control_nodes = node(node.primary_story == 1,:);
    roof_node = control_nodes(control_nodes.y == max(control_nodes.y),:);
    load([push_dir filesep 'node_TH_' num2str(roof_node.id) '.mat'])
    roof_disp_x = abs(nd_TH.disp_x_TH);
%     roof_disp_x_neg = abs(nd_TH.disp_x_neg_TH);
    v_ratio_x = abs(story_TH.base_shear_x_TH)/sum(story.seismic_wt);
%     v_ratio_x_neg = abs(story_TH.base_shear_x_neg_TH)/sum(story.seismic_wt);
    
%     % Calculate p695 pushover parameters
%     % (Assume buildin is symmetric, so just use the x pos direction)
%     Vmax = max(v_ratio_x);
%     Vmax_80 = 0.8*Vmax;
%     delta_vmax = roof_disp_x(Vmax == v_ratio_x); 
%     post_peak_V = v_ratio_x(roof_disp_x > delta_vmax);
%     post_peak_delta = roof_disp_x(roof_disp_x > delta_vmax);
%     delta_u = post_peak_delta(find(post_peak_V < Vmax_80,1)); % Find the first time its less than 0.8Vmax post peak
%     C0 = story.mode_shape_x(end) * sum(story.seismic_wt .* story.mode_shape_x) / sum(story.seismic_wt .* story.mode_shape_x .^2);
%     g = 386.4; % in/s^2
%     Ct = 0.016; % table 12.8-2
%     x = 0.9; % table 12.8-2
%     hn = sum(story.story_ht)/12;
%     Sds = 1.0;
%     Cu = interp1([0.1, 0.15, 0.2, 0.3],... % Coefficient for upper limit on calculated period table 12.8-1
%                  [1.7,1.6,1.5,1.4],...
%                  min(max(Sds,0.1),0.3));
%     Ta = Ct*hn^x;
%     CuTa = Cu*Ta;
%     dy_eff = C0*Vmax * (g/(4*pi^2)) * max(model_analysis.model.T1_x,CuTa);
%     mu_t = delta_u / dy_eff;
%     
%     % Plot parameters for verification
%     hold on
%     plot(roof_disp_x,v_ratio_x,'DisplayName','Pushover')
%     scatter(delta_vmax,Vmax,'DisplayName','Vmax')
%     scatter(delta_u,Vmax_80,'DisplayName','delta_u')
%     plot([dy_eff,dy_eff],[0,Vmax],'--','DisplayName','dy_{eff}')
%     box on
%     grid on
%     legend('location','northeast')
%     saveas(gcf,[write_dir filesep 'pushover_points.png'])
%     close
%     
%     % collect data in table
%     data.model_id(m,1) = model.id;
%     data.period(m,1) = model_analysis.model.T1_x;
%     data.Vmax(m,1) = Vmax;
%     data.delta_u(m,1) = delta_u;
%     data.dy_eff(m,1) = dy_eff;
%     data.mu_t(m,1) = mu_t;
%     
%     % save individual model data
%     writetable(data(m,:),[write_dir filesep 'pushover_data.csv'])
    
    % Load in hinge data
    for i = 1:height(hinge)
        ele = element(element.id == hinge.element_id(i),:);
        if analysis.model_type == 3
            ele_prop = ele;
        else
            ele_prop = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
        end
        if exist([push_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat'],'file')
            load([push_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat']);
            hinge.hinge_deform(i) = max(abs(hin_TH.deformation_TH));
            hinge.ele_type(i) = ele_prop.type;

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
        end
    end
    
    if strcmp(analysis.nonlinear_type,'fiber')
        % Load element properties table
        ele_prop_table = readtable([model.design_sheet_dir{1} filesep model.design_sheet_name{1} '.xlsm'],'Sheet','element'); % for archetype models, the model properties are already in the table

        % Filter out pdelta column elements
        element(element.rigid == 1,:) = [];
        
        for i = 1:height(element)
            % Load element rotation data
            load([push_dir filesep 'element_TH_' num2str(element.id(i)) '.mat']);
            element.total_deform(i) = max(abs(ele_TH.rot));
            
            ele_prop = ele_prop_table(ele_prop_table.id == element.ele_id(i),:);
            
            K_elastic = 6*ele_prop.e*ele_prop.iz/ele_prop.length;
            theta_yeild_total = ele_prop.Mn_pos/K_elastic;
                
            element.plastic_deform(i) = max(element.total_deform(i) - theta_yeild_total,0);
        end
        
        % Filter to just top 5 components with largest plastic deformation
%         element_filt = sortrows(element,'plastic_deform','descend');
        fn_plot_hinge_response_ida_fiber( [push_dir filesep 'Pushover_Plots'], push_dir, element(element.story == 1,:), ele_prop_table, analysis )
        
        % write hinge table to summary directory
        writetable(element,[push_dir filesep 'element_data.csv'])
    else
        % Create some component deformation plots
        fn_curt_plot_2D( hinge, element, node, story, 'Collapse Mechanism', [push_dir filesep 'Pushover_Plots'])

        % Filter to just top 5 components with largest plastic deformation
%         hinge_filt = sortrows(hinge,'plastic_deform','descend');
        fn_plot_hinge_response_ida( [push_dir filesep 'Pushover_Plots'], push_dir, hinge(hinge.story == 1,:), element, element, node, analysis )

        % write hinge table to summary directory
        writetable(hinge,[push_dir filesep 'Pushover_Plots' filesep 'hinge_data.csv'])
    end

end

% Save data table
% writetable(data,[remote_dir filesep 'pushover_data.csv'])
