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
analysis.id = 'Dissertation_Study'; % ID of the analysis for it to create its own directory
analysis.summit = 0; % Write tcl files to be run on summit and change location of opensees call
analysis.skip_2_outputs = 0; % Skip all the way to the plotters

% dynamic analysis inputs
% analysis.gm_seq_id = 16; % Maybe also make this part ot the defaults or model?

% Archetpye Inputs
analysis.additional_elements = 1; % this is the leaning column
analysis.eq_lat_load_factor = 1;

%% Initial Setup
% Import packages
import asce_41.main_ASCE_41

% Define remote directory
remote_dir = ['G:\My Drive\Dissertation Archetype Study\Archetypes RCMF\Archetype Model Responses'];

%% Secondary Inputs
[ analysis ] = fn_analysis_options( analysis );

%% Pull Model Data
% Define models to run
model_data = readtable(['inputs' filesep 'archetype_models.csv'],'ReadVariableNames',true);
model_data = model_data(model_data.num_stories == 12,:);
% model_data = model_data(model_data.num_stories ~= 20,:);
% model_data = model_data(model_data.ie ~= 1.25,:);
% model_data = model_data(~contains(model_data.name,'drift'),:);
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
    
    % Load pushover parameters
    model_analysis = load([opensees_dir filesep 'model_analysis.mat']);
    load([opensees_dir filesep 'story_analysis.mat'])
    load([opensees_dir filesep 'node_analysis.mat']);
    load([push_dir filesep 'story_TH.mat']);

    % Grab pushover THs
    control_nodes = node(node.primary_story == 1,:);
    roof_node = control_nodes(control_nodes.y == max(control_nodes.y),:);
    load([push_dir filesep 'node_TH_' num2str(roof_node.id) '.mat'])
    roof_disp_x = abs(nd_TH.disp_x_TH);
%     roof_disp_x_neg = abs(nd_TH.disp_x_neg_TH);
    v_ratio_x = abs(story_TH.base_shear_x_TH)/sum(story.seismic_wt);
%     v_ratio_x_neg = abs(story_TH.base_shear_x_neg_TH)/sum(story.seismic_wt);
    
    % Calculate p695 pushover parameters
    % (Assume buildin is symmetric, so just use the x pos direction)
    Vmax = max(v_ratio_x);
    Vmax_80 = 0.8*Vmax;
    delta_vmax = roof_disp_x(Vmax == v_ratio_x); 
    post_peak_V = v_ratio_x(roof_disp_x > delta_vmax);
    post_peak_delta = roof_disp_x(roof_disp_x > delta_vmax);
    delta_u = post_peak_delta(find(post_peak_V < Vmax_80,1)); % Find the first time its less than 0.8Vmax post peak
    C0 = story.mode_shape_x(end) * sum(story.seismic_wt .* story.mode_shape_x) / sum(story.seismic_wt .* story.mode_shape_x .^2);
    g = 386.4; % in/s^2
    Ct = 0.016; % table 12.8-2
    x = 0.9; % table 12.8-2
    hn = sum(story.story_ht)/12;
    Sds = 1.0;
    Cu = interp1([0.1, 0.15, 0.2, 0.3],... % Coefficient for upper limit on calculated period table 12.8-1
                 [1.7,1.6,1.5,1.4],...
                 min(max(Sds,0.1),0.3));
    Ta = Ct*hn^x;
    CuTa = Cu*Ta;
    dy_eff = C0*Vmax * (g/(4*pi^2)) * max(model_analysis.model.T1_x,CuTa);
    mu_t = delta_u / dy_eff;
    
%     % Plot parameters for verification
%     hold on
%     plot(roof_disp_x,v_ratio_x,'DisplayName','Pushover')
%     scatter(delta_vmax,Vmax,'DisplayName','Vmax')
%     scatter(delta_u,Vmax_80,'DisplayName','delta_u')
%     plot([dy_eff,dy_eff],[0,Vmax],'--','DisplayName','dy_eff')
%     close
    
    % collect data in table
    data.model_id(m,1) = model.id;
    data.period(m,1) = model_analysis.model.T1_x;
    data.Vmax(m,1) = Vmax;
    data.delta_u(m,1) = delta_u;
    data.dy_eff(m,1) = dy_eff;
    data.mu_t(m,1) = mu_t;
    
    % save individual model data
    write_dir = [remote_dir filesep model.name{1}];
    if ~exist(write_dir)
        mkdir(write_dir)
    end
    writetable(data(m,:),[write_dir filesep 'pushover_data.csv'])
end

% Save data table
writetable(data,[remote_dir filesep 'pushover_data.csv'])
