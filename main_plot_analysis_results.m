function [ ] = main_plot_analysis_results( model, analysis, ele_prop_table )
% Description: Main script that compiles interesting results and creates
% visuals.

% Created By: Dustin Cook
% Date Created: 1/3/2019

% Inputs:

% Outputs:

% Assumptions:


%% Initial Setup
% Import Packages
import plotting_tools.*
import asce_41.*

% Define Read and Write Directories
read_dir = [analysis.out_dir filesep 'asce_41_data'];
read_dir_opensees = [analysis.out_dir filesep 'opensees_data'];
plot_dir = [analysis.out_dir filesep 'analysis_plots'];
if analysis.filter_accel
    read_dir_TH = [analysis.out_dir filesep 'asce_41_data'];
else
    read_dir_TH = [analysis.out_dir filesep 'opensees_data'];
end

% Load Analysis Data
if exist(read_dir,'dir')
    load([read_dir filesep 'element_analysis.mat'])
    load([read_dir filesep 'story_analysis.mat'])
    load([read_dir filesep 'model_analysis.mat'])
    load([read_dir filesep 'node_analysis.mat']);
    load([read_dir filesep 'joint_analysis.mat'])
    if exist([read_dir filesep 'hinge_analysis.mat'],'file')
        load([read_dir filesep 'hinge_analysis.mat'])
    end
else
    load([read_dir_opensees filesep 'story_analysis.mat'])
end

if sum(analysis.type_list == 1) > 0 % Dynamic Analysis was run as part of this proceedure
    temp = load([read_dir_opensees filesep 'gm_data.mat']);
    ground_motion = temp.ground_motion;
    eq = temp.eq;
    eq_analysis_timespace = temp.eq_analysis_timespace;
    dirs_ran = temp.dirs_ran;
end

if exist([analysis.out_dir filesep 'backbones'],'dir')
    backbones = load([analysis.out_dir filesep 'backbones' filesep 'element_analysis.mat']);
end

if exist([analysis.out_dir filesep 'backbones_pushover'],'dir')
    backbones_pushover = load([analysis.out_dir filesep 'backbones_pushover' filesep 'element_analysis.mat']);
end

%% Cyclic Analysis
if sum(analysis.type_list == 3) > 0 % Cyclic was run as part of this proceedure
    if analysis.element_plots
        % Load Pushover Results
        cyclic_read_dir = [analysis.out_dir filesep 'cyclic'];
        cyclic_analysis = load([cyclic_read_dir filesep 'analysis_options.mat']);
        if cyclic_analysis.analysis.nonlinear ~= 0 && exist('backbones','var') % Is the pushover nonlinear?
            % Plot Hinge Response
            cyclic_hinge = load([cyclic_read_dir filesep 'hinge_analysis.mat']);
            fn_plot_hinge_response( cyclic_read_dir, cyclic_read_dir, cyclic_hinge.hinge, backbones_pushover.element, ele_prop_table, node, analysis, joint  )
        end
    end
end

%% Pushover Analysis
if sum(analysis.type_list == 2) > 0 % Pushover was run as part of this proceedure
    pushover_read_dir = [analysis.out_dir filesep 'pushover'];
    pushover_analysis = load([pushover_read_dir filesep 'analysis_options.mat']);
    pushover_story_TH = load([pushover_read_dir filesep 'story_TH.mat']);
    % Plot Building Pushovers
    [ model ] = fn_plot_pushover( pushover_read_dir, story, pushover_story_TH.story_TH, model );

    if analysis.element_plots
        % Plot PM Diagrams for each element
%         fn_plot_PM_response( pushover_read_dir, pushover_read_dir, element, analysis.hinge_stories_2_plot )

        % Plot Hinge Response
        if pushover_analysis.analysis.nonlinear ~= 0 && exist('backbones_pushover','var') % Is the pushover nonlinear?
            pushover_hinge = load([pushover_read_dir filesep 'hinge_analysis.mat']);
            fn_plot_hinge_response( pushover_read_dir, pushover_read_dir, pushover_hinge.hinge, backbones_pushover.element, ele_prop_table, node, analysis, joint  )
        end
    end
end

%% Dynamic Analysis
if sum(analysis.type_list == 1) > 0 % Dynamic Analysis was run as part of this proceedure 
    %% Haselton Mechanism Plots
    fn_curt_plot_2D( hinge, element, node, story, 'Collapse Mechanism', plot_dir )
    
    %% Plot Hinge accpetance
    % Elevation
    fn_plot_elevation( hinge, element, node, 'Frame 1', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'x', 0, 4000, 0, 0 )
%     fn_plot_elevation( hinge, element, node, 'Frame 1', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'x', 71, 1571, 0, 0 )
%     fn_plot_elevation( hinge, element, node, 'Frame 2', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'x', 71, 1571, 300, 300 )
%     fn_plot_elevation( hinge, element, node, 'Frame 3', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'x', 71, 1571, 600, 600 )
%     fn_plot_elevation( hinge, element, node, 'Frame 4', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'x', 71, 1571, 900, 900 )
%     fn_plot_elevation( hinge, element, node, 'West Upper Wall', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'z', 0, 0, 0, 900 )
%     fn_plot_elevation( hinge, element, node, 'Lower Wall 1', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'z', 71, 71, 0, 900 )
%     fn_plot_elevation( hinge, element, node, 'Lower Wall 2', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'z', 671, 671, 0, 900 )
%     fn_plot_elevation( hinge, element, node, 'Lower Wall 3', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'z', 971, 971, 0, 900 )
%     fn_plot_elevation( hinge, element, node, 'Lower Wall 4', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'z', 1271, 1271, 0, 900 )
%     fn_plot_elevation( hinge, element, node, 'East Upper Wall', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'z', 1500, 2000, 0, 900 )
%     fn_plot_elevation( joint, element, node, 'Frame 1 - Joint', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'x', 71, 1571, 0, 0 )
%     fn_plot_elevation( joint, element, node, 'Frame 2 - Joint', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'x', 71, 1571, 300, 300 )
%     fn_plot_elevation( joint, element, node, 'Frame 3 - Joint', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'x', 71, 1571, 600, 600 )
%     fn_plot_elevation( joint, element, node, 'Frame 4 - Joint', [plot_dir filesep 'Acceptance Plots' filesep 'elevations'], 'x', 71, 1571, 900, 900 )
% 
    % Plan View
%     fn_plot_plan_view( hinge, element, ele_prop_table, node, 1, 'Story 1 - Column Bases', [plot_dir filesep 'Acceptance Plots' filesep 'columns'])
%     fn_plot_plan_view( hinge, element, ele_prop_table, node, 2, 'Story 1 - Top of Columns', [plot_dir filesep 'Acceptance Plots' filesep 'columns'])
% 
%     % Plot Element Scatter
    fn_plot_element_scatter( element, 'column', story, hinge, plot_dir )
    fn_plot_element_scatter( element, 'beam', story, hinge, plot_dir )
%     fn_plot_element_scatter( element, 'wall', story, hinge, plot_dir )

    %% ASCE 41 Target Displacement
%     if analysis.run_eigen
%         [ target_disp_in.x ] = fn_spectra_and_target_displ( model, story, ground_motion, 'x' );
%         if strcmp(model.dimension,'3D')
%             if isfield(ground_motion,'z')
%                 [ target_disp_in.z ] = fn_spectra_and_target_displ( model, story, ground_motion, 'z' );
%             else
%               target_disp_in.z = NaN;
%             end
%         end  
%     else
        target_disp_in.x = NaN;
        target_disp_in.z = NaN;
%     end
    
    %% Load in Recordings to compare with EDPs and Time Histories
    if analysis.plot_recordings
        temp = load([pwd filesep 'ground_motions' filesep 'ICSB_recordings' filesep 'recorded_edp_profile.mat']);
        record_edp = temp.record_edp;
        pile_model = 0;
        if model.foundation == 2
            record_edp = temp.ff_edp;
            pile_model = 1;
        end
        
%         % Plot Spectra
%         fn_plot_spectra(node, 'Ground', ground_motion, plot_dir, read_dir_TH, pile_model)
%         fn_plot_spectra(node, 'Roof', ground_motion, plot_dir, read_dir_TH, pile_model)
    else
        record_edp = [];
    end
    
    % Plot EDP Profiles
    fn_plot_edp_profiles( plot_dir, ground_motion, story, target_disp_in, record_edp)
% 
%     % Plot specific TH comparisons
%     fn_plot_main_response_history( plot_dir, read_dir_TH, node, analysis, eq_analysis_timespace, eq, ground_motion.x.eq_dt, record_edp )
%     
%     % Plan view twist
%     fn_plot_plan_twist( node, read_dir_opensees, [plot_dir filesep 'EDP Profiles'], 'Roof Max Twist - Plan', record_edp )
        
    %% Plots Element results for Dynamic analysis
    if analysis.element_plots
        % Plot PM Diagrams for each element
%         fn_plot_PM_response( plot_dir, read_dir, element, analysis.hinge_stories_2_plot )

        % Plot Hinge Response
        if exist('backbones','var')
            load([read_dir filesep 'hinge_analysis.mat'])
            fn_plot_hinge_response( plot_dir, read_dir_opensees, hinge, backbones.element, ele_prop_table, node, analysis, joint )
        end
    end
    
    %% Compile outputs table and write csv file
    idx = 0;
    if analysis.plot_recordings
        for i = 1:length(dirs_ran)
            idx = idx + 1;
            outputs.type{idx,:} = ['Record-' dirs_ran{i}];
            outputs.period(idx,:) = NaN; % Maybe figures out a way to do this automaticly
            outputs.max_roof_drift(idx,:) = record_edp.max_disp.(dirs_ran{i})(end) / sum(story.story_ht);
            outputs.max_roof_drift_center(idx,:) = record_edp.max_disp_center.(dirs_ran{i})(end) / sum(story.story_ht);
            rel_disp = record_edp.max_disp.(dirs_ran{i})(2:end) - record_edp.max_disp.(dirs_ran{i})(1:(end-1));
            outputs.max_drift(idx,:) = max( rel_disp(1:height(story)) ./ story.story_ht);
            outputs.max_first_story_drift_center(idx,:) = record_edp.max_disp_center.(dirs_ran{i})(2) / story.story_ht(1);
            end_of_record = 5/ground_motion.(dirs_ran{i}).eq_dt;
            roof_disp_TH = record_edp.disp_TH_roof.(dirs_ran{i});
            outputs.residual_drift(idx,:) = abs(mean(roof_disp_TH((end-end_of_record):end))) / sum(story.story_ht);
            outputs.max_base_shear(idx,:) = NaN;
            outputs.max_roof_accel(idx,:) = record_edp.max_accel.(dirs_ran{i})(end);   
            outputs.max_roof_accel_center(idx,:) = record_edp.max_accel_center.(dirs_ran{i})(end); 
        end
    end
    
    for i = 1:length(dirs_ran)
        idx = idx + 1;
        outputs.type{idx,:} = ['Analysis-' dirs_ran{i}];
        if analysis.run_eigen
            outputs.period(idx,:) = model.(['T1_' dirs_ran{i}]);
        end
        outputs.max_roof_drift(idx,:) = story.(['max_disp_' dirs_ran{i}])(end) / sum(story.story_ht);
        outputs.max_roof_drift_center(idx,:) = story.(['max_disp_center_' dirs_ran{i}])(end) / sum(story.story_ht);
        outputs.max_drift(idx,:) = max(story.(['max_drift_' dirs_ran{i}]));
        outputs.max_first_story_drift_center(idx,:) = story.(['max_disp_center_' dirs_ran{i}])(1) / story.story_ht(1);
        outputs.residual_drift(idx,:) = story.(['residual_disp_' dirs_ran{i}])(end) / sum(story.story_ht);
        outputs.max_base_shear(idx,:) = story.(['max_reaction_' dirs_ran{i}])(1);
        outputs.max_roof_accel(idx,:) = story.(['max_accel_' dirs_ran{i}])(end);
        outputs.max_roof_accel_center(idx,:) = story.(['max_accel_center_' dirs_ran{i}])(end);
    end
    
    % Save as csv
    outputs_table = struct2table(outputs);
    writetable(outputs_table,[analysis.out_dir filesep 'analysis_summary.csv'])
end



end

