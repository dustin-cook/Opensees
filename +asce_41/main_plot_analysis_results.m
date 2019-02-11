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
fn_make_directory( plot_dir )

% Load Analysis Data
load([read_dir filesep 'element_analysis.mat'])
load([read_dir filesep 'element_TH.mat']);
load([read_dir filesep 'story_analysis.mat'])
load([read_dir filesep 'hinge_analysis.mat'])
load([read_dir filesep 'model_analysis.mat'])

load([read_dir_opensees filesep 'node_analysis.mat']) % Could probably change these to read from model inputs instead of the opensses analysis
if sum(analysis.type_list == 1) > 0 % Dynamic Analysis was run as part of this proceedure
    load([read_dir_opensees filesep 'gm_data.mat'])
end

if analysis.asce_41_post_process
    load([read_dir filesep 'element_PM.mat']) % Maybe move this to where its needed
end

if analysis.plot_recordings
    load([pwd filesep 'ground_motions' filesep 'ICSB_recordings' filesep 'recorded_edp_profile.mat']);
end

if sum(strcmp(analysis.case_list,'backbones')) > 0
    backbones = load([analysis.out_dir filesep 'backbones' filesep 'element_analysis.mat']);
end

%% Cyclic Analysis
if sum(analysis.type_list == 3) > 0 % Cyclic was run as part of this proceedure
    if analysis.element_plots
        % Load Pushover Results
        cyclic_read_dir = [analysis.out_dir filesep 'cyclic'];
        cyclic_analysis = load([cyclic_read_dir filesep 'analysis_options.mat']);
        if cyclic_analysis.analysis.nonlinear ~= 0 % Is the pushover nonlinear?
            % Plot Hinge Response
            cyclic_hinge = load([cyclic_read_dir filesep 'hinge_analysis.mat']);
            fn_plot_hinge_response( cyclic_read_dir, cyclic_hinge.hinge, backbones.element, ele_prop_table, node, analysis.hinge_stories_2_plot )
        end
    end
end

%% Pushover Analysis
if sum(analysis.type_list == 2) > 0 % Pushover was run as part of this proceedure
    pushover_read_dir = [analysis.out_dir filesep 'pushover'];
    pushover_analysis = load([pushover_read_dir filesep 'analysis_options.mat']);
    % Plot Building Pushovers
    fn_plot_pushover( pushover_read_dir, 'x', story.story_dead_load, pushover_analysis.analysis.base_shear_x )
    if strcmp(model.dimension,'3D')
        fn_plot_pushover( pushover_read_dir, 'z', story.story_dead_load, pushover_analysis.analysis.base_shear_z )
    end
    
    if analysis.element_plots
        % Load Pushover Results
        pushover_TH = load([pushover_read_dir filesep 'element_TH.mat']);

        % Plot PM Diagrams for each element
        fn_plot_PM_response( pushover_read_dir, element, pushover_TH.element_TH, element_PM, analysis.hinge_stories_2_plot )

        % Plot Hinge Response
        if pushover_analysis.analysis.nonlinear ~= 0 % Is the pushover nonlinear?
            pushover_hinge = load([pushover_read_dir filesep 'hinge_analysis.mat']);
            fn_plot_hinge_response( pushover_read_dir, pushover_hinge.hinge, backbones.element, ele_prop_table, node, analysis.hinge_stories_2_plot )
        end
    end
end

%% Dynamic Analysis
if sum(analysis.type_list == 1) > 0 % Dynamic Analysis was run as part of this proceedure
    if analysis.asce_41_post_process    
        %% ASCE 41 Acceptance Plots
        if strcmp(analysis.proceedure,'LDP') % Linear Procedures
            %% Plot DCR
            plot_variable = {'all', 'M', 'V', 'P'};
            frame_lines = {'ext', 'int'};
            if strcmp(model.dimension,'3D')% for 3D linear analysis
                for i = 1:length(plot_variable)
                    for j = 1:length(frame_lines)
                        fn_plot_building( element.(['DCR_max_' plot_variable{i}]), element, node, ['DCR_view_' plot_variable{i} '_' frame_lines{j}], plot_dir, '3D', 'linear', frame_lines{j} )
                        fn_plot_building( element.(['DCR_raw_max_' plot_variable{i}]), element, node, ['DCR_view_' plot_variable{i} '_' frame_lines{j} '_raw'], plot_dir, '3D', 'raw', frame_lines{j} )
                    end
                end

            else
                for i = 1:length(plot_variable)
                    fn_plot_building_2D( element.(['DCR_max_' plot_variable{i}]), element, node, ['DCR_view_' plot_variable{i}], plot_dir,  'linear' )
                    fn_plot_building_2D( element.(['DCR_raw_max_' plot_variable{i}]), element, node, ['DCR_view_' plot_variable{i} '_raw'], plot_dir,  'raw' )
                end
            end
        elseif strcmp(analysis.proceedure,'NDP') % Nonlinear Procedures
            %% Plot Hinge accpetance
            fn_plot_building_nl_3d( hinge, element, node, 'Acceptance Plot - Interior Frame', plot_dir, 'int_frame' )
            fn_plot_building_nl_3d( hinge, element, node, 'Acceptance Plot - Exterior Frame', plot_dir, 'ext_frame' )
            fn_plot_building_nl_3d( hinge, element, node, 'Acceptance Plot - East Wall Frame', plot_dir, 'east_wall' )
%             fn_plot_building_nl( hinge, element, node, 'Acceptance Plot', plot_dir)
        end
        
        %% ASCE 41 Target Displacement
        [ target_disp_in.x ] = fn_spectra_and_target_displ( model, story, ground_motion, 'x' );
        if strcmp(model.dimension,'3D')
            if isfield(ground_motion,'z')
                [ target_disp_in.z ] = fn_spectra_and_target_displ( model, story, ground_motion, 'z' );
            else
              target_disp_in.z = NaN;
            end
        end        
    else
        target_disp_in.x = NaN;
        if strcmp(model.dimension,'3D')
            target_disp_in.z = NaN;
        end 
    end

    %% Load in Recordings to compare with EDPs and Time Histories
    if analysis.plot_recordings
        % Plot EDP Profiles
        fn_plot_edp_profiles( plot_dir, ground_motion.x.pga, model, story, target_disp_in.x, 'x', record_edp)
        if strcmp(model.dimension,'3D') && isfield(ground_motion,'z')
            fn_plot_edp_profiles( plot_dir, ground_motion.z.pga, model, story, target_disp_in.z, 'z', record_edp )
        end

        % Plot specific TH comparisons
        fn_plot_main_response_history( plot_dir, node, analysis, eq_analysis_timespace, eq.x, ground_motion.x.eq_dt, 'x', record_edp )
        if strcmp(model.dimension,'3D') && isfield(ground_motion,'z')
            fn_plot_main_response_history( plot_dir, node, analysis, eq_analysis_timespace, eq.z, ground_motion.z.eq_dt, 'z', record_edp )
        end
    else
        % Plot EDP Profiles
        fn_plot_edp_profiles( plot_dir, ground_motion.x.pga, model, story, target_disp_in.x, 'x' )
        if strcmp(model.dimension,'3D') && isfield(ground_motion,'z')
            fn_plot_edp_profiles( plot_dir, ground_motion.z.pga, model, story, target_disp_in.z, 'z' )
        end

        % Plot specific TH comparisons
        fn_plot_main_response_history( plot_dir, node, analysis, eq_analysis_timespace, eq.x, ground_motion.x.eq_dt, 'x' )
        if strcmp(model.dimension,'3D') && isfield(ground_motion,'z')
            fn_plot_main_response_history( plot_dir, node, analysis, eq_analysis_timespace, eq.z, ground_motion.z.eq_dt, 'z' )
        end

    end

    %% Plots Element results for both Dynamic analysis
    if analysis.element_plots
        % Plot PM Diagrams for each element
        fn_plot_PM_response( plot_dir, element, element_TH, element_PM, analysis.hinge_stories_2_plot )

        % Plot Hinge Response
        load([read_dir filesep 'hinge_analysis.mat'])
        fn_plot_hinge_response( plot_dir, hinge, backbones.element, ele_prop_table, node, analysis.hinge_stories_2_plot )
    end
end

end

