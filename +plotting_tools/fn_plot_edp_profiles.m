function [ ] = fn_plot_edp_profiles( plot_dir, ground_motion, story, target_disp_in, record_edp )
% Description: Fn to plot edp profiles of the analysis

% Created By: Dustin Cook
% Date Created: 1/7/2019

% Inputs:

% Outputs:

% Assumptions:

%% Initial Setup
% Plot Directory
edp_plot_dir = [plot_dir filesep 'EDP Profiles'];

%% Caclulate Profiles
% X accels
if any(strcmp('max_accel_x',story.Properties.VariableNames))
    edp_accel_x = [ground_motion.x.pga; story.max_accel_x];
%     edp_accel_x_center = [ground_motion.x.pga; story.max_accel_center_x];
end

% X disps
edp_disp_x = [0; story.max_disp_x];
% edp_disp_x_center = [0; story.max_disp_center_x];
if any(strcmp('max_disp_x_mod',story.Properties.VariableNames))
    edp_mods.max_disp.x = [0; story.max_disp_x_mod];
%     edp_mods.max_disp_center.x = [0; story.max_disp_center_x_mod];
else
    edp_mods.max_disp = [];
%     edp_mods.max_disp_center = [];
end

% Z direction
if isfield(ground_motion,'z')
    % Z accels
    if any(strcmp('max_accel_z',story.Properties.VariableNames))
        edp_accel_z = [ground_motion.z.pga; story.max_accel_z];
%         edp_accel_z_center = [ground_motion.z.pga; story.max_accel_center_z];
    end
    
    % Z disps
    edp_disp_z = [0; story.max_disp_z];
    edp_disp_z_center = [0; story.max_disp_center_z];
    if ismember('max_twist_z',story.Properties.VariableNames)
        edp_disp_z_twist = [0; story.max_twist_z];
    else
        edp_disp_z_twist = [];
    end
    if any(strcmp('max_disp_z_mod',story.Properties.VariableNames))
        edp_mods.max_disp.z = [0; story.max_disp_z_mod];
        edp_mods.max_disp_center.z = [0; story.max_disp_center_z_mod];
    end
else
    edp_accel_z = [];
%     edp_accel_z_center = [];
    edp_disp_z = [];
%     edp_disp_z_center = [];
    edp_disp_z_twist = [];
end

if isempty(record_edp)
    record_edp.max_accel = [];
%     record_edp.max_accel_center = [];
    record_edp.max_disp = [];
%     record_edp.max_disp_center = [];
    record_edp.max_twist = [];
end

%% Plot Profiles
% Acceleration
if any(strcmp('max_accel_x',story.Properties.VariableNames))
    fn_plot_profile( edp_accel_x, edp_accel_z, [min(story.id);(story.id+1)], edp_plot_dir, 'Acceleration Profile', 'Acceleration (g)', NaN, record_edp.max_accel, [])
%     fn_plot_profile( edp_accel_x_center, edp_accel_z_center,  [min(story.id);(story.id+1)], edp_plot_dir, 'Acceleration Profile Center', 'Acceleration (g)', NaN, record_edp.max_accel_center, [])
end
% Displacement
fn_plot_profile( edp_disp_x, edp_disp_z, [min(story.id);(story.id+1)], edp_plot_dir, 'Displacement Profile', 'Displacement (in)', target_disp_in, record_edp.max_disp, edp_mods.max_disp)
% if ~any(isnan(edp_disp_x_center))
%     fn_plot_profile( edp_disp_x_center, edp_disp_z_center, [min(story.id);(story.id+1)], edp_plot_dir, 'Displacement Profile Center', 'Displacement (in)', target_disp_in, record_edp.max_disp_center, edp_mods.max_disp_center)
% end
    
% Drift
if isfield(ground_motion,'z')
    fn_plot_profile( story.max_drift_x, story.max_drift_z, story.id, edp_plot_dir, 'Drift Profile', 'Interstory Drift', NaN, [], [] )
else
    fn_plot_profile( story.max_drift_x, [], story.id, edp_plot_dir, 'Drift Profile', 'Interstory Drift', NaN, [], [] )
end

% Twist
if ~isempty(edp_disp_z_twist)
    fn_plot_profile( [], edp_disp_z_twist, [min(story.id);(story.id+1)], edp_plot_dir, 'Twist Profile', 'Twist (in)', NaN, record_edp.max_twist, [])
end

end



function [ ] = fn_plot_profile( profile_x, profile_z, story_ids, plot_dir, plot_name, name, target_disp, recording, mod_profile )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import Packages
import plotting_tools.*

% Define Colors
matlab_colors =  [0, 0.4470, 0.7410;
                  0.8500, 0.3250, 0.0980;
                  0.9290, 0.6940, 0.1250;
                  0.4940, 0.1840, 0.5560;
                  0.4660, 0.6740, 0.1880;
                  0.3010, 0.7450, 0.9330;
                  0.6350, 0.0780, 0.1840];

%% NS Plot
if ~isempty(profile_z)
    if strcmp(name,'Interstory Drift') || strcmp(name,'Twist (in)')
        subplot(1,2,1)
    else
        subplot(1,3,1)
    end
    hold on
    max_data = ceil(120*max(profile_z))/100;
    
    % Plot Extra data
    if ~isempty(recording)
       plot(recording.z(~isnan(recording.z)),story_ids(~isnan(recording.z)),'k--o','MarkerFaceColor','k','DisplayName','Recorded')
       max_data = ceil(1.2*max([profile_z;recording.z]));
    end
    if ~isempty(mod_profile)
       plot(mod_profile.z,story_ids,'--','color',matlab_colors(1,:),'linewidth',1.5,'DisplayName','Linear Modifications')
    end
    
    % Plot main profile
    plot(profile_z,story_ids,'color',matlab_colors(1,:),'linewidth',1.5,'DisplayName','Analysis')
    % Add target dispalcement
    if isstruct(target_disp)
        scatter(target_disp.z,max(story_ids),75,matlab_colors(1,:),'*','DisplayName','Target Displacement')
    end
    xlabel(['NS ' name])
    if strcmp(name,'Interstory Drift')
        ylabel('Story')
    else
        ylabel('Floor Level')
    end
    xlim([0,max_data])
    ylim([min(story_ids),max(story_ids)])
    yticks(story_ids)
    % legend('location','northwest')
    % legend boxoff
    box on

    % Subplot space for EW plot
    if strcmp(name,'Interstory Drift')
        subplot(1,2,2)
    elseif ~strcmp(name,'Twist (in)')
        subplot(1,3,2)
    end
end

%% EW Plot
if ~isempty(profile_x)
    hold on
    max_data = ceil(120*max(profile_x))/100;

    % Additional data
    if ~isempty(recording)
       plot(recording.x(~isnan(recording.x)),story_ids(~isnan(recording.x)),'k--o','MarkerFaceColor','k','DisplayName','Recorded')
       max_data = ceil(1.2*max([profile_x;recording.x]));
    end
    if ~isempty(mod_profile)
       plot(mod_profile.x,story_ids,'--','color',matlab_colors(2,:),'linewidth',1.5,'DisplayName','Linear Modifications')
    end

    % Main data
    plot(profile_x,story_ids,'color',matlab_colors(2,:),'linewidth',1.5,'DisplayName','Analysis')

    % Add target dispalcement
    if isstruct(target_disp)
        scatter(target_disp.x,max(story_ids),75,matlab_colors(2,:),'*','DisplayName','Target Displacement')
    end
    xlabel(['EW ' name])
    xlim([0,max_data])
    ylim([min(story_ids),max(story_ids)])
    yticks(story_ids)
    box on
end

if ~strcmp(name,'Interstory Drift')
    hL = legend('location','east');
    legend boxoff
    set(hL,'Position', [0.52 0.48 0.52 0.1])
end
fn_format_and_save_plot( plot_dir, plot_name, 0 )

end

