function [ ] = fn_plot_profile( profile_x, profile_z, story_ids, plot_dir, plot_name, name, target_disp, additional_data )
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

%% Plot EDP profile
if ~isempty(profile_z)
    % NS Plot
    if strcmp(name,'Interstory Drift')
        subplot(1,2,1)
    else
        subplot(1,3,1)
    end
    hold on
    max_data = ceil(120*max(profile_z))/100;
    if exist('additional_data','var')
       plot(additional_data.z(~isnan(additional_data.z)),story_ids(~isnan(additional_data.z)),'k--o','MarkerFaceColor','k','DisplayName','Recorded')
       max_data = ceil(1.2*max([profile_z;additional_data.z]));
    end
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
    else
        subplot(1,3,2)
    end
end
% EW Plot
hold on
max_data = ceil(120*max(profile_x))/100;
if exist('additional_data','var')
   plot(additional_data.x(~isnan(additional_data.x)),story_ids(~isnan(additional_data.x)),'k--o','MarkerFaceColor','k','DisplayName','Recorded')
   max_data = ceil(1.2*max([profile_x;additional_data.x]));
end
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

if ~strcmp(name,'Interstory Drift')
    hL = legend('location','east');
    legend boxoff
    set(hL,'Position', [0.52 0.48 0.52 0.1])
end
fn_format_and_save_plot( plot_dir, plot_name, 0 )

end

