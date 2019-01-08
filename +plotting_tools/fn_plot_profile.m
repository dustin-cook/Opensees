function [ ] = fn_plot_profile( profile_x, story_ids, plot_dir, plot_name, name, plot_max, target_disp, num_stories, additional_data )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import plotting_tools.fn_format_and_save_plot

%% Plot EDP profile
hold on
plot(profile_x,story_ids,'DisplayName','Analysis')
xlabel(name)
if strcmp(name,'IDR')
    ylabel('Story')
else
	ylabel('Floor')
end
if exist('additional_data','var')
   scatter(additional_data,[0:(length(additional_data)-1)],'r','filled','DisplayName','Recorded')
end

% Add target dispalcement
if ~isnan(target_disp)
    scatter(target_disp,num_stories,75,'m','*','DisplayName','Target Displacement')
end

xlim([0,plot_max])
fn_format_and_save_plot( plot_dir, plot_name, 1 )

end

