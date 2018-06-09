function [ ] = fn_plot_profile( profile_x, story_ids, plot_dir, plot_name, name, additional_data )
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

fn_format_and_save_plot( plot_dir, plot_name, 1 )


%% Add target dispalcement at a later point
% c0 = 1.3;
% targ_disp = Sd.x*c1.x*c2.x*c0;
% plot([0,targ_disp],[0,num_stories],'--r','DisplayName','Target Displacement')

end

