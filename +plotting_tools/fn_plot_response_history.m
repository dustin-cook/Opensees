function [ ] = fn_plot_response_history( data1, data1_time, data2, eq, eq_dt, plot_dir, name, time_cutoff )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import plotting_tools.fn_format_and_save_plot

%% Plot response history at the roof
figure
hold on
max_data = max(abs([data1,data2']));
plot([7,7],[-max_data,max_data]*2,'--','color',[0.5,0.5,0.5],'LineWidth',1,'DisplayName','Columns Hinge')
plot([9,9],[-max_data,max_data]*2,'--','color',[0.25,0.25,0.25],'LineWidth',1,'DisplayName','Strength Loss')
plot([11.3,11.3],[-max_data,max_data]*1.5,'--','color',[0,0,0],'LineWidth',1,'DisplayName','Column Failure')
plot(data1_time,data1,'b','DisplayName','Analysis')
plot((1:min([length(eq),length(data2)]))*eq_dt,data2(1:min([length(eq),length(data2)])),'r','DisplayName','Recorded')
xlabel('time (s)')
ylabel(name)
xlim([0,time_cutoff])
if max_data > 0
    ylim([-max_data,max_data]*1.2)  
end
fn_format_and_save_plot( plot_dir, name, 2)

end

