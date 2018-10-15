function [ ] = fn_plot_response_history( data1, data2, eq, eq_dt, plot_dir, name )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import plotting_tools.fn_format_and_save_plot

%% Plot response history at the roof
figure
hold on
plot((1:length(data1))*eq_dt,data1,'DisplayName','Analysis')
plot((1:min([length(eq),length(data2)]))*eq_dt,data2(1:min([length(eq),length(data2)])),'DisplayName','Recorded')
xlabel('time (s)')
ylabel(name)
fn_format_and_save_plot( plot_dir, name, 1)

end

