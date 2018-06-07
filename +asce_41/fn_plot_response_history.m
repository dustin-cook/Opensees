function [ ] = fn_plot_response_history( disp, eq, eq_dt, plot_dir, name )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

figure
hold on
plot((1:length(eq))*eq_dt,disp(end,:))
xlabel('time (s)')
ylabel(name)
fn_format_and_save_plot( plot_dir, name, 1)

end

