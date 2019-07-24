function [ ] = fn_plot_response_history( data_dir1, data_dir2, time_vec, eq, eq_dt, plot_dir, name, time_cutoff, data_record)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import plotting_tools.fn_format_and_save_plot

matlab_colors =  [0, 0.4470, 0.7410;
                  0.8500, 0.3250, 0.0980;
                  0.9290, 0.6940, 0.1250;
                  0.4940, 0.1840, 0.5560;
                  0.4660, 0.6740, 0.1880;
                  0.3010, 0.7450, 0.9330;
                  0.6350, 0.0780, 0.1840];

%% Plot response history at the roof
if ~isempty(data_dir2) % Both directions of response
    % NS Direction
    subplot(2,1,1)
    hold on
    max_data = max(abs(data_dir2));
    line([0,time_cutoff],[0,0],'color','k')
    if exist('data_record','var')
        max_data = max(abs([data_dir1,data_record.z']));
    %     plot([7,7],[-max_data,max_data]*2,'--','color',[0.5,0.5,0.5],'LineWidth',1,'DisplayName','Columns Hinge')
    %     plot([9,9],[-max_data,max_data]*2,'--','color',[0.25,0.25,0.25],'LineWidth',1,'DisplayName','Strength Loss')
    %     plot([11.3,11.3],[-max_data,max_data]*1.5,'--','color',[0,0,0],'LineWidth',1,'DisplayName','Column Failure')
        max_eq_idx = min([length(eq),length(data_record.z)]);
        plot((1:max_eq_idx)*eq_dt,data_record.z(1:max_eq_idx),'--k','DisplayName','Recorded')
    end
    plot(time_vec,data_dir2,'color',matlab_colors(1,:),'DisplayName','Analysis')
    title('Plan NS (Y)')
    ylabel(name)
    xlim([0,time_cutoff])
    if max_data > 0
        ylim([-max_data,max_data]*1.2)  
    end
    box off

    % EW Direction
    subplot(2,1,2)
    hold on
    max_data = max(abs(data_dir1));
    line([0,time_cutoff],[0,0],'color','k')
    if exist('data_record','var')
        max_data = max(abs([data_dir1,data_record.x']));
    %     plot([7,7],[-max_data,max_data]*2,'--','color',[0.5,0.5,0.5],'LineWidth',1,'DisplayName','Columns Hinge')
    %     plot([9,9],[-max_data,max_data]*2,'--','color',[0.25,0.25,0.25],'LineWidth',1,'DisplayName','Strength Loss')
    %     plot([11.3,11.3],[-max_data,max_data]*1.5,'--','color',[0,0,0],'LineWidth',1,'DisplayName','Column Failure')
        max_eq_idx = min([length(eq),length(data_record.x)]);
        plot((1:max_eq_idx)*eq_dt,data_record.x(1:max_eq_idx),'--k','DisplayName','Recorded')
    end
    plot(time_vec,data_dir1,'color',matlab_colors(2,:),'DisplayName','Analysis')
    title('Plan EW (X)')
    ylabel(name)
    xlim([0,time_cutoff])
    if max_data > 0
        ylim([-max_data,max_data]*1.2)  
    end
    box off
else
    % Single Response Direction
    hold on
    max_data = max(abs(data_dir1));
    line([0,time_cutoff],[0,0],'color','k')
    if exist('data_record','var')
        max_data = max(abs([data_dir1,data_record']));
    %     plot([7,7],[-max_data,max_data]*2,'--','color',[0.5,0.5,0.5],'LineWidth',1,'DisplayName','Columns Hinge')
    %     plot([9,9],[-max_data,max_data]*2,'--','color',[0.25,0.25,0.25],'LineWidth',1,'DisplayName','Strength Loss')
    %     plot([11.3,11.3],[-max_data,max_data]*1.5,'--','color',[0,0,0],'LineWidth',1,'DisplayName','Column Failure')
        max_eq_idx = min([length(eq),length(data_record)]);
        plot((1:max_eq_idx)*eq_dt,data_record(1:max_eq_idx),'--k','DisplayName','Recorded')
    end
    plot(time_vec,data_dir1,'color',matlab_colors(1,:),'DisplayName','Analysis')
    title('Plan EW (X)')
    ylabel(name)
    xlim([0,time_cutoff])
    if max_data > 0
        ylim([-max_data,max_data]*1.2)  
    end
    box off
end

xlabel('Time (s)')
fn_format_and_save_plot( plot_dir, name, 0)

end

