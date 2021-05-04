function [ ] = fn_format_and_save_plot( plot_dir, plot_name, plot_ops, save_and_close )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Check to see if plot dir exists 
if ~exist(plot_dir,'dir')
  mkdir(plot_dir)  
end

fnt_sz = 10;
set(gca,'FontSize',fnt_sz)
% plot_dims = [100 100 400 300];
plot_dims = [100 100 650 300];

% Format plot
if plot_ops == 1 % Default settings
    legend('Location','eastoutside')
    grid on
    grid minor
    box on
    plot_dims = [100 100 650 300];
elseif  plot_ops == 2
    grid on
    grid minor
    box on
elseif  plot_ops == 3
    legend('Location','northeast')
    box on
elseif  plot_ops == 4
    box on
elseif plot_ops == 5
    legend('Location','southeast')
    grid on
    grid minor
    box on
elseif  plot_ops == 6
    legend('Location','southeast')
    box on 
end

hold off

set(gcf,'Position',plot_dims)

% Save and Close
if ~exist('save_and_close','var') || save_and_close == 1
    savefig([plot_dir, filesep, plot_name '.fig'])
    saveas(gcf, [plot_dir, filesep, plot_name '.png'])
    close
end

end

