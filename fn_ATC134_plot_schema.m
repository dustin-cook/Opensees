function [ ] = fn_ATC134_plot_schema( plt, plt_dir )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Set Colors
colormap jet

% Define Plot Linetypes
hline = findobj(plt, 'type', 'line');
set(hline,'LineWidth',1.5)

% Define Scatter Options
if strcmp(plt.Type,'scatter')
    plt.MarkerFaceColor = 'flat';
end

% Format plot stlyes
grid on
grid minor
box on

% Set Font
set(gca,'FontSize',15)

% Set Figure Size
set(gcf,'Position',[100 100 1100 600]);
saveas(plt,[plt_dir filesep 'test.png'])
close

end

