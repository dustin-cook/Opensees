function [ ] = fn_plot_element_scatter( element, ele_type, story, hinge, write_dir )
% Description: Fn to create scatter plots of hinge acceptance results

% Created By: Dustin Cook
% Date Created: 3/11/2019

% Inputs:

% Outputs:

% Assumptions:


%% Initial Setup
% Import Packages
import plotting_tools.fn_format_and_save_plot

% Define Plot Dir
plot_dir = [write_dir filesep 'Scatter Plots'];

%% Begin Method
for s = 1:height(story)
    ele_story = element(element.story == story.id(s),:); 
    % Elements of one type
    ele = ele_story(strcmp(ele_story.type,ele_type),:);
    hins = hinge(ismember(hinge.element_id,ele.id),:);

    % Primary Direction
    prim_hin = hins(strcmp(hins.direction,'primary'),:);
    if ~isempty(prim_hin)
        side_1 = prim_hin(prim_hin.ele_side == 1,:);
        side_2 = prim_hin(prim_hin.ele_side == 2,:);
        plot_name = ['Primary - ' 'Story - ' num2str(s) ' - ' ele_type];

        % A Ratio
        subplot(3,1,1)
        hold on
        plot([0.5,height(side_1)+height(side_2)+0.5],[1,1],'--k')
        scatter(1:height(side_1),side_1.a_ratio, 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
        scatter(height(side_1)+1:height(side_1)+height(side_2),side_2.a_ratio, 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
        ylabel('Max(\theta)/"a"')
        set(gca,'XTickLabel',[])
        xlim([0.5,height(side_1)+height(side_2)+0.5])
        ylim([0, inf])
        fn_format_and_save_plot( plot_dir, plot_name, 2, 0 )
        
        % B Ratio
        subplot(3,1,2)
        hold on
        plot([0.5,height(side_1)+height(side_2)+0.5],[1,1],'--k')
        scatter(1:height(side_1),side_1.b_ratio, 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
        scatter(height(side_1)+1:height(side_1)+height(side_2),side_2.b_ratio, 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
        ylabel('Max(\theta)/"b"')
        set(gca,'XTickLabel',[])
        xlim([0.5,height(side_1)+height(side_2)+0.5])
        ylim([0, inf])
        fn_format_and_save_plot( plot_dir, plot_name, 2, 0 )
        
        % V Ratio
        subplot(3,1,3)
        hold on
        plot([0.5,height(side_1)+height(side_2)+0.5],[1,1],'--k')
        scatter(1:height(side_1),side_1.V_ratio, 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
        scatter(height(side_1)+1:height(side_1)+height(side_2),side_2.V_ratio, 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
        xlabel([upper(ele_type(1)) ele_type(2:end) 's'])
        ylabel('Max(V)/Vn')
        xlim([0.5,height(side_1)+height(side_2)+0.5])
        ylim([0, inf])
        fn_format_and_save_plot( plot_dir, plot_name, 2, 1)
    end

    % OOP Direction
    oop_hin = hins(strcmp(hins.direction,'oop'),:);
    if ~isempty(oop_hin)
        side_1 = oop_hin(oop_hin.ele_side == 1,:);
        side_2 = oop_hin(oop_hin.ele_side == 1,:);
        plot_name = ['OOP - ' 'Story - ' num2str(s) ' - ' ele_type];

        % A Ratio
        subplot(3,1,1)
        hold on
        plot([0.5,height(side_1)+height(side_2)+0.5],[1,1],'--k')
        scatter(1:height(side_1),side_1.a_ratio, 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
        scatter(height(side_1)+1:height(side_1)+height(side_2),side_2.a_ratio, 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
        ylabel('Max(\theta)/"a"')
        set(gca,'XTickLabel',[])
        xlim([0.5,height(side_1)+height(side_2)+0.5])
        ylim([0, inf])
        fn_format_and_save_plot( plot_dir, plot_name, 2, 0 )
        
        % B Ratio
        subplot(3,1,2)
        hold on
        plot([0.5,height(side_1)+height(side_2)+0.5],[1,1],'--k')
        scatter(1:height(side_1),side_1.b_ratio, 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
        scatter(height(side_1)+1:height(side_1)+height(side_2),side_2.b_ratio, 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
        ylabel('Max(\theta)/"b"')
        set(gca,'XTickLabel',[])
        xlim([0.5,height(side_1)+height(side_2)+0.5])
        ylim([0, inf])
        fn_format_and_save_plot( plot_dir, plot_name, 2, 0 )
        
        % V Ratio
        subplot(3,1,3)
        hold on
        plot([0.5,height(side_1)+height(side_2)+0.5],[1,1],'--k')
        scatter(1:height(side_1),side_1.V_ratio, 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
        scatter(height(side_1)+1:height(side_1)+height(side_2),side_2.V_ratio, 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
        xlabel([upper(ele_type(1)) ele_type(2:end) 's'])
        ylabel('Max(V)/Vn')
        xlim([0.5,height(side_1)+height(side_2)+0.5])
        ylim([0, inf])
        fn_format_and_save_plot( plot_dir, plot_name, 2, 1 )
    end
end



end

