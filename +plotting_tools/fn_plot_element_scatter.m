function [ ] = fn_plot_element_scatter( element, ele_type, story, hinge, write_dir )
% Description: Fn to create scatter plots of hinge acceptance results

% Created By: Dustin Cook
% Date Created: 3/11/2019

% Inputs:

% Outputs:

% Assumptions:


%% Initial Setup
% Define Plot Dir
plot_dir = [write_dir filesep 'Scatter Plots'];

%% Begin Method
for s = 1:height(story)
    ele_story = element(element.story == story.id(s),:); 
    % Elements of one type
    ele = ele_story(strcmp(ele_story.type,ele_type),:);
    ele2use = ele(ele.length >= 89,:); % Omit the nubs on the end of the frame (improve the way I am doing this)
    hins = hinge(ismember(hinge.element_id,ele2use.id),:);

    % Primary Direction
    prim_hin = hins(strcmp(hins.direction,'primary'),:);
    if ~isempty(prim_hin)
        side_1 = prim_hin(prim_hin.ele_side == 1,:);
        side_2 = prim_hin(prim_hin.ele_side == 2,:);
        
        if strcmp(ele_type,'wall')
            plot_name = ['d ratio - Primary - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'd_ratio', 'Max(\Delta)/"d"', ele_type )
        else
            % a ratio
            plot_name = ['a ratio - Primary - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'a_ratio', 'Max(\theta)/"a"', ele_type )

            plot_name = ['a ratio - Primary - Side 2 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_2, plot_dir, plot_name, 'a_ratio', 'Max(\theta)/"a"', ele_type )
            
            % b ratio
            plot_name = ['b ratio - Primary - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'b_ratio', 'Max(\theta)/"b"', ele_type )

            plot_name = ['b ratio - Primary - Side 2 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_2, plot_dir, plot_name, 'b_ratio', 'Max(\theta)/"b"', ele_type )
            
            % V ratio
            plot_name = ['V ratio - Primary - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'V_ratio', 'Max(V)/Vn', ele_type )

            plot_name = ['V ratio - Primary - Side 2 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_2, plot_dir, plot_name, 'V_ratio', 'Max(V)/Vn', ele_type )
            
            % M ratio
            plot_name = ['M ratio - Primary - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'M_ratio_pos', 'Max(M)/Mn', ele_type )

            plot_name = ['M ratio - Primary - Side 2 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_2, plot_dir, plot_name, 'M_ratio_pos', 'Max(M)/Mn', ele_type )
        end
    end

    % OOP Direction
    oop_hin = hins(strcmp(hins.direction,'oop'),:);
    if ~isempty(oop_hin)
        side_1 = oop_hin(oop_hin.ele_side == 1,:);
        side_2 = oop_hin(oop_hin.ele_side == 2,:);
        
        if strcmp(ele_type,'wall')
            % a ratio
            plot_name = ['a ratio - OOP - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'a_ratio', 'Max(\theta)/"a"', ele_type )

            % b ratio
            plot_name = ['b ratio - OOP - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'b_ratio', 'Max(\theta)/"b"', ele_type )

            % V ratio
            plot_name = ['V ratio - OOP - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'V_ratio', 'Max(V)/Vn', ele_type )

        else
            % a ratio
            plot_name = ['a ratio - OOP - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'a_ratio', 'Max(\theta)/"a"', ele_type )

            plot_name = ['a ratio - OOP - Side 2 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_2, plot_dir, plot_name, 'a_ratio', 'Max(\theta)/"a"', ele_type )

            % b ratio
            plot_name = ['b ratio - OOP - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'b_ratio', 'Max(\theta)/"b"', ele_type )

            plot_name = ['b ratio - OOP - Side 2 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_2, plot_dir, plot_name, 'b_ratio', 'Max(\theta)/"b"', ele_type )

            % V ratio
            plot_name = ['V ratio - OOP - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'V_ratio', 'Max(V)/Vn', ele_type )

            plot_name = ['V ratio - OOP - Side 2 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_2, plot_dir, plot_name, 'V_ratio', 'Max(V)/Vn', ele_type )
        end
    end
    
    % SRSS Direction
    srss_hin = hins(strcmp(hins.direction,'primary'),:);
    srss_hin_oop = hins(strcmp(hins.direction,'oop'),:);
    if ~isempty(srss_hin) && ~isempty(srss_hin_oop)
        plot_name = ['SRSS - ' 'Story - ' num2str(s) ' - ' ele_type];
        srss_hin.a_ratio = sqrt(srss_hin.a_ratio.^2 + srss_hin_oop.a_ratio.^2);
        srss_hin.b_ratio = sqrt(srss_hin.b_ratio.^2 + srss_hin_oop.b_ratio.^2);
        srss_hin.V_ratio = sqrt(srss_hin.V_ratio.^2 + srss_hin_oop.V_ratio.^2);
        clear srss_names
        [srss_names{1:height(srss_hin)}] = deal('SRSS');
        srss_hin.direction = srss_names';
        side_1 = srss_hin(srss_hin.ele_side == 1,:);
        side_2 = srss_hin(srss_hin.ele_side == 2,:);
        
        if strcmp(ele_type,'wall')
            % V ratio
            plot_name = ['V ratio - OOP - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'V_ratio', 'Max(V)/Vn', ele_type )
        else
            % a ratio
            plot_name = ['a ratio - SRSS - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'a_ratio', 'Max(\theta)/"a"', ele_type )

            plot_name = ['a ratio - SRSS - Side 2 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_2, plot_dir, plot_name, 'a_ratio', 'Max(\theta)/"a"', ele_type )

            % b ratio
            plot_name = ['b ratio - SRSS - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'b_ratio', 'Max(\theta)/"b"', ele_type )

            plot_name = ['b ratio - SRSS - Side 2 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_2, plot_dir, plot_name, 'b_ratio', 'Max(\theta)/"b"', ele_type )

            % V ratio
            plot_name = ['V ratio - SRSS - Side 1 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_1, plot_dir, plot_name, 'V_ratio', 'Max(V)/Vn', ele_type )

            plot_name = ['V ratio - SRSS - Side 2 - Story - ' num2str(s) ' - ' ele_type];
            fn_plot_results_scatter( side_2, plot_dir, plot_name, 'V_ratio', 'Max(V)/Vn', ele_type )
        end
    end
end

function [ ] = fn_plot_results_scatter( data, plot_dir, plot_name, target, y_lab, ele_type )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Import Packages
import plotting_tools.fn_format_and_save_plot

% Matlab colors
matlab_colors =  [0, 0.4470, 0.7410;
                  0.8500, 0.3250, 0.0980;
                  0.9290, 0.6940, 0.1250;
                  0.4940, 0.1840, 0.5560;
                  0.4660, 0.6740, 0.1880;
                  0.3010, 0.7450, 0.9330;
                  0.6350, 0.0780, 0.1840];

if strcmp(data.direction(1),'SRSS')
    plt_color = [0 0 0];
elseif strcmp(data.direction(1),'primary')
    if strcmp(data.ele_direction(1),'x')
        plt_color = matlab_colors(2,:);
    elseif strcmp(data.ele_direction(1),'z')
        plt_color = matlab_colors(1,:);
    end
elseif strcmp(data.direction(1),'oop')
    if strcmp(data.ele_direction(1),'x')
        plt_color = matlab_colors(1,:);
    elseif strcmp(data.ele_direction(1),'z')
        plt_color = matlab_colors(2,:);
    end
end
ele_nums = 1:height(data);

% Plot Scatter
hold on
plot([0.5,height(data)+0.5],[1,1],'--k')
scatter(ele_nums(data.damage_recorded >= 2),data.(target)(data.damage_recorded >= 2), 50, plt_color, 's', 'filled' )
scatter(ele_nums(data.damage_recorded <= 1),data.(target)(data.damage_recorded <= 1), 50, plt_color, 's' )
ylabel(y_lab)
xlim([0.5,height(data)+0.5])
ylim([0, inf])
xlabel([ele_type ' #'])
box on
set(gcf,'position',[200, 500, 350, 300]);
fn_format_and_save_plot( plot_dir, plot_name, 0 )

% % B or E Ratio
% subplot(3,1,2)
% hold on
% plot([0.5,height(data)+height(side_2)+0.5],[1,1],'--k')
% if strcmp(ele_type,'wall') && strcmp(direction,'primary')
%     scatter(ele_nums(data.damage_recorded == 3),data.e_ratio(data.damage_recorded == 3), 'o', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
%     scatter(ele_nums(data.damage_recorded == 2),data.e_ratio(data.damage_recorded == 2), 'o', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
%     scatter(ele_nums(data.damage_recorded == 1),data.e_ratio(data.damage_recorded == 1), 'o', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
%     scatter(ele_nums(data.damage_recorded == 0),data.e_ratio(data.damage_recorded == 0), 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
%     scatter(side_2_id(side_2.damage_recorded == 3),side_2.e_ratio(side_2.damage_recorded == 3), 'd', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
%     scatter(side_2_id(side_2.damage_recorded == 2),side_2.e_ratio(side_2.damage_recorded == 2), 'd', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
%     scatter(side_2_id(side_2.damage_recorded == 1),side_2.e_ratio(side_2.damage_recorded == 1), 'd', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
%     scatter(side_2_id(side_2.damage_recorded == 0),side_2.e_ratio(side_2.damage_recorded == 0), 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
%     ylabel('Max(\Delta)/"e"')
% else
%     scatter(ele_nums(data.damage_recorded == 3),data.b_ratio(data.damage_recorded == 3), 'o', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
%     scatter(ele_nums(data.damage_recorded == 2),data.b_ratio(data.damage_recorded == 2), 'o', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
%     scatter(ele_nums(data.damage_recorded == 1),data.b_ratio(data.damage_recorded == 1), 'o', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
%     scatter(ele_nums(data.damage_recorded == 0),data.b_ratio(data.damage_recorded == 0), 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
%     scatter(side_2_id(side_2.damage_recorded == 3),side_2.b_ratio(side_2.damage_recorded == 3), 'd', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
%     scatter(side_2_id(side_2.damage_recorded == 2),side_2.b_ratio(side_2.damage_recorded == 2), 'd', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
%     scatter(side_2_id(side_2.damage_recorded == 1),side_2.b_ratio(side_2.damage_recorded == 1), 'd', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
%     scatter(side_2_id(side_2.damage_recorded == 0),side_2.b_ratio(side_2.damage_recorded == 0), 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
%     ylabel('Max(\theta)/"b"')
% end
% set(gca,'XTickLabel',[])
% xlim([0.5,height(data)+height(side_2)+0.5])
% ylim([0, inf])
% fn_format_and_save_plot( plot_dir, plot_name, 2, 0 )
% 
% % V Ratio
% subplot(3,1,3)
% hold on
% plot([0.5,height(data)+height(side_2)+0.5],[1,1],'--k')
% scatter(ele_nums(data.damage_recorded == 3),data.V_ratio(data.damage_recorded == 3), 'o', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
% scatter(ele_nums(data.damage_recorded == 2),data.V_ratio(data.damage_recorded == 2), 'o', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
% scatter(ele_nums(data.damage_recorded == 1),data.V_ratio(data.damage_recorded == 1), 'o', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
% scatter(ele_nums(data.damage_recorded == 0),data.V_ratio(data.damage_recorded == 0), 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
% scatter(side_2_id(side_2.damage_recorded == 3),side_2.V_ratio(side_2.damage_recorded == 3), 'd', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
% scatter(side_2_id(side_2.damage_recorded == 2),side_2.V_ratio(side_2.damage_recorded == 2), 'd', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
% scatter(side_2_id(side_2.damage_recorded == 1),side_2.V_ratio(side_2.damage_recorded == 1), 'd', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
% scatter(side_2_id(side_2.damage_recorded == 0),side_2.V_ratio(side_2.damage_recorded == 0), 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
% xlabel([upper(ele_type(1)) ele_type(2:end) 's'])
% ylabel('Max(V)/Vn')
% xlim([0.5,height(data)+height(side_2)+0.5])
% ylim([0, inf])
% fn_format_and_save_plot( plot_dir, plot_name, 2, 1 )
end

end

