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
    ele2use = ele(ele.length >= 100,:); % Omit the nubs on the end of the frame (improve the way I am doing this)
    hins = hinge(ismember(hinge.element_id,ele2use.id),:);

    % Primary Direction
    prim_hin = hins(strcmp(hins.direction,'primary'),:);
    if ~isempty(prim_hin)
        side_1 = prim_hin(prim_hin.ele_side == 1,:);
        side_2 = prim_hin(prim_hin.ele_side == 2,:);
        plot_name = ['Primary - ' 'Story - ' num2str(s) ' - ' ele_type];

        fn_plot_results_scatter( side_1, side_2, plot_dir, plot_name, ele_type, 'primary' )
    end

    % OOP Direction
    oop_hin = hins(strcmp(hins.direction,'oop'),:);
    if ~isempty(oop_hin)
        side_1 = oop_hin(oop_hin.ele_side == 1,:);
        side_2 = oop_hin(oop_hin.ele_side == 2,:);
        plot_name = ['OOP - ' 'Story - ' num2str(s) ' - ' ele_type];
        
        fn_plot_results_scatter( side_1, side_2, plot_dir, plot_name, ele_type, 'oop' )
    end
    
    % SRSS Direction
    srss_hin = hins(strcmp(hins.direction,'primary'),:);
    srss_hin_oop = hins(strcmp(hins.direction,'oop'),:);
    if ~isempty(srss_hin) && ~isempty(srss_hin_oop)
        plot_name = ['SRSS - ' 'Story - ' num2str(s) ' - ' ele_type];
        srss_hin.a_ratio = sqrt(srss_hin.a_ratio.^2 + srss_hin_oop.a_ratio.^2);
        srss_hin.b_ratio = sqrt(srss_hin.b_ratio.^2 + srss_hin_oop.b_ratio.^2);
        srss_hin.V_ratio = sqrt(srss_hin.V_ratio.^2 + srss_hin_oop.V_ratio.^2);
        srss_hin.d_ratio = sqrt(srss_hin.d_ratio.^2 + srss_hin_oop.d_ratio.^2);
        srss_hin.e_ratio = sqrt(srss_hin.e_ratio.^2 + srss_hin_oop.e_ratio.^2);
        side_1 = srss_hin(srss_hin.ele_side == 1,:);
        side_2 = srss_hin(srss_hin.ele_side == 2,:);
        
        fn_plot_results_scatter( side_1, side_2, plot_dir, plot_name, ele_type, 'oop' )
    end
end

function [ ] = fn_plot_results_scatter( side_1, side_2, plot_dir, plot_name, ele_type, direction )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Import Packages
import plotting_tools.fn_format_and_save_plot

side_1_id = 1:height(side_1);
side_2_id = height(side_1)+1:height(side_1)+height(side_2);
% A of D Ratio
subplot(3,1,1)
hold on
plot([0.5,height(side_1)+height(side_2)+0.5],[1,1],'--k')
if strcmp(ele_type,'wall') && strcmp(direction,'primary')
    scatter(side_1_id(side_1.damage_recorded == 3),side_1.d_ratio(side_1.damage_recorded == 3), 'o', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
    scatter(side_1_id(side_1.damage_recorded == 2),side_1.d_ratio(side_1.damage_recorded == 2), 'o', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
    scatter(side_1_id(side_1.damage_recorded == 1),side_1.d_ratio(side_1.damage_recorded == 1), 'o', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
    scatter(side_1_id(side_1.damage_recorded == 0),side_1.d_ratio(side_1.damage_recorded == 0), 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
    scatter(side_2_id(side_2.damage_recorded == 3),side_2.d_ratio(side_2.damage_recorded == 3), 'd', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
    scatter(side_2_id(side_2.damage_recorded == 2),side_2.d_ratio(side_2.damage_recorded == 2), 'd', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
    scatter(side_2_id(side_2.damage_recorded == 1),side_2.d_ratio(side_2.damage_recorded == 1), 'd', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
    scatter(side_2_id(side_2.damage_recorded == 0),side_2.d_ratio(side_2.damage_recorded == 0), 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
    ylabel('Max(\Delta)/"d"')
else
    scatter(side_1_id(side_1.damage_recorded == 3),side_1.a_ratio(side_1.damage_recorded == 3), 'o', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
    scatter(side_1_id(side_1.damage_recorded == 2),side_1.a_ratio(side_1.damage_recorded == 2), 'o', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
    scatter(side_1_id(side_1.damage_recorded == 1),side_1.a_ratio(side_1.damage_recorded == 1), 'o', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
    scatter(side_1_id(side_1.damage_recorded == 0),side_1.a_ratio(side_1.damage_recorded == 0), 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
    scatter(side_2_id(side_2.damage_recorded == 3),side_2.a_ratio(side_2.damage_recorded == 3), 'd', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
    scatter(side_2_id(side_2.damage_recorded == 2),side_2.a_ratio(side_2.damage_recorded == 2), 'd', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
    scatter(side_2_id(side_2.damage_recorded == 1),side_2.a_ratio(side_2.damage_recorded == 1), 'd', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
    scatter(side_2_id(side_2.damage_recorded == 0),side_2.a_ratio(side_2.damage_recorded == 0), 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
    ylabel('Max(\theta)/"a"')
end
set(gca,'XTickLabel',[])
xlim([0.5,height(side_1)+height(side_2)+0.5])
ylim([0, inf])
fn_format_and_save_plot( plot_dir, plot_name, 2, 0 )

% B or E Ratio
subplot(3,1,2)
hold on
plot([0.5,height(side_1)+height(side_2)+0.5],[1,1],'--k')
if strcmp(ele_type,'wall') && strcmp(direction,'primary')
    scatter(side_1_id(side_1.damage_recorded == 3),side_1.e_ratio(side_1.damage_recorded == 3), 'o', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
    scatter(side_1_id(side_1.damage_recorded == 2),side_1.e_ratio(side_1.damage_recorded == 2), 'o', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
    scatter(side_1_id(side_1.damage_recorded == 1),side_1.e_ratio(side_1.damage_recorded == 1), 'o', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
    scatter(side_1_id(side_1.damage_recorded == 0),side_1.e_ratio(side_1.damage_recorded == 0), 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
    scatter(side_2_id(side_2.damage_recorded == 3),side_2.e_ratio(side_2.damage_recorded == 3), 'd', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
    scatter(side_2_id(side_2.damage_recorded == 2),side_2.e_ratio(side_2.damage_recorded == 2), 'd', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
    scatter(side_2_id(side_2.damage_recorded == 1),side_2.e_ratio(side_2.damage_recorded == 1), 'd', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
    scatter(side_2_id(side_2.damage_recorded == 0),side_2.e_ratio(side_2.damage_recorded == 0), 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
    ylabel('Max(\Delta)/"e"')
else
    scatter(side_1_id(side_1.damage_recorded == 3),side_1.b_ratio(side_1.damage_recorded == 3), 'o', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
    scatter(side_1_id(side_1.damage_recorded == 2),side_1.b_ratio(side_1.damage_recorded == 2), 'o', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
    scatter(side_1_id(side_1.damage_recorded == 1),side_1.b_ratio(side_1.damage_recorded == 1), 'o', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
    scatter(side_1_id(side_1.damage_recorded == 0),side_1.b_ratio(side_1.damage_recorded == 0), 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
    scatter(side_2_id(side_2.damage_recorded == 3),side_2.b_ratio(side_2.damage_recorded == 3), 'd', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
    scatter(side_2_id(side_2.damage_recorded == 2),side_2.b_ratio(side_2.damage_recorded == 2), 'd', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
    scatter(side_2_id(side_2.damage_recorded == 1),side_2.b_ratio(side_2.damage_recorded == 1), 'd', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
    scatter(side_2_id(side_2.damage_recorded == 0),side_2.b_ratio(side_2.damage_recorded == 0), 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
    ylabel('Max(\theta)/"b"')
end
set(gca,'XTickLabel',[])
xlim([0.5,height(side_1)+height(side_2)+0.5])
ylim([0, inf])
fn_format_and_save_plot( plot_dir, plot_name, 2, 0 )

% V Ratio
subplot(3,1,3)
hold on
plot([0.5,height(side_1)+height(side_2)+0.5],[1,1],'--k')
scatter(side_1_id(side_1.damage_recorded == 3),side_1.V_ratio(side_1.damage_recorded == 3), 'o', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
scatter(side_1_id(side_1.damage_recorded == 2),side_1.V_ratio(side_1.damage_recorded == 2), 'o', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
scatter(side_1_id(side_1.damage_recorded == 1),side_1.V_ratio(side_1.damage_recorded == 1), 'o', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
scatter(side_1_id(side_1.damage_recorded == 0),side_1.V_ratio(side_1.damage_recorded == 0), 'o', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
scatter(side_2_id(side_2.damage_recorded == 3),side_2.V_ratio(side_2.damage_recorded == 3), 'd', 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r' )
scatter(side_2_id(side_2.damage_recorded == 2),side_2.V_ratio(side_2.damage_recorded == 2), 'd', 'filled', 'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0] )
scatter(side_2_id(side_2.damage_recorded == 1),side_2.V_ratio(side_2.damage_recorded == 1), 'd', 'filled', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c' )
scatter(side_2_id(side_2.damage_recorded == 0),side_2.V_ratio(side_2.damage_recorded == 0), 'd', 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b' )
xlabel([upper(ele_type(1)) ele_type(2:end) 's'])
ylabel('Max(V)/Vn')
xlim([0.5,height(side_1)+height(side_2)+0.5])
ylim([0, inf])
fn_format_and_save_plot( plot_dir, plot_name, 2, 1 )
end

end

