function [ ] = fn_plot_plan_view( hinge, element, node, ele_side, plot_label, plot_dir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import plotting_tools.*

%% Begin Method
% Filter Elements
ele_2_use = element(element.story == 1 & strcmp(element.type,'column'),:);
node_2_use = node(ismember(node.id,ele_2_use.node_1),:);
hinge_2_use = hinge(ismember(hinge.element_id,ele_2_use.id) & hinge.ele_side == ele_side,:);

%% Plot Recorded Damage
fn_plot_plan_scatter( node_2_use, hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - Recorded Damage'], 1, plot_label )

%% Confusion Matrix
analysis_match_fail = min(abs((hinge_2_use.accept(strcmp(hinge_2_use.direction,'primary')) - 1) - hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary'))),1);
fn_plot_plan_scatter( node_2_use, analysis_match_fail, plot_dir, [plot_label, ' - Confusion Matrix'], 3, plot_label, 'ASCE 41 Acceptance Criteria', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')) )

%% Acceptance Criteria
fn_plot_plan_scatter( node_2_use, hinge_2_use.accept(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - ASCE 41 Acceptance Criteria'], 2, plot_label, 'ASCE 41 Acceptance Criteria', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')) )

%% Plot Linear DCRs
if any(strcmp('V_ratio_mod',hinge_2_use.Properties.VariableNames))
    % Modified Forces
    fn_plot_plan_scatter( node_2_use, hinge_2_use.V_ratio_mod(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - V ratio mod'], 0, plot_label, 'Max(V)/Vn', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
    fn_plot_plan_scatter( node_2_use, hinge_2_use.M_ratio_mod(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - M ratio mod'], 0, plot_label, 'Max(M)/Mn', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
    fn_plot_plan_scatter( node_2_use, hinge_2_use.P_ratio_asce41_mod(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - P ratio asce 41 mod'], 0, plot_label, 'P_{max} / P_{n}', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
    fn_plot_plan_scatter( node_2_use, hinge_2_use.P_ratio_nominal_mod(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - P ratio nominal mod'], 0, plot_label, 'P_{max} / f''_{c,n}A_g', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
    fn_plot_plan_scatter( node_2_use, hinge_2_use.P_ratio_expected_mod(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - P ratio expected mod'], 0, plot_label, 'P_{max} / f''_{c,e}A_g', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
    
    % cp acceptance
    fn_plot_plan_scatter( node_2_use, hinge_2_use.V_ratio_accept_cp(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - DCR - V accept'], 0, plot_label, 'Max(V)/V_n', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
    fn_plot_plan_scatter( node_2_use, hinge_2_use.M_ratio_accept_cp(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - DCR - M accept'], 0, plot_label, 'Max(M)/M_n', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
    fn_plot_plan_scatter( node_2_use, hinge_2_use.P_ratio_accept_cp(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - DCR - P accept'], 0, plot_label, 'Max(P)/P_n', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
end

%% Plot A Ratio
if any(strcmp('a_ratio',hinge_2_use.Properties.VariableNames))
    fn_plot_plan_scatter( node_2_use, hinge_2_use.a_ratio(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - A ratio'], 0, plot_label, 'Max(\theta)/"a"', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
    fn_plot_plan_scatter( node_2_use, hinge_2_use.a_ratio(strcmp(hinge_2_use.direction,'oop')), plot_dir, [plot_label, ' - A ratio OOP'], 0, plot_label, 'Max(\theta)/"a"', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'oop')), 2 )
    srss_value = sqrt(hinge_2_use.a_ratio(strcmp(hinge_2_use.direction,'primary')).^2 + hinge_2_use.a_ratio(strcmp(hinge_2_use.direction,'oop')).^2); % Update to 1.5 for linear value
    fn_plot_plan_scatter( node_2_use, srss_value, plot_dir, [plot_label, ' - A ratio SRSS'], 0, plot_label, 'Max(\theta)/"a"', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
end

%% Plot B Ratio
if any(strcmp('b_ratio',hinge_2_use.Properties.VariableNames))
    fn_plot_plan_scatter( node_2_use, hinge_2_use.b_ratio(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - B ratio'], 0, plot_label, 'Max(\theta)/"b"', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
    fn_plot_plan_scatter( node_2_use, hinge_2_use.b_ratio(strcmp(hinge_2_use.direction,'oop')), plot_dir, [plot_label, ' - B ratio OOP'], 0, plot_label, 'Max(\theta)/"b"', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'oop')), 2 )
    srss_value = sqrt(hinge_2_use.b_ratio(strcmp(hinge_2_use.direction,'primary')).^2 + hinge_2_use.b_ratio(strcmp(hinge_2_use.direction,'oop')).^2);
    fn_plot_plan_scatter( node_2_use, srss_value, plot_dir, [plot_label, ' - B ratio SRSS'], 0, plot_label, 'Max(\theta)/"b"', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
end

%% Plot V Ratio
fn_plot_plan_scatter( node_2_use, hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - V ratio'], 0, plot_label, 'Max(V)/Vn', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 1 )
if ~isempty(hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'oop')))
    fn_plot_plan_scatter( node_2_use, hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'oop')), plot_dir, [plot_label, ' - V ratio OOP'], 0, plot_label, 'Max(V)/Vn', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'oop')), 1 )
    srss_value = sqrt(hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'primary')).^2 + hinge_2_use.V_ratio(strcmp(hinge_2_use.direction,'oop')).^2);
    fn_plot_plan_scatter( node_2_use, srss_value, plot_dir, [plot_label, ' - V ratio SRSS'], 0, plot_label, 'Max(V)/Vn', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 1 )
end

%% Plot P Ratio
fn_plot_plan_scatter( node_2_use, hinge_2_use.P_ratio_expected(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - P ratio expected'], 0, plot_label, 'P_{max} / f''_{c,e}A_g', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 1 )
fn_plot_plan_scatter( node_2_use, hinge_2_use.P_ratio_nominal(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - P ratio nominal'], 0, plot_label, 'P_{max} / f''_{c,n}A_g', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 1 )
fn_plot_plan_scatter( node_2_use, hinge_2_use.P_ratio_asce41(strcmp(hinge_2_use.direction,'primary')), plot_dir, [plot_label, ' - P ratio asce41'], 0, plot_label, 'P_{max} / P_n', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 1 )

%% Plot TAR
story_nodes = node(node.story == 1 & node.record_accel == 1 & node.mass > 0 & node.z ~= 450,:);
fn_plot_plan_scatter( story_nodes, story_nodes.TAR_x, plot_dir, [plot_label, ' - TAR x'], 0, plot_label, 'TAR - EW Shaking', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
fn_plot_plan_scatter( story_nodes, story_nodes.TAR_z, plot_dir, [plot_label, ' - TAR z'], 0, plot_label, 'TAR - NS Shaking', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )
fn_plot_plan_scatter( story_nodes, story_nodes.TAR_srss, plot_dir, [plot_label, ' - TAR srrs'], 0, plot_label, 'TAR - SRSS', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), 2 )

%% B Value
max_value = ceil(max(ele_2_use.(['b_hinge_' num2str(ele_side)]))*100)/100;
fn_plot_plan_scatter( node_2_use, ele_2_use.(['b_hinge_' num2str(ele_side)]), plot_dir, [plot_label, ' - B Value'], 0, plot_label, 'b Hinge Rotation', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), max_value )

%% CP Value
max_value = ceil(max(ele_2_use.(['cp_' num2str(ele_side)]))*100)/100;
fn_plot_plan_scatter( node_2_use, ele_2_use.(['cp_' num2str(ele_side)]), plot_dir, [plot_label, ' - CP Value'], 0, plot_label, 'CP Hinge Rotation', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), max_value )

%% LS Value
max_value = ceil(max(ele_2_use.(['ls_' num2str(ele_side)]))*100)/100;
fn_plot_plan_scatter( node_2_use, ele_2_use.(['ls_' num2str(ele_side)]), plot_dir, [plot_label, ' - LS Value'], 0, plot_label, 'LS Hinge Rotation', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), max_value )

%% IO Value
max_value = ceil(max(ele_2_use.(['io_' num2str(ele_side)]))*100)/100;
fn_plot_plan_scatter( node_2_use, ele_2_use.(['io_' num2str(ele_side)]), plot_dir, [plot_label, ' - IO Value'], 0, plot_label, 'IO Hinge Rotation', hinge_2_use.damage_recorded(strcmp(hinge_2_use.direction,'primary')), max_value )

end

