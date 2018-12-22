function [ ] = fn_plot_backbone( ele, ele_props, hinge_props, output_dir, plot_name, plot_style, hinge_disp_to_plot, hinge_force_to_plot)
% Plot the backbone curve from an ASCE 41 analysis in terms of normalized 
% moment and rotation or drift and normalized shear force

% Created By: Dustin Cook
% Data Created: December 17, 2018

% Inputs
% plot_style: 1 = single backbone with no response, 2 = double backbone with response

% Assumptions:
% 1) hinge stiffness modification, n = 10

%% Begin Method
% Import Tools
import plotting_tools.*
import asce_41.*

% For Beams and Columns and (walls controlled by flexure), plot rotational hinge
if strcmp(ele.type,'beam') || strcmp(ele.type,'column') || (strcmp(ele.type,'wall') && strcmp(ele.critical_mode,'flexure'))
    n = 10;
    [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'full', ele.Mn_aci_pos, ele.Mn_aci_neg, ele.Mp_pos, ele.Mp_neg, ele.length, ele_props.e, ele_props.iz, hinge_props, n );
    if plot_style == 1
        plot([0,rot_vec_pos],[0,moment_vec_pos/moment_vec_pos(1)],'k','LineWidth',2) % Don't need to worry about negative bending because this plot is normalized by Qy
        xlabel('Total Rotation (rad)')
        ylabel('Q/Qy')
    elseif plot_style == 2
        hold on
        elastic_element_disp_pos = rot_vec_pos(1)*(n/(n+1)); % removes contribution from elastic beam/column. assumes n = 10;
        elastic_element_disp_neg = rot_vec_neg(1)*(n/(n+1)); % removes contribution from elastic beam/column. assumes n = 10;
        plot([fliplr(-(rot_vec_neg-elastic_element_disp_neg)),0,(rot_vec_pos-elastic_element_disp_pos)],[fliplr(-moment_vec_neg),0,moment_vec_pos]/1000,'k','LineWidth',2,'DisplayName','ASCE 41 Backone') % Converted to K-in
        plot(hinge_disp_to_plot,-hinge_force_to_plot/1000,'b','LineWidth',1,'DisplayName','Analysis'); % transform from lb-in to kip-in
        xlabel('Hinge Rotation (rad)')
        ylabel('Moment (k-in)')
        xlim([-max(rot_vec_neg)*1.25,max(rot_vec_pos)*1.25])
        ylim([-max(moment_vec_neg)*1.25/1000,max(moment_vec_pos)*1.25/1000])
    end
    
% For Walls Contolled by shearPlot shear springs
elseif strcmp(ele.type,'wall') && strcmp(ele.critical_mode,'shear')
    [ force_vec, disp_vec ] = fn_define_backbone_shear( ele.Vn_aci, ele.length, ele_props.g, ele_props.av, hinge_props );
    if plot_style == 1
        plot([0,disp_vec/ele.length],[0,force_vec/force_vec(2)],'k','LineWidth',2)
        xlabel('Drift')
        ylabel('Q/Qy')
    elseif plot_style == 2
        hold on
        plot([fliplr(-disp_vec),0,disp_vec],[fliplr(-force_vec),0,force_vec]/1000,'k','LineWidth',2,'DisplayName','ASCE 41 Backone') % transform from lbs to kips
        plot(hinge_disp_to_plot,-hinge_force_to_plot/1000,'b','LineWidth',1,'DisplayName','Analysis'); % transform from lbs to kips
        xlabel('Hinge Displacement (in)')
        ylabel('Shear Force (k)')
        xlim([-max(disp_vec)*1.25,max(disp_vec)*1.25])
        ylim([-max(force_vec)*1.25/1000,max(force_vec)*1.25/1000])
    end    
end

% Format and save plot      
if plot_style == 1
    fn_format_and_save_plot( [output_dir filesep 'hinge_plots' filesep] , plot_name, 2 )
elseif plot_style == 2
    fn_format_and_save_plot( [output_dir filesep 'hinge_plots' filesep] , plot_name, 1 )
end

end

