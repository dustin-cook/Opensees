function [ ] = fn_plot_backbone( ele, ele_side, ele_props, read_dir, output_dir, plot_name, plot_style, hinge_disp_to_plot, hinge_force_to_plot, crit_mode, hin_dir, line_color )
% Plot the backbone curve from an ASCE 41 analysis in terms of normalized 
% moment and rotation or drift and normalized shear force

% Created By: Dustin Cook
% Data Created: December 17, 2018

% Inputs
% plot_style: 1 = single backbone with no response, 2 = double backbone with response

% Assumptions:
% 1) hinge stiffness modification, n = 10

%% Begin Method
matlab_colors =  [0, 0.4470, 0.7410;
                  0.8500, 0.3250, 0.0980;
                  0.9290, 0.6940, 0.1250;
                  0.4940, 0.1840, 0.5560;
                  0.4660, 0.6740, 0.1880;
                  0.3010, 0.7450, 0.9330;
                  0.6350, 0.0780, 0.1840];
              
if strcmp(ele.direction,'x')
    plt_color = matlab_colors(2,:);
elseif strcmp(ele.direction,'z')
    plt_color = matlab_colors(1,:);
end
              
% Import Tools
import plotting_tools.*
import asce_41.*

% For Beams and Columns and (walls controlled by flexure), plot rotational hinge
if strcmp(ele.type,'beam') || strcmp(ele.type,'column') || (strcmp(ele.type,'wall') && strcmp(crit_mode,'flexure'))
    n = 10;
    if strcmp(hin_dir,'oop')
        [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'full', ele.(['Mn_oop_' ele_side]), ele.(['Mn_oop_' ele_side]), ele.(['Mp_oop_' ele_side]), ele.(['Mp_oop_' ele_side]), ele.length, ele_props.e, ele_props.iy, ele.(['a_hinge_oop_' ele_side]), ele.(['b_hinge_oop_' ele_side]), ele.(['c_hinge_oop_' ele_side]), n, 0.1, ele.(['critical_mode_oop_' ele_side]) );
        if strcmp(ele.direction,'x')
            ele_dir = 'z';
        elseif strcmp(ele.direction,'z')
            ele_dir = 'x';
        end
        mom_of_I = ele_props.iy;
    else
        [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'full', ele.(['Mn_pos_' ele_side]), ele.(['Mn_neg_' ele_side]), ele.(['Mp_pos_' ele_side]), ele.(['Mp_neg_' ele_side]), ele.length, ele_props.e, ele_props.iz, ele.(['a_hinge_' ele_side]), ele.(['b_hinge_' ele_side]), ele.(['c_hinge_' ele_side]), n, 0.1, ele.(['critical_mode_' ele_side]) );
        ele_dir = ele.direction{1};
        mom_of_I = ele_props.iz;
    end
    if plot_style == 1
        plot([0,rot_vec_pos],[0,moment_vec_pos/moment_vec_pos(1)],'Color',line_color,'LineWidth',2) % Don't need to worry about negative bending because this plot is normalized by Qy
        xlabel('Total Rotation (rad)')
        ylabel('Q/Qy')
    elseif plot_style == 2
        hold on
        yax = plot([0,0],[-1e6,1e6],'color',[0.5, 0.5, 0.5]);
        set(get(get(yax,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        xax = plot([-1e6,1e6],[0,0],'color',[0.5, 0.5, 0.5]);
        set(get(get(xax,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        plot([fliplr(-(rot_vec_neg)),0,rot_vec_pos],[fliplr(-moment_vec_neg),0,moment_vec_pos]/1000,'Color','k','LineWidth',2,'DisplayName','ASCE 41 Backone') % Don't need to worry about negative bending because this plot is normalized by Qy
        elastic_rotation = (hinge_force_to_plot*ele.length / (6*ele_props.e*mom_of_I))*(10/11);
        plot(hinge_disp_to_plot + elastic_rotation,hinge_force_to_plot/1000,'color', plt_color, 'LineWidth',1.25,'DisplayName','Analysis'); % transform from lb-in to kip-in
        xlabel('Total Rotation (rad)')
        ylabel('Moment (k-in)')
%         xlim([-max(rot_vec_neg)*1.25,max(rot_vec_pos)*1.25])
%         ylim([-max(moment_vec_neg)*1.25/1000,max(moment_vec_pos)*1.25/1000])
        xlim([-0.05,0.05])
        ylim([-15000,15000])
    end
    
% For Walls Contolled by shear, Plot shear springs
elseif strcmp(ele.type,'wall') && strcmp(crit_mode,'shear')
    [ force_vec, disp_vec ] = fn_define_backbone_shear( ele.(['Vn_' ele_side]), ele.length, ele_props.g, ele_props.av, ele.(['c_hinge_' ele_side]), ele.(['d_hinge_' ele_side]), ele.(['e_hinge_' ele_side]), ele.(['f_hinge_' ele_side]), ele.(['g_hinge_' ele_side])  );
    if plot_style == 1
        plot([0,disp_vec/ele.length],[0,force_vec/force_vec(2)],'Color',line_color,'LineWidth',2)
        xlabel('Drift')
        ylabel('Q/Qy')
    elseif plot_style == 2
        hold on
        yax = plot([0,0],[1e6,1e6],'color',[0.5, 0.5, 0.5]);
        set(get(get(yax,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        xax = plot([1e6,1e6],[0,0],'color',[0.5, 0.5, 0.5]);
        set(get(get(xax,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        plot([fliplr(-disp_vec),0,disp_vec],[fliplr(-force_vec),0,force_vec]/1000,'k','LineWidth',1.5,'DisplayName','ASCE 41 Backone') % transform from lbs to kips
        plot(hinge_disp_to_plot,hinge_force_to_plot/1000,'color', plt_color,'LineWidth',1.25,'DisplayName','Analysis'); % transform from lbs to kips
        xlabel('Hinge Displacement (in)')
        ylabel('Shear Force (k)')
        xlim([-max(disp_vec)*1.25,max(disp_vec)*1.25])
        ylim([-max(force_vec)*1.25/1000,max(force_vec)*1.25/1000])
    end    
end

% Format and save plot
box on
% legend('location','northeast')
% legend boxoff
set(gcf,'Position',[300 400 300 300]);
if plot_style == 1
    fn_format_and_save_plot( output_dir , plot_name, 2 )
elseif plot_style == 2
    fn_format_and_save_plot( output_dir , plot_name, 0 )
end

end

