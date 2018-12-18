function [ ] = fn_plot_backbone( ele, ele_props, hinge_props, output_dir, plot_name, plot_style, hinge_disp_to_plot, hinge_force_to_plot)
% Plot the backbone curve from an ASCE 41 analysis in terms of normalized 
% moment and rotation or drift and normalized shear force

% Created By: Dustin Cook
% Data Created: December 17, 2018

% Inputs
% plot_style: 1 = single backbone with no response, 2 = double backbone with response

% Assumptions:
% 1) Post yeild hardening slope ...
% 2) Yeild point is equal to nominal capacity
% 3) Ultimate point is equal to the plastic or probabl capacity (not sure which one we are reffering to it as)
% 4) Elastic rotational stiffness from ibbarra et al 2005
% 5) Post peak slope goes straight from point C to E (peak to end)
% 6) Elastic shear stiffness from delta = PL/GA

%% Begin Method
% Import Plotting Tools
import plotting_tools.*

% For Beams and Columns plot rotational hinge
if strcmp(ele.type,'beam') || strcmp(ele.type,'column')
    elastic_stiffness = 6*ele_props.e*ele_props.iz/ele.length;
    theta_yeild_pos = ele.Mn_aci_pos/elastic_stiffness;
    theta_yeild_neg = ele.Mn_aci_pos/elastic_stiffness;
    Q_y_pos = ele.Mn_aci_pos;
    Q_ult_pos = ele.Mp_pos;
    Q_y_neg = ele.Mn_aci_neg;
    Q_ult_neg = ele.Mp_neg;
    post_yeild_slope_pos = min([((Q_ult_pos-Q_y_pos)/Q_y_pos)/hinge_props.a_hinge,0.1*(1/theta_yeild_pos)]);
    post_yeild_slope_neg = min([((Q_ult_neg-Q_y_neg)/Q_y_neg)/hinge_props.a_hinge,0.1*(1/theta_yeild_neg)]);
    force_vector_pos = [1, post_yeild_slope_pos*hinge_props.a_hinge+1,hinge_props.c_hinge];
    force_vector_neg = [1, post_yeild_slope_neg*hinge_props.a_hinge+1,hinge_props.c_hinge];
    disp_vector_pos = theta_yeild_pos + [0, hinge_props.a_hinge, hinge_props.b_hinge];
    disp_vector_neg = theta_yeild_neg + [0, hinge_props.a_hinge, hinge_props.b_hinge];
    if plot_style == 1
        plot([0,disp_vector_pos],[0,force_vector_pos],'k','LineWidth',2,'DisplayName','ASCE 41 Backone') % Don't need to worry about negative bending because this plot is normalized by Qy
        xlabel('Total Rotation (rad)')
        ylabel('Q/Qy')
    elseif plot_style == 2
        hold on
        % -hinge_rotation_to_plot+(10/11)*theta_yeild
        plot([fliplr(-(disp_vector_neg-theta_yeild_neg)),0,(disp_vector_pos-theta_yeild_pos)],[fliplr(-force_vector_neg*Q_y_neg),0,force_vector_pos*Q_y_pos]/1000,'k','LineWidth',2,'DisplayName','ASCE 41 Backone')
        plot(hinge_disp_to_plot,-hinge_force_to_plot/1000,'b','LineWidth',1,'DisplayName','Analysis');
        xlabel('Hinge Rotation (rad)')
        ylabel('Moment (k-in)')
    end
    
% For Walls Plot shear springs
elseif strcmp(ele.type,'wall')
    elastic_stiffness = 0.4*ele_props.g*ele_props.a/ele.length; 
    f1 = hinge_props.f_hinge*ele.Vn_aci;
    u1 = f1/elastic_stiffness;
    f2 = ele.Vn_aci;
    u2 = (hinge_props.g_hinge/100)*ele.length;
    f3 = ele.Vn_aci*1.001;
    u3 = (hinge_props.d_hinge/100)*ele.length;
    f4 = ele.Vn_aci*ele.c_hinge;
    u4 = (hinge_props.e_hinge/100)*ele.length;
    force_vector_pos = [f1,f2,f3,f4]/f1;
    disp_vector_pos = [u1,u2,u3,u4]/ele.length;
    if plot_style == 1
        plot([0,disp_vector_pos],[0,force_vector_pos],'k','LineWidth',2,'DisplayName','ASCE 41 Backone')
    elseif plot_style == 2
        hold on
        plot([fliplr(-disp_vector_pos),0,disp_vector_pos],[fliplr(-force_vector_pos),0,force_vector_pos],'k','LineWidth',2,'DisplayName','ASCE 41 Backone')
        plot(hinge_disp_to_plot/ele.length,-hinge_force_to_plot/Q_y,'b','LineWidth',1,'DisplayName','Analysis');
    end
    xlabel('Drift')
end

% Format and save plot      
ylabel('Q/Qy')
if plot_style == 1
    xlim([0,max(disp_vector_pos)*1.25])
    ylim([0,max(force_vector_pos)*1.25])
elseif plot_style == 2
    xlim([-max(disp_vector_neg)*1.25,max(disp_vector_pos)*1.25])
    ylim([-max(force_vector_neg*Q_y_neg)*1.25/1000,max(force_vector_pos*Q_y_pos)*1.25/1000])
end
fn_format_and_save_plot( [output_dir filesep 'hinge_plots' filesep] , plot_name, 1 )
end

