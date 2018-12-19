function [ force_vector_pos, force_vector_neg, disp_vector_pos, disp_vector_neg ] = fn_define_backbone( ele, ele_props, hinge_props )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


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
end

