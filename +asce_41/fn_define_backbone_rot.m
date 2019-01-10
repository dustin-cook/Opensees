function [ moment_vec_pos, moment_vec_neg, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( type, Mn_pos, Mn_neg, Mp_pos, Mp_neg, length, e, iz, hinge_props, n, strain_harden_ratio )
% Takes points from ASCE 41-17 chapter 10 tables and Define backbone curve
% as a vector of displacements and forces for both the positive and
% negative directions. Units are based on the basic units that are input.

% Created By: Dustin Cook
% Data Created: December 20, 2018

% Inputs
% 

% Outputs
% moment_vec: Vector of moment points along the backbone curve for the
%               positive or negative direction. Does not include the
%               origin.
% rot_vec: Vector of rotation points along the backbone curve for the
%               positive or negative direction. Does not include the
%               origin.

% Assumptions:
% 1) Post yeild hardening slope is the min of 10% elastic or the full increase to probable moment after a_hinge length
% 2) Yeild point is equal to nominal capacity
% 3) Ultimate point is equal to the probable capacity
% 4) Elastic rotational stiffness from ibbarra et al 2005 based on stiffness matrix
% 5) Post peak slope goes straight from point C to E (peak to end)
% 6) Backbone shape based on figure 10-1a of ASCE 41-17
% 7) If a_hinge and b_hinge values are the same, reduce a by 5% to aid convergence (unless a is already zero)

%% Begin Method
% reduce a_hinge by 5% if it is the same as b_hinge
if (hinge_props.b_hinge-hinge_props.a_hinge) < 0.05*hinge_props.a_hinge
    hinge_props.a_hinge = 0.95*hinge_props.a_hinge;
end

% make sure a_hinge and b_hinge has some minor value
a_hinge = max([1e-6, hinge_props.a_hinge]);
b_hinge = max([2e-6, hinge_props.b_hinge]);

% Define stiffness of lumped plasticisty model based on Ibarra 2005
k_mem = 6*e*iz/length;
if isnan(n) % Ignore Amplification
    k_ele = k_mem; 
    k_spring = k_ele;
else
    k_ele = ((n+1)/n)*k_mem; 
    k_spring = n*k_ele;
end
% Select which stiffness to use for this backbone
if strcmp(type,'full')
    k_elastic = k_mem;
elseif strcmp(type,'hinge')
    k_elastic = k_spring;
end

% Define Post Yeild Strain Hardening Slope
post_yeild_slope_pos = min([(Mp_pos-Mn_pos)/a_hinge,strain_harden_ratio*k_mem]);
post_yeild_slope_neg = min([(Mp_neg-Mn_neg)/a_hinge,strain_harden_ratio*k_mem]);

% Define Backbone Curves
moment_vec_pos = [Mn_pos, post_yeild_slope_pos*a_hinge+Mn_pos,hinge_props.c_hinge*Mn_pos];
moment_vec_neg = [Mn_neg, post_yeild_slope_neg*a_hinge+Mn_neg,hinge_props.c_hinge*Mn_neg];
rot_vec_pos = Mn_pos/k_elastic + [0, a_hinge, b_hinge];
rot_vec_neg = Mn_neg/k_elastic + [0, a_hinge, b_hinge];

end

