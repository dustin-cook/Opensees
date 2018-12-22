function [ force_vec, disp_vec ] = fn_define_backbone_shear( Vn, length, g, av, hinge_props )
% Takes points from ASCE 41-17 chapter 10 tables and Define backbone curve
% as a vector of displacements and forces for both the positive and
% negative directions. Units are based on the basic units that are input.

% Created By: Dustin Cook
% Data Created: December 20, 2018

% Inputs
% 

% Outputs
% force_vec: Vector of shear force points along the backbone curve for the
%               positive or negative direction. Does not include the
%               origin.
% disp_vec: Vector of lateral displacement points along the backbone curve for the
%               positive or negative direction. Does not include the
%               origin.

% Assumptions:
% 1) Backbone shape based on figure 10-1c of ASCE 41-17
% 2) Elastic stiffness derived from delta = PL/GA
% 3) Shear backbone is the same in both directions
% 4) Post peak slope goes straight from point C to E (peak to end)
% 5) If a_hinge and b_hinge values are the same, reduce a by 5% to aid convergence (unless a is already zero)

%% Begin Method
% Define Elastic Stiffness
elastic_stiffness = g*av/length;

% Define Force vector
f1 = hinge_props.f_hinge*Vn;
f2 = Vn;
f3 = Vn + 1e-6; % add neglibable amount to aid with interpolation
f4 = Vn*hinge_props.c_hinge;
force_vec = [f1,f2,f3,f4];

% Define Displacement Vector
u1 = f1/elastic_stiffness; 
u2 = (hinge_props.g_hinge/100)*length; % change to percent and displacemnt
u3 = (hinge_props.d_hinge/100)*length;
u4 = (hinge_props.e_hinge/100)*length;
disp_vec = [u1,u2,u3,u4];

end

