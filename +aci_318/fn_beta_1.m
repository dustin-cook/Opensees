function [ beta_1 ] = fn_beta_1( fc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Assumes F'c is in psi

beta1_min = 0.65;
beta1_max = 0.85;
beta_1_interp = 0.85 - (fc - 4000)*(0.05/1000);

beta_1 = min([max([beta_1_interp,beta1_min]),beta1_max]);

end

