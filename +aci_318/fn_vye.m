function [ vye ] = fn_vye( ele_type, Mn_pos, Mn_neg, ele_length, gravity_load )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

%% Calculate the Shear at flexural yield, Vye, according to ACI 318
if strcmp(ele_type,'beam')
    vye = (Mn_pos + Mn_neg)/ele_length;
elseif strcmp(ele_type,'column')
    vye = (Mn_pos + Mn_neg)/ele_length + gravity_load/2;
elseif strcmp(ele_type,'wall')
    vye = NaN;
end

end

