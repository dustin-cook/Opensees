function [ ele ] = fn_element_critical_mode( ele )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% % Method 1 - based on actual loading from analysis
% demand = ele.Mmax/ele.Vmax;
% capacity = ele.Mn_aci/ele.Vn_aci;
% if demand > capacity
%     ele.critical_mode = {'flexure'};
% else
%     ele.critical_mode = {'shear'};
% end

%% Method 2 - Based on Stiffness Matrix
shear_at_flexure_yeild = 2*max([ele.Mn_aci_c,ele.Mn_aci_t])/ele.length;
if ele.Vn_aci > shear_at_flexure_yeild
    ele.critical_mode = {'flexure'};
else
    ele.critical_mode = {'shear'};
end

end

