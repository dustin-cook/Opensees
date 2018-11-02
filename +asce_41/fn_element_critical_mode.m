function [ ele ] = fn_element_critical_mode( ele, d )
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

%% Determine Controlling Factors
if strcmp(ele.type,'beam')
    % Method 2 - From ACI / Based on Stiffness Matrix
    shear_at_flexure_yeild = (ele.Mn_aci_pos + ele.Mn_aci_neg)/ele.length + ele.gravity_load/2;
    if ele.Vu_aci > shear_at_flexure_yeild
        ele.critical_mode = {'flexure'};
    else
        ele.critical_mode = {'shear'};
    end
    
    % Save shear at moment yeild as Vye
    ele.vye = shear_at_flexure_yeild;
elseif strcmp(ele.type,'column')
    % Method 2 - From ACI / Based on Stiffness Matrix
    shear_at_flexure_yeild = (ele.Mn_aci_pos + ele.Mn_aci_neg)/ele.length;
    if ele.Vu_aci > shear_at_flexure_yeild
        ele.critical_mode = {'flexure'};
    else
        ele.critical_mode = {'shear'};
    end
    
    % Save shear at moment yeild as Vye
    ele.vye = shear_at_flexure_yeild;
elseif strcmp(ele.type,'wall')
    % Method 3 - Based on Shear Span Ratio (as defined by Pugh, Lowes, Lehman)
    shear_span_ratio = ele.Mmax/(ele.Vmax*d);
    if shear_span_ratio > 2
        ele.critical_mode = {'flexure'};
    else
        ele.critical_mode = {'shear'};
    end
    
    % Save shear at moment yeild as Vye
    ele.vye = NaN;
end

% %% Method 4 - Based on Table 10-11 from ASCE 41-13 for Columns
% shear_ratio = shear_at_flexure_yeild/ele.V0;
% if shear_ratio <= 0.6
%     mode_4 = {'flexure'};
% elseif shear_ratio <= 1.0
%     mode_4 = {'flexure-shear'};
% else
%     mode_4 = {'shear'};
% end

end

