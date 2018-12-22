function [ ele ] = fn_element_critical_mode( ele, ele_prop )
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
    % Start by checking the aspect ratio criteria as defiend by the
    % appendix of ASCE 41-17 (A10.7)
    aspect_ratio = ele.length/ele_prop.d;
    
    if aspect_ratio > 3
        ele.critical_mode = {'flexure'};
    elseif aspect_ratio < 1.5
        ele.critical_mode = {'shear'};
    else
        % Based on Shear Span Ratio (as defined by Pugh, Lowes, Lehman)
        shear_span_ratio = ele.Mmax/(ele.Vmax*ele_prop.d);
        if shear_span_ratio > 2
            ele.critical_mode = {'flexure'};
        else
            ele.critical_mode = {'shear'};
        end
    end
    
    % Save shear at moment yeild as Vye
    ele.vye = NaN;
end

% Check if shear deformations need to be considered
flexural_delta = (1*ele.length^3)/(12*ele_prop.e*ele_prop.iz); % assume a unit lateral force
shear_delta = 1*ele.length/(ele_prop.g*ele_prop.av); % assume a unit lateral force

if shear_delta > 0.1*flexural_delta
    ele.model_shear_deform = true;
else
    ele.model_shear_deform = false;
end



end

