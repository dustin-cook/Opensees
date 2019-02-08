function [ ele ] = fn_element_critical_mode( ele, ele_prop )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% %% Import Packages
% import aci_318.fn_fye

% % Method 1 - based on actual loading from analysis
% demand = ele.Mmax/ele.Vmax;
% capacity = ele.Mn_pos/ele.Vn;
% if demand > capacity
%     ele.critical_mode = {'flexure'};
% else
%     ele.critical_mode = {'shear'};
% end

%% Calculate Shear at Flexural Yield
% [ ele.vye ] = fn_vye( ele.type, ele.Mn_pos, ele.Mn_neg, ele.length, ele.gravity_load );

%% Determine Controlling Factors
% Method 2 - From ACI / Based on Stiffness Matrix
if strcmp(ele.type,'beam') || strcmp(ele.type,'column')
    if ele.Vn > ele.vye
        ele.critical_mode = {'flexure'};
    else
        ele.critical_mode = {'shear'};
    end
    
    % OOP
    if ele.Vn > ele.vye_oop % Assumes shear strength is the same for OOP
        ele.critical_mode_oop = {'flexure'};
    else
        ele.critical_mode_oop = {'shear'};
    end
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
        if sum(strcmp('Mmax',ele.Properties.VariableNames)) == 0
            shear_span_ratio = 1;
        else
            shear_span_ratio = ele.Mmax/(ele.Vmax*ele_prop.d);
        end
        if shear_span_ratio > 2
            ele.critical_mode = {'flexure'};
        else
            ele.critical_mode = {'shear'};
        end
    end
    
    % OOP
    ele.critical_mode_oop = {'flexure'}; % assumes wall out of plane is flexure controlled
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

