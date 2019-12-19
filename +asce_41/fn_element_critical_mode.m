function [ critical_mode, critical_mode_oop, model_shear_deform ] = fn_element_critical_mode( ele, ele_prop, Vn, vye, vye_oop )
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
    if Vn > vye
        critical_mode = {'flexure'};
    else
        critical_mode = {'shear'};
    end
    
    % OOP
    if Vn > vye_oop % Assumes shear strength is the same for OOP
        critical_mode_oop = {'flexure'};
    else
        critical_mode_oop = {'shear'};
    end
elseif strcmp(ele.type,'wall')
    % Start by checking the aspect ratio criteria as defiend by the
    % appendix of ASCE 41-17 (A10.7)
    aspect_ratio = ele.length/ele_prop.h;
    
    if aspect_ratio > 3
        critical_mode = {'flexure'};
    elseif aspect_ratio < 1.5
        critical_mode = {'shear'};
    else
        % Based on Shear Span Ratio (as defined by Pugh, Lowes, Lehman)
        if sum(strcmp('Mmax',ele.Properties.VariableNames)) == 0
            shear_span_ratio = 1;
        else
            shear_span_ratio = ele.Mmax/(ele.Vmax*ele_prop.h);
        end
        if shear_span_ratio > 2
            critical_mode = {'flexure'};
        else
            critical_mode = {'shear'};
        end
    end
    
    % OOP
    critical_mode_oop = {'flexure'}; % assumes wall out-of-plane are flexure controlled
    
else % Not a beam column or wall
    critical_mode = {'NA'};
    critical_mode_oop = {'NA'}; 
end

% Check if shear deformations need to be considered
flexural_delta = (1*ele.length^3)/(12*ele_prop.e*ele_prop.iz); % assume a unit lateral force
shear_delta = 1*ele.length/(ele_prop.g*ele_prop.av); % assume a unit lateral force

if shear_delta > 0.1*flexural_delta
    model_shear_deform = true;
else
    model_shear_deform = false;
end



end

