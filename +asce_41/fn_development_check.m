function [ pass_aci_dev_length, ld_avail, ld_req ] = fn_development_check( ele, ele_prop )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Assumptions
% 1. Bars are not coated
% 2. Bars are deformed and not plain
% 3. C_d is just the distance to the edge
% 4. Ignores Hook clause of 25.4.3.3 (when cover is less than 2.5 in)
% 5. Ignores 108 deg hooks
% 6. Hooks in tension do not have high detailed confinement (table 25.4.3.2)
% 7. ASCE 41-17 eq 10-1a will always equal Fy
% 8. Splices do not occur in excpected inelastic zones
% 9. DCR as a placeholder for Mmax/My and Vmax/Vy calc
% 10. Embedment is not an issue.

%% Import Packages
import aci_318.*

%% Maniluptae String Vector Properties to Arrays
As_d = str2double(strsplit(strrep(strrep(ele_prop.As_d{1},'[',''),']',''),','));
d_b = str2double(strsplit(strrep(strrep(ele_prop.d_b{1},'[',''),']',''),','));
l_d_min = str2double(strsplit(strrep(strrep(ele_prop.l_d_min{1},'[',''),']',''),','));
l_d_hook_min = str2double(strsplit(strrep(strrep(ele_prop.l_d_hook_min{1},'[',''),']',''),','));
l_d_hook_ext_min = str2double(strsplit(strrep(strrep(ele_prop.l_d_hook_ext_min{1},'[',''),']',''),','));
l_s_min = str2double(strsplit(strrep(strrep(ele_prop.l_s_min{1},'[',''),']',''),','));

for i = 1:length(d_b)
    %% Calculate Development Length Factors According to ACI 318-14
    c_b = ele_prop.clear_cover; % Just assume distance to edge in min for now

    % Table 25.4.2.4
    psi_e = 1.0; % Assume Non Coated for now
    if d_b(i) >= 0.875
        psi_s = 1.0; % for no. 7 bars or greater
    else
        psi_s = 0.8; % for no. 6 bars or less
    end
    if strcmp(ele_prop.type,'beam') && (ele_prop.d - As_d(i)) > 12
        psi_t = 1.3; % for top bars in members deeper than 12 in
    else
        psi_t = 1.0;
    end

    % Table 25.4.3.2
    if d_b(i) <= 1.41 && ele_prop.clear_cover >= 2.5
        psi_c = 0.7; % for enough cover on bar
    else
        psi_c = 1.0; % not enough cover on bar
    end
    psi_r = 1.0; % Assume not enough confinement reinforcing for now

    % Table 25.4.9.3
    if ele_prop.S <= 4
        psi_rc = 0.75; % transverse reinformcent spacing less than 4"
    else
        psi_rc = 1.0; % transverse reinformcent spacing greater than 4"
    end

    %% Calculate Development Length According to ACI 318-14
    [ l_d(i), l_dt, l_dt_raw, l_dht ] = fn_aci_development_length( ele_prop.fy_e, ele_prop.fc_e, ele_prop.Av, ele_prop.S, d_b(i), c_b, ele_prop.lambda, psi_t, psi_e, psi_s, psi_c, psi_r, psi_rc );

    %% Calculate splice length
    [ l_s ] = fn_aci_splice_length( ele_prop.fy_e, ele_prop.fc_e, d_b(i), l_dt_raw );

    %% Check to see if development and embedment lengths are long enough
    if l_d_min == 999
        test1 = 1;
    else
        test1 = l_d_min(i) >= l_d(i);
    end
    if l_d_hook_min == 999
        test2 = 1;
    else
        test2 = l_d_hook_min(i) >= l_dht;
    end
    if l_d_hook_ext_min == 999
        test3 = 1;
    else
        test3 = l_d_hook_ext_min(i) >= 12*d_b(i);
    end
    if l_s_min == 999
        test4 = 1;
    else
        test4 = l_s_min(i) >= min([l_s,l_dt]); % Per 10.3.5 of ASCE 41-17
    end
    if test1 && test2 && test3 && test4
        pass_aci_dev_length(i) = 1; 
    else
        pass_aci_dev_length(i) = 0;
    end
    
    % Check Hooked Embedments on First Story Columns
    if strcmp(ele.type,'column') && ele.story == 1
        if ~(test2 && test3)
            error('Columns are not embeded enough into foundations, rework foundation model')
        end
    end
end

pass_aci_dev_length = min(pass_aci_dev_length);

% ASCE 41-17 10.3.5 point 2 check on if max applied stress exceeds Fs
if sum(strcmp('DCR_raw_max_V',ele.Properties.VariableNames)) == 1 && ele.DCR_raw_max_V < 1.0 % Should replace with My and Vy check
    pass_aci_dev_length = 1; % Not enough demand on memeber to be controlled by inadequate development length
end

% Save development length info for ATC 134 P-58 data collection
[~, idx] = min(l_d_min - l_d);
ld_avail = l_d_min(idx); % Use Min diff between avail and requred as representative value
ld_req = l_d(idx);

end

