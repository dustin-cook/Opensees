function [ hinge, trans_rien ] = fn_beam_hinge( ele, ele_props, ele_side )
% Description: Find beam hinge properties based on Table 10-7 of ASCE 41-17
% Created by: Dustin Cook
% Date Created: 7-1-18

%% Import Packages
import asce_41.fn_filter_asce41_table

% Define vars if they do not exist
if sum(strcmp('Vmax',ele.Properties.VariableNames)) == 0
    ele.Vmax = 0;
end

% Load Beam Hinge Table 10-7 from ASCE 41-17
hinge_table = readtable(['+asce_41' filesep 'beam_hinge.csv'],'ReadVariableNames',true);
hinge_table.id = []; % Omit id 

%% Calculate condition
if strcmp(ele.(['critical_mode_' num2str(ele_side)]),'flexure')
    condition(1) = 1; 
elseif strcmp(ele.(['critical_mode_' num2str(ele_side)]),'shear')
    condition(1) = 2; 
end
if ele.pass_aci_dev_length == 0 % Controlled by inadequate development (assuming no embedment issues)
    condition(2) = 3; 
end

for i = 1:length(condition)
    %% Filter Based on Condition
    hinge_filt = hinge_table(hinge_table.condition == condition(i),:);
    
    % Calculate if the transverse reinforcement
    if ele_props.(['S_' num2str(ele_side)]) <= ele_props.d_eff/3
        if ele.(['disp_duct_' num2str(ele_side)]) < 2
            % Conforming Transverse Reinforcement
            trans_rien = 'C';
        elseif ele.(['Vs_' num2str(ele_side)]) > 0.75*ele.Vmax
            % Conforming Transverse Reinforcement
            trans_rien = 'C'; 
        else
            % Non-conforming Transverse Reinforcement
            trans_rien = 'NC';
        end
    else
        % Non-conforming Transverse Reinforcement
        trans_rien = 'NC';
    end
    if condition(i) == 1 
        %% Fitler Table based on Transverse Rienforcement 
        hinge_filt = hinge_filt(strcmp(hinge_filt.trans_rien,trans_rien),:);

        %% Filter table based on row ratio
        row_ratio = (ele_props.row - ele_props.row_prime) / ele.row_bal;
        [ hinge_filt ] = fn_filter_asce41_table( hinge_filt, row_ratio, 'row_ratio', {'a_hinge','b_hinge','c_hinge','io','ls','cp'} );

        %% Filter table based on V/bd*sqrt(fc)
        if sum(isnan(hinge_filt.v_ratio)) == 0
            v_ratio = ele.Vmax/(ele_props.w*ele_props.d_eff*sqrt(ele_props.fc_e));
            [ hinge_filt ] = fn_filter_asce41_table( hinge_filt, v_ratio, 'v_ratio', {'a_hinge','b_hinge','c_hinge','io','ls','cp'} );
        end
    elseif condition(i) == 2 || condition(i) == 3
        %% Filter table based on S
        [ hinge_filt ] = fn_filter_asce41_table( hinge_filt, ele_props.S, 's', {'a_hinge','b_hinge','c_hinge','io','ls','cp'} );
    end

    %% Double Check only 1 row of the hinge table remains
    if length(hinge_filt.a_hinge) ~= 1
        error('Hinge table filtering failed to find unique result')
    end
    
    %% Save filtered hinge to final table structure
    hinge_temp(i,:) = hinge_filt;
end

%% If multiple conditions existed pick the smallest values
if length(condition) > 1
    hinge.a_hinge = min(hinge_temp.a_hinge);
    hinge.b_hinge = min(hinge_temp.b_hinge);
    hinge.c_hinge = min(hinge_temp.c_hinge);
    hinge.io = min(hinge_temp.io);
    hinge.ls = min(hinge_temp.ls);
    hinge.cp = min(hinge_temp.cp);
    hinge = struct2table(hinge);
else
    hinge = hinge_temp(:,6:11);
end
    

end

