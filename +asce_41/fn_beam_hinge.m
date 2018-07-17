function [ hinge ] = fn_beam_hinge( ele, ele_props  )
% Description: Find beam hinge properties based on Table 10-7 of ASCE 41-17
% Created by: Dustin Cook
% Date Created: 7-1-18

%% Import Packages
import asce_41.*

% Load Beam Hinge Table 10-7 from ASCE 41-17
hinge_table = readtable(['+asce_41' filesep 'beam_hinge.csv'],'ReadVariableNames',true);
hinge_table.id = []; % Omit id 

%% Calculate condition
condition = 1; % UPDATE THIS TO READ FROM TABLE 10-11
hinge = hinge_table(hinge_table.condition == condition,:);

%% Fitler Table based on Transverse Rienforcement 
if ele_props.S <= ele_props.d/3
    if ele.DCR_total_raw < 2
        % Conforming Transverse Reinforcement
        trans_rien = 'C';
    elseif ele.Vs_aci > 0.75*ele.Vmax
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
hinge = hinge(strcmp(hinge.trans_rien,trans_rien),:);

%% Filter table based on row ratio
row_ratio = 0.5; % Assumptions for now, NEED TO UPDATE
[ hinge ] = fn_filter_asce41_table( hinge, row_ratio, 'row_ratio', {'a_hinge','b_hinge','c_hinge','io','ls','cp'} );

%% Filter table based on V/bd*sqrt(fc)
if sum(isnan(hinge.v_ratio)) == 0
    v_ratio = ele.Vmax/(ele_props.w*ele_props.d*sqrt(ele_props.fc_e));
    [ hinge ] = fn_filter_asce41_table( hinge, v_ratio, 'v_ratio', {'a_hinge','b_hinge','c_hinge','io','ls','cp'} );
end

%% Double Check only 1 row of the hinge table remains
if length(hinge.a_hinge) ~= 1
    error('Hinge table filtering failed to find unique result')
end

end

