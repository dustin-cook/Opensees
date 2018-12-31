function [ ele ] = fn_m_factors_beams( m_factor_table, ele, ele_props )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import asce_41.*

%% Find M factors 
% Calculate condition
if strcmp(ele.critical_mode,'flexure')
    condition(1) = 1; 
elseif strcmp(ele.critical_mode,'shear')
    condition(1) = 2; 
end
if ele.pass_aci_dev_length == 0 % Controlled by inadequate development (assuming no embedment issues)
    condition(2) = 3; 
end


for i = 1:length(condition)
    %% Filter Based on Condition
    m_factor = m_factor_table(m_factor_table.condition == condition(i),:);
    
    if condition(i) == 1
        % Filter table for transverse reinforcement
        if ele_props.S <= ele_props.d/3
            if ele.DCR_total_raw < 2
                % Conforming Transverse Reinforcement
                trans_rien = 'C';
            elseif ele.Vs > 0.75*ele.Vmax
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
        m_factor = m_factor(strcmp(m_factor.tranv_rein, trans_rien),:);

        % Fitler Table based on row-row'/row_bal
        row_ratio = (ele_props.row - ele_props.row_prime) / ele.row_bal;
        [ m_factor ] = fn_filter_asce41_table( m_factor, row_ratio, 'row_ratio', {'m_io','m_ls','m_cp'} );

        % Filter table based on V/bw*d*sqrt(f'ce)
        v_ratio = ele.Vmax/(ele_props.w*ele_props.d*sqrt(ele_props.fc_e)); % NEED TO FIGURE OUT WHAT VcolOE is and UPDATE THIS
        [ m_factor ] = fn_filter_asce41_table( m_factor, v_ratio, 'v_ratio', {'m_io','m_ls','m_cp'} );
    elseif condition(i) == 2 || condition(i) == 3
        %% Filter table based on S
        [ m_factor ] = fn_filter_asce41_table( m_factor, ele_props.S, 's', {'a_hinge','b_hinge','c_hinge','io','ls','cp'} );
    end

    % Double Check only 1 row of the table remains
    if length(m_factor{:,1}) ~= 1
        error('Table filtering failed to find unique result')
    end
    
    %% Save filtered hinge to final table structure
    m_factor_temp(i,:) = m_factor;
end

%% If multiple conditions existed pick the smallest values
if length(condition) > 1
    ele.m_io = min(m_factor_temp.m_io);
    ele.m_ls = min(m_factor_temp.m_ls);
    ele.m_cp = min(m_factor_temp.m_cp);
else
    ele.m_io = m_factor.m_io;
    ele.m_ls = m_factor.m_ls;
    ele.m_cp = m_factor.m_cp;
end

end

