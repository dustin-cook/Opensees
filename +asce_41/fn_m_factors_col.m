function [ ele ] = fn_m_factors_col( m_factor_table, ele, ele_props )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import asce_41.*

%% Find M factors 
% Calculate condition
if ele.pass_aci_dev_length == 1
    condition = 1; % Not Controlled by development length
else
    condition = [1,2]; % Controlled by development length (run both to find min)
end

for i = 1:length(condition)
    m_factor = m_factor_table(m_factor_table.condition == condition(i),:);

    % Fitler Table based on P/Asfc
    p_ratio = ele.Pmax/(ele_props.a*ele_props.fc_e);
    [ m_factor ] = fn_filter_asce41_table( m_factor, p_ratio, 'p_ratio', {'m_io','m_ls','m_cp'} );

    % Filter table based on row
    row_1 = ele_props.Av_1/(ele_props.w*ele_props.S_1);
    [ m_factor_1 ] = fn_filter_asce41_table( m_factor, row_1, 'row', {'m_io','m_ls','m_cp'} );
    row_2 = ele_props.Av_2/(ele_props.w*ele_props.S_2);
    [ m_factor_2 ] = fn_filter_asce41_table( m_factor, row_2, 'row', {'m_io','m_ls','m_cp'} );

    % Filter table based on Vye/VcolOE
    v_ratio_1 = ele.Vn_1/ele.V0_1;
    [ m_factor_1 ] = fn_filter_asce41_table( m_factor_1, v_ratio_1, 'v_ratio', {'m_io','m_ls','m_cp'} );
    v_ratio_2 = ele.Vn_2/ele.V0_2;
    [ m_factor_2 ] = fn_filter_asce41_table( m_factor_2, v_ratio_2, 'v_ratio', {'m_io','m_ls','m_cp'} );

    % Double Check only 1 row of the table remains
    if length(m_factor_1{:,1}) ~= 1 || length(m_factor_2{:,1}) ~= 1
        error('Table filtering failed to find unique result')
    end

    % save to elements table
    m_factor_temp_1.m_io(i) = m_factor_1.m_io;
    m_factor_temp_1.m_ls(i) = m_factor_1.m_ls;
    m_factor_temp_1.m_cp(i) = m_factor_1.m_cp;
    m_factor_temp_2.m_io(i) = m_factor_2.m_io;
    m_factor_temp_2.m_ls(i) = m_factor_2.m_ls;
    m_factor_temp_2.m_cp(i) = m_factor_2.m_cp;
end

% Find Min values of cases run
ele.m_io_1 = min(m_factor_temp_1.m_io);
ele.m_ls_1 = min(m_factor_temp_1.m_ls);
ele.m_cp_1 = min(m_factor_temp_1.m_cp);
ele.m_io_2 = min(m_factor_temp_2.m_io);
ele.m_ls_2 = min(m_factor_temp_2.m_ls);
ele.m_cp_2 = min(m_factor_temp_2.m_cp);
end

