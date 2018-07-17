function [ ele ] = fn_m_factors_col( m_factor, ele, ele_props )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import asce_41.*

%% Find M factors 
% Fitler Table based on P/Asfc
p_ratio = ele.Pmax/(ele_props.a*ele_props.fc_e);
[ m_factor ] = fn_filter_asce41_table( m_factor, p_ratio, 'p_ratio', {'m_io','m_ls','m_cp'} );

% Filter table based on row
row = ele_props.Av/(ele_props.w*ele_props.S);
[ m_factor ] = fn_filter_asce41_table( m_factor, row, 'row', {'m_io','m_ls','m_cp'} );

% Filter table based on Vye/VcolOE
v_ratio = ele.Vmax/ele.V0; % NEED TO FIGURE OUT WHAT Vye is and UPDATE THIS
[ m_factor ] = fn_filter_asce41_table( m_factor, v_ratio, 'v_ratio', {'m_io','m_ls','m_cp'} );

% Double Check only 1 row of the table remains
if length(m_factor{:,1}) ~= 1
    error('Table filtering failed to find unique result')
end

% save to elements table
ele.m_io = m_factor.m_io;
ele.m_ls = m_factor.m_ls;
ele.m_cp = m_factor.m_cp;

end

