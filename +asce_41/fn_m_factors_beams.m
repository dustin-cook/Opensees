function [ ele ] = fn_m_factors_beams( m_factor, ele, ele_props )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import asce_41.*

%% Find M factors 
% Filter table for transverse reinforcement
tranv_rein = 'C'; % ASSUME CONFORMING FOR NOW (update later)
m_factor = m_factor(strcmp(m_factor.tranv_rein, tranv_rein),:);

% Fitler Table based on row-row'/row_bal
row_ratio = 0; % ASSUME EQUAL AMOUNT OF TENSION AND COMPRESSION REINFORCEMENT FOR NOW
[ m_factor ] = fn_filter_asce41_table( m_factor, row_ratio, 'row_ratio', {'m_io','m_ls','m_cp'} );

% Filter table based on V/bw*d*sqrt(f'ce)
v_ratio = ele.Vmax/(ele_props.w*ele_props.d*sqrt(ele_props.fc_e)); % NEED TO FIGURE OUT WHAT VcolOE is and UPDATE THIS
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

