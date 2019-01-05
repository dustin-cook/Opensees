function [ ele ] = fn_m_factors( ele, ele_props )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
% Import Packages
import asce_41.*

% Load M Factor Tables
m_table.col = readtable(['+asce_41' filesep 'linear_col_m.csv'],'ReadVariableNames',true);
m_table.beam = readtable(['+asce_41' filesep 'linear_beam_m.csv'],'ReadVariableNames',true);

%% Find M factors 
% Omit ID on m-factors table
m_table.col.id = [];
m_table.beam.id = [];

% Preallocate table
ele.m_io = 1;
ele.m_ls = 1;
ele.m_cp = 1;

% Primary or secondary component
m_table.col = m_table.col(strcmp(m_table.col.comp_type, 'primary'),:); % ASSUME ONLY PRIMARY FOR NOW
m_table.beam = m_table.beam(strcmp(m_table.beam.comp_type, 'primary'),:); % ASSUME ONLY PRIMARY FOR NOW

if strcmp(ele_props.type,'column')
    % M factors for Columns
    [ ele ] = fn_m_factors_col( m_table.col, ele, ele_props );
elseif strcmp(ele_props.type,'beam')
    % M factors for Beams
    [ ele ] = fn_m_factors_beams( m_table.beam, ele, ele_props );
end



end

