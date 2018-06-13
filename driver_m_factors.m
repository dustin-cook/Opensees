clear
close
clc

%% Import Packages
import asce_41.*

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 1;
analysis.name = 'test';

%% Read in element and hinge data tables
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
element = readtable([output_dir filesep 'element.csv'],'ReadVariableNames',true);
m_table = readtable(['+asce_41' filesep 'linear_col_m.csv'],'ReadVariableNames',true);
m_table.id = []; % Omit id 

%% Go through each element and calculate the hinge properties
% CURRENTLY JUST FOR COLUMS, UPDATE TO BE FOR ALL MEMBERS
for i = 1:length(element.id)
    ele = element(i,:);
    ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
    
    % Primary or secondary component
    m_factor = m_table(strcmp(m_table.comp_type, 'primary'),:); % ASSUME ONLY PRIMARY FOR NOW
    
    % Calculate condition
    condition = 1; % UPDATE THIS TO BE DYNAMIC
    m_factor = m_factor(m_factor.condition == condition,:);
    
    % Fitler Table based on P/Asfc
    p_ratio = ele.Pmax/(ele_props.a*ele_props.fc_e);
    [ m_factor ] = fn_filter_asce41_table( m_factor, p_ratio, 'p_ratio', {'m_io','m_ls','m_cp'} );
    
    % Filter table based on row
    row = ele_props.Av/(ele_props.w*ele_props.S);
    [ m_factor ] = fn_filter_asce41_table( m_factor, row, 'row', {'m_io','m_ls','m_cp'} );
    
    % Filter table based on Vye/VcolOE
    v_ratio = ele.Vmax/ele.V0; % NEED TO FIGURE OUT WHAT VcolOE is and UPDATE THIS
    [ m_factor ] = fn_filter_asce41_table( m_factor, v_ratio, 'v_ratio', {'m_io','m_ls','m_cp'} );
    
    % Double Check only 1 row of the table remains
    if length(m_factor.row) ~= 1
        error('Table filtering failed to find unique result')
    end
    
    % save to elements table
    element.m_io(i) = m_factor.m_io;
    element.m_ls(i) = m_factor.m_ls;
    element.m_cp(i) = m_factor.m_cp;
end

%% Save capacities to element database
writetable(element,[output_dir filesep 'element.csv'])