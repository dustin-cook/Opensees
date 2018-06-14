%%% Driver to Calculate the Capacities of Elements in model %%%
clear all
close all
clc

%% Import Packages
import asce_41.*

%% Read in element table
analysis.model_id = 3;
analysis.name = '09DL';
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
element = readtable([output_dir filesep 'element.csv'],'ReadVariableNames',true);
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

% First Element
ele = element(1,:);
ele_id = ele.ele_id;
ele_prop = ele_prop_table(ele_id,:);

% Axial Capacity per ACI
[ ~, Pn_max ] = fn_aci_axial_capacity( ele_prop.fc_e, ele_prop.a, ele_prop.As, ele_prop.fy_e );
As = str2double(strsplit(strrep(strrep(ele_prop.As{1},'[',''),']',''),','));
tension_cap = -sum(As)*ele_prop.fy_e;

% Percent Pmax
p_vals = linspace(0.75*tension_cap,0.9*Pn_max,15);

% Go through range of P values and calc moment capacity
for i = 1:length(p_vals)
    p_vals(i)
    [ ~, Mn_aci(i) ] = fn_aci_moment_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.clear_cover, ele_prop.Es, p_vals(i) );
end

plot([0,Mn_aci,0],[tension_cap, p_vals,Pn_max])
grid on
box on
xlabel('Moment, M')
ylabel('Axial Load, P')