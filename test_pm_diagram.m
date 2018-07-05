%%% Driver to Calculate the Capacities of Elements in model %%%
clear all
close all
clc

%% Import Packages
import asce_41.*

%% Read in element table
analysis.model_id = 4;
analysis.name = '11DL11LL';
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
element = readtable([output_dir filesep 'element.csv'],'ReadVariableNames',true);
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);

% Pick Element
element_id = 6;
ele = element(element_id,:);
ele_id = ele.ele_id;
ele_prop = ele_prop_table(ele_id,:);

% Axial Capacity per ACI
[ ~, Pn_max ] = fn_aci_axial_capacity( ele_prop.fc_n, ele_prop.a, ele_prop.As, ele_prop.fy_n );
As = str2double(strsplit(strrep(strrep(ele_prop.As{1},'[',''),']',''),','));
tension_cap = -sum(As)*ele_prop.fy_n;

% Percent Pmax
p_vals = linspace(0.7*tension_cap,0.99*Pn_max,20);

% Go through range of P values and calc moment capacity
for i = 1:length(p_vals)
    p_vals(i)
    [ ~, Mn_aci(i) ] = fn_aci_moment_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.clear_cover, ele_prop.Es, p_vals(i) );
end

plot([0,Mn_aci,0]/1000,[tension_cap, p_vals,Pn_max]/1000)
grid on
box on
xlabel('Moment, M (kip-in)')
ylabel('Axial Load, P (kips)')