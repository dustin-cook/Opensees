%%% Driver to Calculate the Capacities of Elements in model %%%
clear all
close all
clc

%% Import Packages
import aci_318.*

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

% Or Manually Set Element Proprs
% ele_prop.w = 14;
% ele_prop.d = 24;
% ele_prop.a = ele_prop.w*ele_prop.d;
% ele_prop.As = {'[3,3]'};
% ele_prop.As_d = {'[2.5,21.5]'};
% ele_prop.fc_n = 4000;
% ele_prop.fc_e = 4000;
% ele_prop.fy_n = 60000;
% ele_prop.fy_e = 60000;
% ele_prop.Es = 29000000;

% Axial Capacity per ACI
[ ~, Pn_aci_c, ~, ~ ] = fn_aci_axial_capacity( ele_prop.fc_n, ele_prop.a, ele_prop.As, ele_prop.fy_n );
[ ~, ~, ~, Pn_aci_t ] = fn_aci_axial_capacity( ele_prop.fc_e, ele_prop.a, ele_prop.As, ele_prop.fy_e );
As = str2double(strsplit(strrep(strrep(ele_prop.As{1},'[',''),']',''),','));
compression_cap_exp = ele_prop.fc_e*(ele_prop.a-sum(As))+ele_prop.fy_e*sum(As);

% Percent Pmax
p_vals_1 = linspace(-0.95*Pn_aci_t,Pn_aci_c,100);
p_vals_2 = linspace(-0.7*Pn_aci_t,0.7*compression_cap_exp,100);
% p_vals = [0,504400,623700,1000000];

% Go through range of P values and calc moment capacity
for i = 1:length(p_vals_1)
    [ ~, Mn_1(i) ] = fn_aci_moment_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, p_vals_1(i) );
end

for i = 1:length(p_vals_2)
    [ ~, Mn_2(i) ] = fn_aci_moment_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.As, ele_prop.As_d, ele_prop.fy_e, ele_prop.Es, p_vals_2(i) );
end

hold on
plot([0,Mn_2/12,0]/1000,[-Pn_aci_t, p_vals_2,compression_cap_exp]/1000,'--r')
plot([0,Mn_1/12,0]/1000,[-Pn_aci_t, p_vals_1,Pn_aci_c]/1000)
grid on
box on
xlabel('Moment, M (kip-ft)')
ylabel('Axial Load, P (kips)')