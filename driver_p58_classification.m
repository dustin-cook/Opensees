%% FN to calculate P58 fragility selection inputs
clear all
close all
clc

%% User Inputs
analysis.model_id = 12;
analysis.name = 'output_fix_polly';

%% Load in element table
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
mf_joint_table = readtable(['outputs' filesep model.name{1} filesep 'model data' filesep 'mf_joint.csv'],'ReadVariableNames',true);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
plot_dir = [output_dir filesep 'plots'];
load([output_dir filesep 'node_analysis.mat'])
load([output_dir filesep 'element_analysis.mat'])

%% Go through each element and calculate variables
output = [];
for i = 1:height(mf_joint_table)
    mf = mf_joint_table(i,:);
    output.id(i,1) = i;
    
    % Collect each member properties
    column_low = [];
    column_high = [];
    beam_left = [];
    beam_right = [];
    if mf.column_low > 0
        column_low = element(element.id == mf.column_low,:);
    end
    if mf.column_high > 0
        column_high = element(element.id == mf.column_high,:);
    end
    if mf.beam_left > 0
        beam_left = element(element.id == mf.beam_left,:);
    end
    if mf.beam_right > 0
        beam_right = element(element.id == mf.beam_right,:);
    end
    
    %?Mc/?Mb
    if isempty(column_high)
        sum_col = column_low.Mn_aci_pos;
    else
        sum_col = column_low.Mn_aci_pos + column_high.Mn_aci_pos;
    end
    if isempty(beam_left)
        sum_bm = max([beam_right.Mn_aci_pos,beam_right.Mn_aci_neg]);
    elseif isempty(beam_right)
        sum_bm = max([beam_left.Mn_aci_pos,beam_left.Mn_aci_neg]);
    else
        sum_bm = max([beam_left.Mn_aci_pos + beam_right.Mn_aci_neg,beam_left.Mn_aci_neg + beam_right.Mn_aci_pos]);
    end
    output.scwb(i,1) = sum_col/sum_bm;
    
    % Joint V/Vn
    output.joint_v(i,1) = 0; % Assume adequeate joints for now

    % Inadequate Joint Transverse Reinforcment	
    output.joint_trans_rein(i,1) = 0; % Assume adequeate joints for now

    % ASCE 41-06 High Beam Ductility
    output.beam_duct(i,1) = 1; % Assume flexure control and let later logic update
    if ~isempty(beam_left) && strcmp(beam_left.critical_mode,'shear')
        output.beam_duct(i,1) = 0;
    end
    if ~isempty(beam_right) && strcmp(beam_right.critical_mode,'shear') 
        output.beam_duct(i,1) = 0;
    end
    
    % Inadequate Beam Transverse Reinforcment
    output.beam_trans_rein(i,1) = 0; % Assume adequate transverse for now but need to update (I think others with control)

    % Beam V/bwd(f'c)^0.5	
    beam_vc = 0;
    if ~isempty(beam_left)
        ele_prop = ele_prop_table(ele_prop_table.id == beam_left.ele_id,:);
        beam_vc = [beam_vc, beam_left.Vmax/(ele_prop.w*ele_prop.d*sqrt(ele_prop.fc_e))];
    end
    if ~isempty(beam_right)
        ele_prop = ele_prop_table(ele_prop_table.id == beam_right.ele_id,:);
        beam_vc = [beam_vc, beam_right.Vmax/(ele_prop.w*ele_prop.d*sqrt(ele_prop.fc_e))];
    end
    output.beam_vc(i,1) = max(beam_vc);
    
    %Beam Veq/Vn	
    beam_vn = 0;
    if ~isempty(beam_left)
        beam_vn = [beam_vn, beam_left.Vmax/beam_left.Vn_aci];
    end
    if ~isempty(beam_right)
        beam_vn = [beam_vn, beam_right.Vmax/beam_right.Vn_aci];
    end
    output.beam_vn(i,1) = max(beam_vn);
    
    % Inadequate Column Transverse Reinforcment	
    output.col_trans_rein(i,1) = 0; % Assume adequate transverse for now but need to update (I think others with control)
    
    % Pu/0.6f'cAg
    ele_prop = ele_prop_table(ele_prop_table.id == column_low.ele_id,:);
    output.col_pu(i,1) = column_low.P_grav/(0.6*ele_prop.fc_e*ele_prop.a);
    
    % Column Veq/Vn	
    col_vn = 0;
    if ~isempty(column_high)
        col_vn = [col_vn, column_high.Vmax/column_high.Vn_aci];
    end
    col_vn = [col_vn, column_low.Vmax/column_low.Vn_aci];
    output.col_vn(i,1) = max(col_vn);
    
    %Inadequate Reinforcement Development per ASCE 41-06
    output.dev_length(i,1) = ~column_low.pass_aci_dev_length;
    if ~isempty(column_high)
        output.dev_length(i,1) = ~column_high.pass_aci_dev_length;
    end
    if ~isempty(beam_left)
        output.dev_length(i,1) = ~beam_left.pass_aci_dev_length;
    end
    if ~isempty(beam_right)
        output.dev_length(i,1) = ~beam_right.pass_aci_dev_length;
    end
end

%% Save data to csv
output_table = struct2table(output);
writetable(output_table,[output_dir filesep 'fragility_select_data.csv'])