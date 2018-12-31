%% Build P-58 analysis Inputs
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
frag_select = readtable([output_dir filesep 'fragility_select_data.csv'],'ReadVariableNames',true);
plot_dir = [output_dir filesep 'plots'];
load([output_dir filesep 'node_analysis.mat'])
load([output_dir filesep 'element_analysis.mat'])
load([output_dir filesep 'story_analysis.mat'])
load([output_dir filesep 'hinge_analysis.mat'])

%% For each element pull together inputs
id = 0;
for i = 1:height(element)
    ele = element(i,:);
    if strcmp(ele.type,'beam') || strcmp(ele.type,'column')
        if ele.story == 1 %|| ele.story == 2
            ele_prop = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
            joint_id = mf_joint_table.id((mf_joint_table.column_low == ele.id) | (mf_joint_table.column_high == ele.id) | (mf_joint_table.beam_left == ele.id) | (mf_joint_table.beam_right == ele.id));
            for j = 1:length(joint_id)
                id = id + 1;
                inputs.id(id,1) = id;
                inputs.element_type{id,1} = ele.type{1};
                inputs.story(id,1) = ele.story;
                inputs.description(id,1) = ele.id;
                inputs.p58_id{id,1} = frag_select.fragility{frag_select.id == joint_id(j)};
                inputs.edp_type{id,1} = 'sdr';
                inputs.max_edp(id,1) = story.max_drift_x(1); % NEED TO UPDATE THIS
                hin = hinge(hinge.element_id == ele.id,:);
                hinge_demand = [];
                for k = 1:height(hin)
                    hinge_demand = [hinge_demand, max(hin.rotation_TH{k})];
                end
                inputs.max_element_demand(id,1) = max(hinge_demand);
                inputs.e(id,1) = ele_prop.e;
                inputs.i(id,1) = ele_prop.iz;
                inputs.l(id,1) = ele.length;
                inputs.m_n(id,1) = max([ele.Mn_pos,ele.Mn_neg]);
                inputs.m_ult(id,1) = max([ele.Mp_pos,ele.Mp_neg]);
                inputs.a_hinge(id,1) = ele.a_hinge;
                inputs.b_hinge(id,1) = ele.b_hinge;
                inputs.c_hinge(id,1) = ele.c_hinge;
                inputs.io(id,1) = ele.io;
                inputs.ls(id,1) = ele.ls;
                inputs.cp(id,1) = ele.cp;
            end
        end
    end
end

%% Save data to csv
input_table = struct2table(inputs);
writetable(input_table,[output_dir filesep 'p58_inputs.csv'])