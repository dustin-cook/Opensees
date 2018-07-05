clear
close
clc

%% Import Packages
import asce_41.*
import plotting_tools.*

%% Define Analysis and Model parameters
analysis.model_id = 3;
analysis.gm_id = 1;
analysis.name = 'NL_10DL10LL';

%% Read in element and hinge data tables
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs' filesep model.name{1} filesep analysis.name];
ele_prop_table = readtable(['inputs' filesep 'element.csv'],'ReadVariableNames',true);
element = readtable([output_dir filesep 'element_linear.csv'],'ReadVariableNames',true);
hinge_table = readtable(['+asce_41' filesep 'col_hinge.csv'],'ReadVariableNames',true);
hinge_table.id = []; % Omit id 

%% Go through each element and calculate the hinge properties
for i = 1:length(element.id)
    ele = element(i,:);
    ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
    
    % Calculate condition
    condition = 2; % UPDATE THIS TO READ FROM TABLE 10-11
    hinge = hinge_table(hinge_table.condition == condition,:);
    
    % Fitler Table based on P/Asfc
    p_ratio = ele.Pmax/(ele_props.a*ele_props.fc_e);
    if p_ratio<=min(hinge.p_ratio)
        hinge = hinge(hinge.p_ratio == min(hinge.p_ratio),:);
    elseif p_ratio>=max(hinge.p_ratio)
        hinge = hinge(hinge.p_ratio == max(hinge.p_ratio),:);
    else
        hinge_1 = sortrows(hinge(hinge.p_ratio == max(hinge.p_ratio),:));
        hinge_2 = sortrows(hinge(hinge.p_ratio == min(hinge.p_ratio),:));
        hinge = hinge_1;
        hinge.a_hinge = ((hinge_1.a_hinge-hinge_2.a_hinge)/(max(hinge_1.p_ratio)-min(hinge_2.p_ratio)))*(p_ratio-min(hinge_2.p_ratio)) + hinge_2.a_hinge;
        hinge.b_hinge = ((hinge_1.b_hinge-hinge_2.b_hinge)/(max(hinge_1.p_ratio)-min(hinge_2.p_ratio)))*(p_ratio-min(hinge_2.p_ratio)) + hinge_2.b_hinge;
        hinge.c_hinge = ((hinge_1.c_hinge-hinge_2.c_hinge)/(max(hinge_1.p_ratio)-min(hinge_2.p_ratio)))*(p_ratio-min(hinge_2.p_ratio)) + hinge_2.c_hinge;
        hinge.io = ((hinge_1.io-hinge_2.io)/(max(hinge_1.p_ratio)-min(hinge_2.p_ratio)))*(p_ratio-min(hinge_2.p_ratio)) + hinge_2.io;
        hinge.ls = ((hinge_1.ls-hinge_2.ls)/(max(hinge_1.p_ratio)-min(hinge_2.p_ratio)))*(p_ratio-min(hinge_2.p_ratio)) + hinge_2.ls;
        hinge.cp = ((hinge_1.cp-hinge_2.cp)/(max(hinge_1.p_ratio)-min(hinge_2.p_ratio)))*(p_ratio-min(hinge_2.p_ratio)) + hinge_2.cp;
        hinge.p_ratio(:) = p_ratio;
    end
    
    % Filter table based on row
    row = ele_props.Av/(ele_props.w*ele_props.S);
    if row<=min(hinge.row)
        hinge = hinge(hinge.row == min(hinge.row),:);
    elseif row>=max(hinge.row)
        hinge = hinge(hinge.row == max(hinge.row),:);
    else
        hinge_1 = sortrows(hinge(hinge.row == max(hinge.row),:));
        hinge_2 = sortrows(hinge(hinge.row == min(hinge.row),:));
        hinge = hinge_1;
        hinge.a_hinge = ((hinge_1.a_hinge-hinge_2.a_hinge)/(max(hinge_1.row)-min(hinge_2.row)))*(row-min(hinge_2.row)) + hinge_2.a_hinge;
        hinge.b_hinge = ((hinge_1.b_hinge-hinge_2.b_hinge)/(max(hinge_1.row)-min(hinge_2.row)))*(row-min(hinge_2.row)) + hinge_2.b_hinge;
        hinge.c_hinge = ((hinge_1.c_hinge-hinge_2.c_hinge)/(max(hinge_1.row)-min(hinge_2.row)))*(row-min(hinge_2.row)) + hinge_2.c_hinge;
        hinge.io = ((hinge_1.io-hinge_2.io)/(max(hinge_1.row)-min(hinge_2.row)))*(row-min(hinge_2.row)) + hinge_2.io;
        hinge.ls = ((hinge_1.ls-hinge_2.ls)/(max(hinge_1.row)-min(hinge_2.row)))*(row-min(hinge_2.row)) + hinge_2.ls;
        hinge.cp = ((hinge_1.cp-hinge_2.cp)/(max(hinge_1.row)-min(hinge_2.row)))*(row-min(hinge_2.row)) + hinge_2.cp;
        hinge.row(:) = row;
    end
    
    % Filter table based on V/bd*sqrt(fc)
    if sum(isnan(hinge.v_ratio)) == 0
        v_ratio = ele.Vmax/(ele_props.w*ele_props.d*sqrt(ele_props.fc_e));
        if v_ratio<=min(hinge.v_ratio)
            hinge = hinge(hinge.v_ratio == min(hinge.v_ratio),:);
        elseif v_ratio>=max(hinge.v_ratio)
            hinge = hinge(hinge.v_ratio == max(hinge.v_ratio),:);
        else
            hinge_1 = sortrows(hinge(hinge.v_ratio == max(hinge.v_ratio),:));
            hinge_2 = sortrows(hinge(hinge.v_ratio == min(hinge.v_ratio),:));
            hinge = hinge_1;
            hinge.a_hinge = ((hinge_1.a_hinge-hinge_2.a_hinge)/(max(hinge_1.v_ratio)-min(hinge_2.v_ratio)))*(v_ratio-min(hinge_2.v_ratio)) + hinge_2.a_hinge;
            hinge.b_hinge = ((hinge_1.b_hinge-hinge_2.b_hinge)/(max(hinge_1.v_ratio)-min(hinge_2.v_ratio)))*(v_ratio-min(hinge_2.v_ratio)) + hinge_2.b_hinge;
            hinge.c_hinge = ((hinge_1.c_hinge-hinge_2.c_hinge)/(max(hinge_1.v_ratio)-min(hinge_2.v_ratio)))*(v_ratio-min(hinge_2.v_ratio)) + hinge_2.c_hinge;
            hinge.io = ((hinge_1.io-hinge_2.io)/(max(hinge_1.v_ratio)-min(hinge_2.v_ratio)))*(v_ratio-min(hinge_2.v_ratio)) + hinge_2.io;
            hinge.ls = ((hinge_1.ls-hinge_2.ls)/(max(hinge_1.v_ratio)-min(hinge_2.v_ratio)))*(v_ratio-min(hinge_2.v_ratio)) + hinge_2.ls;
            hinge.cp = ((hinge_1.cp-hinge_2.cp)/(max(hinge_1.v_ratio)-min(hinge_2.v_ratio)))*(v_ratio-min(hinge_2.v_ratio)) + hinge_2.cp;
            hinge.v_ratio(:) = v_ratio;
        end
    end
    
    % Double Check only 1 row of the hinge table remains
    if length(hinge.a_hinge) ~= 1
        error('Hinge table filtering failed to find unique result')
    end
    
%     % Plot Hinges
%     theta_yeild = ele.Mn_aci*(ele_props.d/2)/(ele_props.e*ele_props.iz);
%     Q_y = ele.Mn_aci;
%     Q_ult = ele.Mp;
%     post_yeild_strength = min(ele.Mp,1.1);
%     force_vector = [0,1,post_yeild_strength,hinge.c_hinge,hinge.c_hinge];
%     disp_vector = [0, theta_yeild, theta_yeild+hinge.a_hinge, theta_yeild+hinge.a_hinge+(hinge.b_hinge-hinge.a_hinge)/2, theta_yeild+hinge.b_hinge]; % ASSUMING y = d/2 NEED TO UPDATE
%     plot(disp_vector,force_vector)
%     ylabel('Q/Qy')
%     xlabel('Total Rotation')
%     fn_format_and_save_plot( [output_dir filesep 'hinge_plots' filesep] , ['element_' num2str(ele.id)], 2 )

    % save as element hinge table
    element.condition(i) = hinge.condition;
    element.p_ratio(i) = hinge.p_ratio;
    element.row(i) = hinge.v_ratio;
    element.v_ratio(i) = hinge.v_ratio;
    
    element.a_hinge(i) = hinge.a_hinge;
    element.b_hinge(i) = hinge.b_hinge;
    element.c_hinge(i) = hinge.c_hinge;
    element.io(i) = hinge.io;
    element.ls(i) = hinge.ls;
    element.cp(i) = hinge.cp;
end

%% Save capacities to element database
writetable(element,[output_dir filesep 'element_linear.csv'])