clear
close
clc

%% Import Packages
import asce_41.*

%% Read in element table
element_table = readtable(['element.csv'],'ReadVariableNames',true);
hinge_table = readtable(['+asce_41' filesep 'col_hinge.csv'],'ReadVariableNames',true);
hinge_table.id = []; % Omit id 

%% Go through each element and calculate the hinge properties
for i = 1:length(element_table.id)
    ele = element_table(i,:);
    
    % Calculate condition
    condition = 2; % UPDATE THIS TO READ FROM TABLE 10-11
    hinge = hinge_table(hinge_table.condition == condition,:);
    
    % Fitler Table based on P/Asfc
    p_ratio = ele.P/(ele.Ag*ele.fc_e);
    if p_ratio<=min(hinge.p_ratio)
        hinge = hinge(hinge.p_ratio == min(hinge.p_ratio),:);
    elseif p_ratio>=max(hinge.p_ratio)
        hinge = hinge(hinge.p_ratio == max(hinge.p_ratio),:);
    else
        hinge_1 = sortrows(hinge(hinge.p_ratio == max(hinge.p_ratio),:));
        hinge_2 = sortrows(hinge(hinge.p_ratio == min(hinge.p_ratio),:));
        hinge = hinge_1;
        hinge.a_hinge = ((hinge_1.a_hinge-hinge_2.a_hinge)/(max(hinge_1.p)-min(hinge_2.p)))*(p_ratio-min(hinge_2.p)) + hinge_2.a_hinge;
        hinge.b_hinge = ((hinge_1.b_hinge-hinge_2.b_hinge)/(max(hinge_1.p)-min(hinge_2.p)))*(p_ratio-min(hinge_2.p)) + hinge_2.b_hinge;
        hinge.c_hinge = ((hinge_1.c_hinge-hinge_2.c_hinge)/(max(hinge_1.p)-min(hinge_2.p)))*(p_ratio-min(hinge_2.p)) + hinge_2.c_hinge;
        hinge.io = ((hinge_1.io-hinge_2.io)/(max(hinge_1.p)-min(hinge_2.p)))*(p_ratio-min(hinge_2.p)) + hinge_2.io;
        hinge.ls = ((hinge_1.ls-hinge_2.ls)/(max(hinge_1.p)-min(hinge_2.p)))*(p_ratio-min(hinge_2.p)) + hinge_2.ls;
        hinge.cp = ((hinge_1.cp-hinge_2.cp)/(max(hinge_1.p)-min(hinge_2.p)))*(p_ratio-min(hinge_2.p)) + hinge_2.cp;
        hinge.p_ratio(:) = p_ratio;
    end
    
    % Filter table based on row
    row = ele.Av/(ele.b*ele.S);
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
        v_ratio = ele.V/(ele.b*ele.d*sqrt(ele.fc_e));
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
    
    % Plot Hinges
%     force_vector = [0,1,1.1,hinge.c_hinge];
%     disp_vector = [0
    
    % save as element hinge table
    ele_hinge(i,:) = hinge;
end

%% Save capacities to element database
element = [element_table ele_hinge];
writetable(element,['element.csv'])