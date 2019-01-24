function [ element, element_TH, element_PM, joint ] = main_element_capacity( story, ele_prop_table, element, element_TH, analysis, joint )
% Description: Main script that calculates the strength of each element in 
% the model according to ASCE 41 

% Created By: Dustin Cook
% Date Created: 1/3/2019

% Inputs:

% Outputs:

% Assumptions:


%% Import Packages
import asce_41.*

%% Begin Method
for i = 1:length(element.id)
    ele = element(i,:);
    ele_id = ele.ele_id;
    ele_prop = ele_prop_table(ele_prop_table.id == ele_id,:);
    ele_TH = element_TH.(['ele_' num2str(element.id(i))]);
    
    %% Calculate Element Capacties
    [ ele, element_TH.(['ele_' num2str(element.id(i))]), element_PM.(['ele_' num2str(element.id(i))]) ] = fn_element_capacity( story, ele, ele_prop, ele_TH, analysis.nonlinear );
    disp([num2str(i), ' out of ', num2str(length(element.id)) ' elements complete' ])
    
    %% Caculate required development length and make sure there is enough
    [ ele.pass_aci_dev_length ] = fn_development_check( ele, ele_prop );
    
    % Save capcity for convergence check
    if strcmp(ele.type,'wall')
        ele.capacity = ele.Vn;
    else
        ele.capacity = ele.Mn_pos;
    end
    
    % Save data to main table
    ele_to_save(i,:) = ele;
end

element = ele_to_save;

%% Calculate Beam Column Strength 
for i =1:length(joint.id)
    jnt = joint(i,:);
    
    % Grab Beam and Column element properties
    beam_left = element(element.id == jnt.beam_left,:);
    if isempty(beam_left)
         beam_left_props = ele_prop_table(false,:);
    else
        beam_left_props = ele_prop_table(ele_prop_table.id == beam_left.ele_id,:);
    end
    beam_right = element(element.id == jnt.beam_right,:); 
    if isempty(beam_right)
         beam_right_props = ele_prop_table(false,:);
    else
        beam_right_props = ele_prop_table(ele_prop_table.id == beam_right.ele_id,:);
    end
    column_low = element(element.id == jnt.column_low,:); 
    if isempty(column_low)
         column_low_props = ele_prop_table(false,:);
    else
        column_low_props = ele_prop_table(ele_prop_table.id == column_low.ele_id,:);
    end
    column_high = element(element.id == jnt.column_high,:);
    if isempty(column_high)
         column_high_props = ele_prop_table(false,:);
    else
        column_high_props = ele_prop_table(ele_prop_table.id == column_high.ele_id,:);
    end
    
    % Specifiy joint properties based on surrounding element properties
    jnt.lambda = mean([column_low_props.lambda,column_high_props.lambda]); % Average of the two columns
    jnt.fc_e = mean([column_low_props.fc_e,column_high_props.fc_e]);
    jnt.e = mean([column_low_props.e,column_high_props.e]);
    jnt.iz = mean([column_low_props.iz,column_high_props.iz,beam_left_props.iz,beam_right_props.iz]);
    jnt.Pmax = max([column_high.Pmax,0]);
    
    % Calculate the joint area according to ASCE 41-17 10.4.2.3.2
    jnt.d = mean([column_low_props.d,column_high_props.d]); % Average of the two columns depths
    opt1 = mean([column_low_props.w,column_high_props.w]); % Average of the two columns widths
    opt2 = mean([beam_left_props.w,beam_right_props.w]) + jnt.d; % Average of the two beam widths
    opt3 = opt1; % Assumes beam axis frames to the center of the column axis
    jnt.w = min([opt1,opt2,opt3]);
    jnt.a = jnt.w*jnt.d;
    jnt.h = mean([beam_left_props.d,beam_right_props.d]) ; % Average of the two beam heights
    
    % Calculate whether the joint has conforming or non conforming reinforcement
    jnt.S = min([beam_left_props.S,beam_right_props.S,column_low_props.S,column_high_props.S]);
    h_c = mean([column_low_props.d,column_high_props.d]);
    if jnt.S <= h_c/2
        jnt.trans_rien = 'C';
    else
        jnt.trans_rien = 'NC';
    end
    
    % Calculate SCWB ratio
    beam_strength_1 = sum([beam_left.Mn_pos,beam_right.Mn_neg]); % case 1: 1 pos and 1 neg bending
    beam_strength_2 = sum([beam_left.Mn_neg,beam_right.Mn_pos]); % case 2: 1 pos and 1 neg bending the other way
    jnt.beam_strength = mean([beam_strength_1,beam_strength_2]); % Average of two cases
    col_strength_1 = sum([column_low.Mn_pos,column_high.Mn_pos]); % case 1: positive bending for both columns
    col_strength_2 = sum([column_low.Mn_neg,column_high.Mn_neg]); % case 2: negative bending for both columns
    jnt.column_strength = mean([col_strength_1,col_strength_2]); % Avearge of two cases
    jnt.col_bm_ratio = jnt.column_strength/jnt.beam_strength;
    
    % Caclulate the Joint Shear Strength according to ASCE 41-17 eq 10-4
    gamma_table = readtable(['+asce_41' filesep 'table_10_12_joint_gamma.csv'],'ReadVariableNames',true);
    jnt.gamma = gamma_table.gamma(strcmp(gamma_table.trans_rien,jnt.trans_rien) & strcmp(gamma_table.classification,jnt.class));
    jnt.Vj = jnt.lambda*jnt.gamma*sqrt(jnt.fc_e)*jnt.a;
    
    % Calculate Joint Moment Capacity
    jnt.Mn = jnt.Vj*jnt.h/2; % Comes from figure 9.3 of Moehle book. I think assumes beam yeild though, therefore should check and modify. I think also assumes there is a column above and column inflection points are at the center of the column height.
    
    % Calculate Joint Shear Demand
    Vcol = max([column_high.Vmax,0]); % Shear from the column above
    Ts1 = max([beam_left.Mmax / (beam_left_props.d*0.75),0]);
    C2 = max([beam_right.Mmax / (beam_right_props.d*0.75),0]); % Assuming jd is 75% of d, rough assumption for now until I take this further, also need to consider both directions
    jnt.Vmax = Ts1 + C2 - Vcol; % From Moehle Book EQ 9.2 Assumes the same as above.
    
    joint_to_save(i,:) = jnt;
end

joint = joint_to_save;

end

