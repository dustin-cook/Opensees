function [ joint ] = fn_joint_capacity( joint, element, ele_prop_table )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


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
    if sum(strcmp('Pmax',column_high.Properties.VariableNames)) == 0
        jnt.Pmax = 0;
    else
        jnt.Pmax = max([column_high.Pmax,0]);
    end
    
    % Calculate the joint area according to ASCE 41-17 10.4.2.3.2
    jnt.h = mean([column_low_props.h,column_high_props.h]); % Average of the two columns depths
    opt1 = mean([column_low_props.w,column_high_props.w]); % Average of the two columns widths
    opt2 = mean([beam_left_props.w,beam_right_props.w]) + jnt.h; % Average of the two beam widths
    opt3 = opt1; % Assumes beam axis frames to the center of the column axis
    jnt.w = opt1; % Use column props
    jnt.a = jnt.w*jnt.h;
    jnt.height = mean([beam_left_props.h,beam_right_props.h]) ; % Average of the two beam heights
    
    % Calculate whether the joint has conforming or non conforming reinforcement
    jnt.S = 999; % Set to non conforming for ISCB (update this to be database driven)
    h_c = mean([column_low_props.h,column_high_props.h]);
    if jnt.S <= h_c/2
        jnt.trans_rien = 'C';
    else
        jnt.trans_rien = 'NC';
    end
    
    % Calculate SCWB ratio
    beam_strength_1 = sum([beam_left.Mn_pos_2,beam_right.Mn_neg_1]); % case 1: 1 pos and 1 neg bending
    beam_strength_2 = sum([beam_left.Mn_neg_2,beam_right.Mn_pos_1]); % case 2: 1 pos and 1 neg bending the other way
    jnt.beam_strength = mean([beam_strength_1,beam_strength_2]); % Average of two cases
    col_strength_1 = sum([column_low.Mn_grav_pos_2,column_high.Mn_grav_pos_1]); % case 1: positive bending for both columns
    col_strength_2 = sum([column_low.Mn_grav_neg_2,column_high.Mn_grav_neg_1]); % case 2: negative bending for both columns
    jnt.column_strength = mean([col_strength_1,col_strength_2]); % Avearge of two cases
    jnt.col_bm_ratio = jnt.column_strength/jnt.beam_strength;

    % Determine implicit joint sfiffness condition according to ASCE 41-17 10.4.2.2
    if jnt.col_bm_ratio > 1.2
        jnt.implicit_stiff = 1;
    elseif jnt.col_bm_ratio < 0.8
        jnt.implicit_stiff = 2;
    else
        jnt.implicit_stiff = 3;
    end
    
    % Caclulate the Joint Shear Strength according to ASCE 41-17 eq 10-4
    gamma_table = readtable(['+asce_41' filesep 'table_10_12_joint_gamma.csv'],'ReadVariableNames',true);
    jnt.gamma = gamma_table.gamma(strcmp(gamma_table.trans_rien,jnt.trans_rien) & strcmp(gamma_table.classification,jnt.class));
    jnt.Vj = jnt.lambda*jnt.gamma*sqrt(jnt.fc_e)*jnt.a;
    
    % Calculate Joint Moment Capacity
    L = mean([beam_left.length,beam_right.length])/2; % Length from inflection point of beam to the column centerline, approximated as half neam centerline span,
    H = mean([column_low.length,column_high.length]); % Column Height measured bewteen column inflection points, approximated as the story height.
    jnt.Mn = jnt.Vj*(L/((L-jnt.w/2)/(0.9*jnt.height) - L/H)); % Equation 2.1 of Hassan and Moehle 2012 (with toa*A taken as Vy)
    
    % Calculate Joint Shear Demand
    if sum(strcmp('Vmax',column_high.Properties.VariableNames)) == 0 % Use capacities when demands do not yet exist
        Vcol = max([column_high.Vn_1,0]); % Shear from the column above
        Ts1 = max([beam_left.Mn_pos_2 / (beam_left_props.d_eff),0]);
        C2 = max([beam_right.Mn_pos_1 / (beam_right_props.d_eff),0]); % Assuming jd is 75% of d, rough assumption for now until I take this further, also need to consider both directions
        jnt.Vmax = Ts1 + C2 - Vcol; % From Moehle Book EQ 9.2 Assumes the same as above.
    else
        Vcol = max([column_high.Vmax,0]); % Shear from the column above
        Ts1 = max([beam_left.Mmax / (beam_left_props.d_eff),0]);
        C2 = max([beam_right.Mmax / (beam_right_props.d_eff),0]); % Assuming jd is 75% of d, rough assumption for now until I take this further, also need to consider both directions
        jnt.Vmax = Ts1 + C2 - Vcol; % From Moehle Book EQ 9.2 Assumes the same as above.
    end

    % Joint Drift
    if sum(strcmp('drift_x',column_low.Properties.VariableNames)) == 1
        jnt.drift_x = column_low.drift_x;
    end
    if sum(strcmp('drift_z',column_low.Properties.VariableNames)) == 1
        jnt.drift_z = column_low.drift_z;
    end
    
    % Check if joint is expected to yeild
    if jnt.Vmax > jnt.Vj
        jnt.yield = 1;
    else
        jnt.yield = 0;
    end
    
    % Save joint data to table
    joint_to_save(i,:) = jnt;
end

joint = joint_to_save;

end

