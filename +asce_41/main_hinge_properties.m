function [ element, joint ] = main_hinge_properties( ele_prop_table, element, joint )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import asce_41.*

%% Go through each element and calculate the hinge properties
for e_idx = 1:height(element)
    ele = element(e_idx,:);
    ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
    element.trans_rein_check{e_idx,1} = 'NA'; % Preallocate Transverse Reinforcement Compliance

    for ele_side = 1:2
        if strcmp(ele.type,'beam')
            [ hinge_props, element.trans_rein_check{e_idx,1} ] = fn_beam_hinge( ele, ele_props, ele_side );
            hinge_props_oop = hinge_props;
        elseif strcmp(ele.type,'column')
            [ hinge_props, element.rho_t(e_idx,1) ] = fn_col_hinge( ele, ele_props, ele_side );
            hinge_props_oop = hinge_props;
        elseif strcmp(ele.type,'wall')
            [ hinge_props ] = fn_wall_hinge( ele, ele_props, 0, ele_side );
            [ hinge_props_oop ] = fn_wall_hinge( ele, ele_props, 1, ele_side );
        end

        % save props into hinge table
        if strcmp(ele.type,'wall') && strcmp(ele.(['critical_mode_' num2str(ele_side)]),'shear')
            element.(['c_hinge_' num2str(ele_side)])(e_idx,1) = hinge_props.c_hinge;
            element.(['d_hinge_' num2str(ele_side)])(e_idx,1) = hinge_props.d_hinge;
            element.(['e_hinge_' num2str(ele_side)])(e_idx,1) = hinge_props.e_hinge;
            element.(['f_hinge_' num2str(ele_side)])(e_idx,1) = hinge_props.f_hinge;
            element.(['g_hinge_' num2str(ele_side)])(e_idx,1) = hinge_props.g_hinge;
        else
            element.(['a_hinge_' num2str(ele_side)])(e_idx,1) = hinge_props.a_hinge;
            element.(['b_hinge_' num2str(ele_side)])(e_idx,1) = hinge_props.b_hinge;
            element.(['c_hinge_' num2str(ele_side)])(e_idx,1) = hinge_props.c_hinge;
        end
        element.(['a_hinge_oop_' num2str(ele_side)])(e_idx,1) = hinge_props_oop.a_hinge;
        element.(['b_hinge_oop_' num2str(ele_side)])(e_idx,1) = hinge_props_oop.b_hinge;
        element.(['c_hinge_oop_' num2str(ele_side)])(e_idx,1) = hinge_props_oop.c_hinge;
        element.(['io_' num2str(ele_side)])(e_idx,1) = hinge_props.io;
        element.(['ls_' num2str(ele_side)])(e_idx,1) = hinge_props.ls;
        element.(['cp_' num2str(ele_side)])(e_idx,1) = hinge_props.cp;
    end
end

%% Go through each joint and calculate the nonlinear properties
for j_idx = 1:height(joint)
    jnt = joint(j_idx,:);
    
    [ hinge_props ] = fn_joint_hinge( jnt );
    
    % save as hinge parameters in joint table
    joint.a_hinge(j_idx,1) = hinge_props.a_hinge;
    joint.b_hinge(j_idx,1) = hinge_props.b_hinge;
    joint.c_hinge(j_idx,1) = hinge_props.c_hinge;
    joint.io(j_idx,1) = hinge_props.io;
    joint.ls(j_idx,1) = hinge_props.ls;
    joint.cp(j_idx,1) = hinge_props.cp;
end

end

