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

    if strcmp(ele.type,'beam')
        [ hinge_props, element.trans_rein_check{e_idx,1} ] = fn_beam_hinge( ele, ele_props );
        hinge_props_oop = hinge_props;
    elseif strcmp(ele.type,'column')
        [ hinge_props, element.rho_t(e_idx,1) ] = fn_col_hinge( ele, ele_props );
        hinge_props_oop = hinge_props;
    elseif strcmp(ele.type,'wall')
        [ hinge_props ] = fn_wall_hinge( ele, ele_props, 0 );
        [ hinge_props_oop ] = fn_wall_hinge( ele, ele_props, 1 );
    end

    % save props into hinge table
    if strcmp(ele.type,'wall') && strcmp(ele.critical_mode,'shear')
        element.c_hinge(e_idx,1) = hinge_props.c_hinge;
        element.d_hinge(e_idx,1) = hinge_props.d_hinge;
        element.e_hinge(e_idx,1) = hinge_props.e_hinge;
        element.f_hinge(e_idx,1) = hinge_props.f_hinge;
        element.g_hinge(e_idx,1) = hinge_props.g_hinge;
    else
        element.a_hinge(e_idx,1) = hinge_props.a_hinge;
        element.b_hinge(e_idx,1) = hinge_props.b_hinge;
        element.c_hinge(e_idx,1) = hinge_props.c_hinge;
    end
    element.a_hinge_oop(e_idx,1) = hinge_props_oop.a_hinge;
    element.b_hinge_oop(e_idx,1) = hinge_props_oop.b_hinge;
    element.c_hinge_oop(e_idx,1) = hinge_props_oop.c_hinge;
    element.io(e_idx,1) = hinge_props.io;
    element.ls(e_idx,1) = hinge_props.ls;
    element.cp(e_idx,1) = hinge_props.cp;
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

