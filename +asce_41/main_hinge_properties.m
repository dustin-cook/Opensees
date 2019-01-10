function [ element, joint ] = main_hinge_properties( ele_prop_table, element, joint )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import asce_41.*

%% Go through each element and calculate the hinge properties
for i = 1:height(element)
    ele = element(i,:);
    ele_props = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
    
    if strcmp(ele.type,'beam')
        [ hinge ] = fn_beam_hinge( ele, ele_props );
    elseif strcmp(ele.type,'column')
        [ hinge ] = fn_col_hinge( ele, ele_props );
    elseif strcmp(ele.type,'wall')
        [ hinge ] = fn_wall_hinge( ele, ele_props );
    end

    % save as element hinge table
    if strcmp(ele.type,'wall') && strcmp(ele.critical_mode,'shear')
        element.c_hinge(i,1) = hinge.c_hinge;
        element.d_hinge(i,1) = hinge.d_hinge;
        element.e_hinge(i,1) = hinge.e_hinge;
        element.f_hinge(i,1) = hinge.f_hinge;
        element.g_hinge(i,1) = hinge.g_hinge;
    else
        element.a_hinge(i,1) = hinge.a_hinge;
        element.b_hinge(i,1) = hinge.b_hinge;
        element.c_hinge(i,1) = hinge.c_hinge;
    end
    element.io(i) = hinge.io;
    element.ls(i) = hinge.ls;
    element.cp(i) = hinge.cp;
end

%% Go through each joint and calculate the nonlinear properties
for i = 1:height(joint)
    jnt = joint(i,:);
    
    [ hinge ] = fn_joint_hinge( jnt );
    
    % save as hinge parameters in joint table
    joint.a_hinge(i,1) = hinge.a_hinge;
    joint.b_hinge(i,1) = hinge.b_hinge;
    joint.c_hinge(i,1) = hinge.c_hinge;
    joint.io(i,1) = hinge.io;
    joint.ls(i,1) = hinge.ls;
    joint.cp(i,1) = hinge.cp;
end

end

