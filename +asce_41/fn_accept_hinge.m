function [ hinge ] = fn_accept_hinge( element, ele_prop_table, hinge )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Import
import asce_41.fn_define_backbone_rot

%% Begin Method
for i = 1:length(hinge.id)
    ele_side = num2str(hinge.ele_side(i));
    ele = element(element.id == hinge.element_id(i),:);
    ele_prop = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
    if strcmp(hinge.direction{i},'oop')
        [ ~, ~, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', ele.(['Mn_oop_' ele_side]), ele.(['Mn_oop_' ele_side]), ele.(['Mp_oop_' ele_side]), ele.(['Mp_oop_' ele_side]), ele.length, ele_prop.e, ele_prop.iz, ele.(['a_hinge_oop_' ele_side]), ele.(['b_hinge_oop_' ele_side]), ele.(['c_hinge_oop_' ele_side]), 10, 0.1, ele.(['critical_mode_oop_' ele_side]) );
    elseif strcmp(hinge.direction{i},'primary')
        [ ~, ~, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', ele.(['Mn_pos_' ele_side]), ele.(['Mn_neg_' ele_side]), ele.(['Mp_pos_' ele_side]), ele.(['Mp_neg_' ele_side]), ele.length, ele_prop.e, ele_prop.iz, ele.(['a_hinge_' ele_side]), ele.(['b_hinge_' ele_side]), ele.(['c_hinge_' ele_side]), 10, 0.1, ele.(['critical_mode_' ele_side]) );
    end
    if abs(min(hinge.deformation_TH{i})) > abs(max(hinge.deformation_TH{i}))
        max_hinge_deform = max(abs(hinge.deformation_TH{i})) - rot_vec_neg(1); % Negative bending
    else
        max_hinge_deform = max(abs(hinge.deformation_TH{i})) - rot_vec_pos(1); % Positive Bending
    end
    if  max_hinge_deform <= ele.(['io_' ele_side])
        hinge.accept(i) = 1; % Passes IO
    elseif max_hinge_deform <= ele.(['ls_' ele_side])
        hinge.accept(i) = 2; % Passes LS
    elseif max_hinge_deform <= ele.(['cp_' ele_side])
        hinge.accept(i) = 3; % Passes CP
    else
        hinge.accept(i) = 4; % Fails all performance levels
    end
end

end

