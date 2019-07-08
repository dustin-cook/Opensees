function [ hinge ] = fn_accept_hinge( element, ele_prop_table, hinge, read_dir, node )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Import
import asce_41.fn_define_backbone_rot
import asce_41.fn_define_backbone_shear

%% Begin Method
for i = 1:height(hinge)
    if ~strcmp(hinge.type{i},'foundation')
        ele_side = num2str(hinge.ele_side(i));
        ele = element(element.id == hinge.element_id(i),:);
        ele_prop = ele_prop_table(ele_prop_table.id == ele.ele_id,:);
        node_x = node.x(node.id == ele.(['node_' ele_side]));
        load([read_dir filesep 'hinge_TH_' num2str(hinge.id(i)) '.mat'])
        if strcmp(ele.type,'wall') && strcmp(hinge.direction{i},'primary')
            [ ~, disp_vec ] = fn_define_backbone_shear( ele.(['Vn_' ele_side]), ele.length, ele_prop.g, ele_prop.av, ele.(['c_hinge_' ele_side]), ele.(['d_hinge_' ele_side]), ele.(['e_hinge_' ele_side]), ele.(['f_hinge_' ele_side]), ele.(['g_hinge_' ele_side])  );
        elseif strcmp(hinge.direction{i},'oop')
            [ ~, ~, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', ele.(['Mn_oop_' ele_side]), ele.(['Mn_oop_' ele_side]), ele.(['Mp_oop_' ele_side]), ele.(['Mp_oop_' ele_side]), ele.length, ele_prop.e, ele_prop.iz, ele.(['a_hinge_oop_' ele_side]), ele.(['b_hinge_oop_' ele_side]), ele.(['c_hinge_oop_' ele_side]), 10, 0.1, ele.(['critical_mode_oop_' ele_side]) );
            [ ~, ~, ele_rot_vec_pos, ele_rot_vec_neg ] = fn_define_backbone_rot( 'full', ele.(['Mn_oop_' ele_side]), ele.(['Mn_oop_' ele_side]), ele.(['Mp_oop_' ele_side]), ele.(['Mp_oop_' ele_side]), ele.length, ele_prop.e, ele_prop.iz, ele.(['a_hinge_oop_' ele_side]), ele.(['b_hinge_oop_' ele_side]), ele.(['c_hinge_oop_' ele_side]), 10, 0.1, ele.(['critical_mode_oop_' ele_side]) );
        elseif strcmp(hinge.direction{i},'primary')
            [ ~, ~, rot_vec_pos, rot_vec_neg ] = fn_define_backbone_rot( 'hinge', ele.(['Mn_pos_' ele_side]), ele.(['Mn_neg_' ele_side]), ele.(['Mp_pos_' ele_side]), ele.(['Mp_neg_' ele_side]), ele.length, ele_prop.e, ele_prop.iz, ele.(['a_hinge_' ele_side]), ele.(['b_hinge_' ele_side]), ele.(['c_hinge_' ele_side]), 10, 0.1, ele.(['critical_mode_' ele_side]) );
            [ ~, ~, ele_rot_vec_pos, ele_rot_vec_neg ] = fn_define_backbone_rot( 'full', ele.(['Mn_pos_' ele_side]), ele.(['Mn_neg_' ele_side]), ele.(['Mp_pos_' ele_side]), ele.(['Mp_neg_' ele_side]), ele.length, ele_prop.e, ele_prop.iz, ele.(['a_hinge_' ele_side]), ele.(['b_hinge_' ele_side]), ele.(['c_hinge_' ele_side]), 10, 0.1, ele.(['critical_mode_' ele_side]) );
        end

        % For rotation controlled elements
        if abs(min(hin_TH.deformation_TH)) > abs(max(hin_TH.deformation_TH))
            max_elastic_hinge_deform = rot_vec_neg(1); % Negative bending
            max_elastic_ele_deform = ele_rot_vec_neg(1);
        else
            max_elastic_hinge_deform = rot_vec_pos(1); % Positive bending
            max_elastic_ele_deform = ele_rot_vec_pos(1);
        end
        max_plastic_deform = max([max(abs(hin_TH.deformation_TH)) - max_elastic_hinge_deform,0]); 
        max_hinge_deform = max(abs(hin_TH.deformation_TH));
        if max_plastic_deform == 0
            max_ele_deform = max_hinge_deform*(max_elastic_ele_deform/max_elastic_hinge_deform); % Not Yeilding, deform is the max hinge deform times the ratio of elastic to hinge yeild point
        else
            max_ele_deform = max_plastic_deform + max_elastic_ele_deform; % If yielding take the plastic deform and add it to the element max elastic deform
        end

        % For Shear Walls
        max_ele_disp = max(abs(hin_TH.deformation_TH));


        % calculate which acceptance criteria it passes
        if  max_plastic_deform <= ele.(['io_' ele_side])
            hinge.accept(i) = 1; % Passes IO
        elseif max_plastic_deform <= ele.(['ls_' ele_side])
            hinge.accept(i) = 2; % Passes LS
        elseif max_plastic_deform <= ele.(['cp_' ele_side])
            hinge.accept(i) = 3; % Passes CP
        else
            hinge.accept(i) = 4; % Fails all performance levels
        end

        % calculate the ratio of the a, b, d, and e values and shear
        if strcmp(hinge.direction{i},'primary')
            if strcmp(ele.type,'wall')
                hinge.a_ratio(i) = NaN;
                hinge.b_ratio(i) = NaN;
                hinge.d_ratio(i) = max_ele_disp/disp_vec(3);
                hinge.e_ratio(i) = max_ele_disp/disp_vec(4);
            else
                hinge.a_ratio(i) = max_ele_deform/(ele.(['a_hinge_' ele_side]) + max_elastic_ele_deform);
                hinge.b_ratio(i) = max_ele_deform/(ele.(['b_hinge_' ele_side]) + max_elastic_ele_deform);
                hinge.d_ratio(i) = NaN;
                hinge.e_ratio(i) = NaN;
            end
            if strcmp(ele_side,'1') && node_x == 1571
                hinge.damage_recorded(i) = ele_prop.(['damage_' ele_side '_east']);
            else
                hinge.damage_recorded(i) = ele_prop.(['damage_' ele_side]);
            end
        elseif strcmp(hinge.direction{i},'oop')
            hinge.a_ratio(i) = max_ele_deform/(ele.(['a_hinge_oop_' ele_side]) + max_elastic_ele_deform);
            hinge.b_ratio(i) = max_ele_deform/(ele.(['b_hinge_oop_' ele_side]) + max_elastic_ele_deform);
            hinge.d_ratio(i) = NaN;
            hinge.e_ratio(i) = NaN;
            if strcmp(ele_side,'1') && node_x == 1571
                hinge.damage_recorded(i) = ele_prop.(['damage_' ele_side '_east']);
            else
                hinge.damage_recorded(i) = ele_prop.(['damage_oop_' ele_side]);
            end
        end
        hinge.V_ratio(i) = ele.(['Vmax_' ele_side])/ele.(['Vn_' ele_side]);
        hinge.P_ratio_expected(i) = ele.Pmax/(ele_prop.fc_e*ele_prop.a);
        hinge.P_ratio_nominal(i) = ele.Pmax/(ele_prop.fc_n*ele_prop.a);
        hinge.P_ratio_asce41(i) = ele.Pmax/ele.Pn_c;
    end
end

end

