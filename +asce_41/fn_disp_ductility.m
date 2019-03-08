function [ disp_duct, vye, vye_oop ] = fn_disp_ductility( Mn_pos, Mn_neg, Mn_oop, ele, ele_prop, story, eff_fyt_e, Av, S )
%% Import Packages
import aci_318.fn_vye
import aci_318.fn_aci_shear_capacity

% Calculate ASCE 41 chapter 10 displacement ductility for each element
% for i = 1:length(element.id)
%     ele = element(i,:);
%     ele_id = ele.ele_id;
%     ele_prop = ele_prop_table(ele_prop_table.id == ele_id,:);
%     

%%% THIS ASSUMES THE COLUMNS IS FEXURALLY YIELDING ACCORDING TO C10.4.2.3.
%%% SHOULD CHECK THIS.
    [ vye ] = fn_vye( ele.type, Mn_pos, Mn_neg, ele.length, ele.gravity_load );
    [ vye_oop ] = fn_vye( ele.type, Mn_oop, Mn_oop, ele.length, ele.gravity_load );
    [ ~, Vn, ~ ] = fn_aci_shear_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.h, Av, eff_fyt_e, S, ele_prop.lambda, ele_prop.a, ele_prop.hw, ele.type, ele_prop.d_eff, ele.Pmax );
    shear2use = min(vye,Vn);
    yeild_disp = (shear2use*ele.length^3) / (12*ele_prop.e*ele_prop.iz);
    if sum(strcmp((['max_disp_' ele.direction{1}]),story.Properties.VariableNames)) == 1 && ele.story > 0
        max_disp_demand = story.(['max_disp_' ele.direction{1}])(ele.story); % Should update this to a diplacement demand specific to the element itself instead of just the story
    else
        max_disp_demand = 0;
    end
    disp_duct = max_disp_demand / yeild_disp;
%     ele_to_save(i,:) = ele;
% % end
% 
% element = ele_to_save;
end

