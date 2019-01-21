function [ ele ] = fn_disp_ductility( ele, ele_prop, story, eff_fyt_e )
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
    [ ele.vye ] = fn_vye( ele.type, ele.Mn_pos, ele.Mn_neg, ele.length, ele.gravity_load );
    [ ~, Vn, ~ ] = fn_aci_shear_capacity( ele_prop.fc_e, ele_prop.w, ele_prop.d, ele_prop.Av, eff_fyt_e, ele_prop.S, ele_prop.lambda, ele_prop.a, ele_prop.hw, ele.type, ele_prop.As_d, ele.P_grav );
    shear2use = min(ele.vye,Vn);
    yeild_disp = (shear2use*ele.length^3) / (12*ele_prop.e*ele_prop.iz);
    if sum(strcmp((['max_disp_' ele.direction{1}]),story.Properties.VariableNames)) == 1 && ele.story > 0
        try
        max_disp_demand = story.(['max_disp_' ele.direction{1}])(ele.story); % Should update this to a diplacement demand specific to the element itself instead of just the story
        catch
            test = 5;
        end
    else
        max_disp_demand = 0;
    end
    ele.disp_duct = max_disp_demand / yeild_disp;
%     ele_to_save(i,:) = ele;
% % end
% 
% element = ele_to_save;
end

