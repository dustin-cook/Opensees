function [ ele ] = fn_disp_ductility( ele, ele_prop, story )
%% Import Packages
import aci_318.fn_vye

% Calculate ASCE 41 chapter 10 displacement ductility for each element
% for i = 1:length(element.id)
%     ele = element(i,:);
%     ele_id = ele.ele_id;
%     ele_prop = ele_prop_table(ele_prop_table.id == ele_id,:);
%     
    [ ele.vye ] = fn_vye( ele.type, ele.Mn_pos, ele.Mn_neg, ele.length, ele.gravity_load );
    yeild_disp = (ele.vye*ele.length^3) / (12*ele_prop.e*ele_prop.iz);
    if isfeild(story,['max_disp_' ele.direction{1}])
        max_disp_demand = story.(['max_disp_' ele.direction{1}])(ele.story); % Should update this to a diplacement demand specific to the element itself instead of just the story
    else
        max_disp_demand = 0;
    end
    ele.disp_duct = max_disp_demand / yeild_disp;
%     ele_to_save(i,:) = ele;
% % end
% 
% element = ele_to_save;
end

