function [ node, element, hinge ] = fn_define_hinge( hinge, element, node, foundation_nodes_id, calc_capacity, ele_props )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% INITIAL SETUP
% Import Packages
import build_model.*
import aci_318.fn_aci_moment_capacity


%% Define Hinges
hinge_id = 0;
for i = 1:length(element.id)
    if strcmp(element.type{i},'column') || strcmp(element.type{i},'beam') % For all columns and beams
        hinge_id = hinge_id+1;
        % Define hinge at start of element
        [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_1', i, hinge_id, foundation_nodes_id ); 
        hinge_id = hinge_id+1;
        % Define hinge at end of element
        [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_2', i, hinge_id, foundation_nodes_id );
    end

    if calc_capacity
        % Calculate the intial (zero axial load) Moment Capacity of Each element in the model
        ele = ele_props(ele_props.id == element.ele_id(i),:);
        if strcmp(ele.type,'truss') || contains(ele.description,'rigid')
            element.Mn_pos(i,1) = 9999999999;
            element.Mn_neg(i,1) = 9999999999;
        else
            % Moment Capcity per ACI
            [ ~, element.Mn_pos(i,1) ] = fn_aci_moment_capacity( 'pos', ele.fc_e, ele.w, ele.d, ele.As, ele.As_d, ele.fy_e, ele.Es, 0 ); % change to be based on the gravity load instead?
            [ ~, element.Mn_neg(i,1) ] = fn_aci_moment_capacity( 'neg', ele.fc_e, ele.w, ele.d, ele.As, ele.As_d, ele.fy_e, ele.Es, 0 );
        end
    end
end

end

