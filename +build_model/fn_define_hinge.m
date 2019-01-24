function [ node, element, hinge ] = fn_define_hinge( analysis, model, hinge, element, node, foundation_nodes_id )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% INITIAL SETUP
% Import Packages
import build_model.fn_create_hinge

%% Begin Method
hinge_id = 0;
if analysis.nonlinear ~= 0
    % Define Hinges
    for i = 1:length(element.id)
        if element.story(i) <= analysis.stories_nonlinear
            if strcmp(element.type{i},'column') || strcmp(element.type{i},'beam') % For all columns and beams
                % Define hinge at start of element
                hinge_id = hinge_id+1;
                [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_1', i, hinge_id, foundation_nodes_id, 'rotational', 'primary' ); 
                if strcmp(model.dimension,'3D') && strcmp(element.type{i},'column') % OOP hinges for 3D columns
                    hinge_id = hinge_id+1;
                    hinge.id(hinge_id,1) = hinge_id;
                    hinge.type{hinge_id,1} = hinge.type{hinge_id-1,1};
                    hinge.direction{hinge_id,1} = 'oop';
                    hinge.node_1(hinge_id,1) = hinge.node_1(hinge_id-1,1);
                    hinge.node_2(hinge_id,1) = hinge.node_2(hinge_id-1,1);
                    hinge.element_id(hinge_id,1) = hinge.element_id(hinge_id-1,1);
                end
                    
                % Define hinge at end of element
                hinge_id = hinge_id+1;
                [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_2', i, hinge_id, foundation_nodes_id, 'rotational', 'primary' );
                if strcmp(model.dimension,'3D') && strcmp(element.type{i},'column') % OOP hinges for 3D columns
                    hinge_id = hinge_id+1;
                    hinge.id(hinge_id,1) = hinge_id;
                    hinge.type{hinge_id,1} = hinge.type{hinge_id-1,1};
                    hinge.direction{hinge_id,1} = 'oop';
                    hinge.node_1(hinge_id,1) = hinge.node_1(hinge_id-1,1);
                    hinge.node_2(hinge_id,1) = hinge.node_2(hinge_id-1,1);
                    hinge.element_id(hinge_id,1) = hinge.element_id(hinge_id-1,1);
                end
            elseif strcmp(element.type{i},'wall') % For walls
                % Check if wall is the bottom wall
                walls_with_same_node = element.id(strcmp(element.type,'wall') & element.node_2 == element.node_1(i));
                if isempty(walls_with_same_node)
                    hinge_id = hinge_id+1;
                    % Define hinge at start of element
                    [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_1', i, hinge_id, foundation_nodes_id, 'shear', 'primary' ); 
                    if strcmp(model.dimension,'3D') % OOP hinges for 3D walls
                        hinge_id = hinge_id+1;
                        hinge.id(hinge_id,1) = hinge_id;
                        hinge.type{hinge_id,1} = 'rotational';
                        hinge.direction{hinge_id,1} = 'oop';
                        hinge.node_1(hinge_id,1) = hinge.node_1(hinge_id-1,1);
                        hinge.node_2(hinge_id,1) = hinge.node_2(hinge_id-1,1);
                        hinge.element_id(hinge_id,1) = hinge.element_id(hinge_id-1,1);
                    end
                end
            end
        end
    end
else
    % Define Shear Spings on Walls
    for i = 1:length(element.id)
        if strcmp(element.type{i},'wall')
            walls_with_same_node = element.id(strcmp(element.type,'wall') & element.node_2 == element.node_1(i));
            if isempty(walls_with_same_node)
                hinge_id = hinge_id+1;
                % Define hinge at start of element
                [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_1', i, hinge_id, foundation_nodes_id, 'shear', 'primary' ); 
            end
        end
    end
end

end

