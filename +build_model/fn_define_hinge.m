function [ node, element, hinge ] = fn_define_hinge( analysis, model, hinge, element, node, foundation_nodes_filter )
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
        if element.story(i) <= analysis.stories_nonlinear && ~element.rigid(i) && ~element.elastic(i)
            if strcmp(element.type{i},'column') || (strcmp(element.type{i},'beam') && ~analysis.elastic_beams) % For all columns and beams
                % Define hinge at start of element
                hinge_id = hinge_id+1;
                [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_1', i, hinge_id, foundation_nodes_filter, 'rotational', 'primary', 1 ); 
                if strcmp(model.dimension,'3D') && strcmp(element.type{i},'column') % OOP hinges for 3D columns
                    hinge_id = hinge_id+1;
                    hinge.id(hinge_id,1) = hinge_id;
                    hinge.type{hinge_id,1} = hinge.type{hinge_id-1,1};
                    hinge.direction{hinge_id,1} = 'oop';
                    hinge.ele_side(hinge_id,1) = 1;
                    hinge.node_1(hinge_id,1) = hinge.node_1(hinge_id-1,1);
                    hinge.node_2(hinge_id,1) = hinge.node_2(hinge_id-1,1);
                    hinge.element_id(hinge_id,1) = hinge.element_id(hinge_id-1,1);
                    hinge.story(hinge_id,1) = element.story(i);
                    hinge.ele_direction{hinge_id,1} = element.direction{i};
                end
                    
                % Define hinge at end of element
                hinge_id = hinge_id+1;
                [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_2', i, hinge_id, foundation_nodes_filter, 'rotational', 'primary', 2 );
                if strcmp(model.dimension,'3D') && strcmp(element.type{i},'column') % OOP hinges for 3D columns
                    hinge_id = hinge_id+1;
                    hinge.id(hinge_id,1) = hinge_id;
                    hinge.type{hinge_id,1} = hinge.type{hinge_id-1,1};
                    hinge.direction{hinge_id,1} = 'oop';
                    hinge.ele_side(hinge_id,1) = 2;
                    hinge.node_1(hinge_id,1) = hinge.node_1(hinge_id-1,1);
                    hinge.node_2(hinge_id,1) = hinge.node_2(hinge_id-1,1);
                    hinge.element_id(hinge_id,1) = hinge.element_id(hinge_id-1,1);
                    hinge.story(hinge_id,1) = element.story(i);
                    hinge.ele_direction{hinge_id,1} = element.direction{i};
                end
            elseif strcmp(element.type{i},'wall') && ~analysis.fiber_walls % Walls: Only for lumped plasticity models
                % Check if wall is the bottom wall
                walls_with_same_node = element.id(strcmp(element.type,'wall') & element.node_2 == element.node_1(i));
                if isempty(walls_with_same_node)
                    hinge_id = hinge_id+1;
                    % Define hinge at start of element
                    [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_1', i, hinge_id, foundation_nodes_filter, 'shear', 'primary', 1 ); 
                    if strcmp(model.dimension,'3D') % OOP hinges for 3D walls
                        hinge_id = hinge_id+1;
                        hinge.id(hinge_id,1) = hinge_id;
                        hinge.type{hinge_id,1} = 'rotational';
                        hinge.direction{hinge_id,1} = 'oop';
                        hinge.ele_side(hinge_id,1) = 1;
                        hinge.node_1(hinge_id,1) = hinge.node_1(hinge_id-1,1);
                        hinge.node_2(hinge_id,1) = hinge.node_2(hinge_id-1,1);
                        hinge.element_id(hinge_id,1) = hinge.element_id(hinge_id-1,1);
                        hinge.story(hinge_id,1) = element.story(i);
                        hinge.ele_direction{hinge_id,1} = element.direction{i};
                    end
                end
            end
        end
    end
else
    % Define Linear Shear Spings on Walls to account for shear deflections
    for i = 1:length(element.id)
        if strcmp(element.type{i},'wall') && ~analysis.fiber_walls % Only for lumped plasticity models
            walls_with_same_node = element.id(strcmp(element.type,'wall') & element.node_2 == element.node_1(i));
            if isempty(walls_with_same_node)
                hinge_id = hinge_id+1;
                % Define hinge at start of element
                [ node, element, hinge ] = fn_create_hinge( node, element, hinge, 'node_1', i, hinge_id, foundation_nodes_filter, 'shear', 'primary', 1 ); 
            end
        end
    end
end

end

