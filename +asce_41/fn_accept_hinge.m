function [ hinge ] = fn_accept_hinge( element, hinge )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(hinge.id)
    ele_id = hinge.element_id(i);
    if max(hinge.deformation_TH{i}) <= element.io(element.id == ele_id)  
        hinge.accept(i) = 1; % Passes IO
    elseif max(hinge.deformation_TH{i}) <= element.ls(element.id == ele_id)
        hinge.accept(i) = 2; % Passes LS
    elseif max(hinge.deformation_TH{i}) <= element.cp(element.id == ele_id)
        hinge.accept(i) = 3; % Passes CP
    else
        hinge.accept(i) = 4; % Fails all performance levels
    end
end

end

