function [ hinge ] = fn_accept_hinge( hinge, output_dir )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Load in element table with hinge properties (UPDATE THIS TO RUN HINGE CALCS ON NEW LOADS)
element = readtable([output_dir filesep 'element_linear.csv'],'ReadVariableNames',true);

for i = 1:length(hinge.id)
    ele_id = hinge.element_id(i);
    if hinge.rotation(i) <= element.io(element.id == ele_id)  
        hinge.accept(i) = 1; % Passes IO
    elseif hinge.rotation(i) <= element.ls(element.id == ele_id)
        hinge.accept(i) = 2; % Passes LS
    elseif hinge.rotation(i) <= element.cp(element.id == ele_id)
        hinge.accept(i) = 3; % Passes CP
    else
        hinge.accept(i) = 4; % Fails all performance levels
    end
end

end

