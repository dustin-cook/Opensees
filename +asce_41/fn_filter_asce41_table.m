function [ table ] = fn_filter_asce41_table( table, value, col_name, output_names )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if value<=min(table.(col_name)) % Less than the given column min
    table = table(table.(col_name) == min(table.(col_name)),:); % filter based on the min criteria 
elseif value>=max(table.(col_name)) % more than the given column max
    table = table(table.(col_name) == max(table.(col_name)),:); % filter based on the max criteria
else % Linear interpolation between the values
    min_values = sortrows(table(table.(col_name) == max(table.(col_name)),:));
    max_values = sortrows(table(table.(col_name) == min(table.(col_name)),:));
    table = min_values;
    for i = 1:length(output_names)
        table.(output_names{i}) = ((min_values.(output_names{i})-max_values.(output_names{i}))/(max(min_values.(col_name))-min(max_values.(col_name))))*(value-min(max_values.(col_name))) + max_values.(output_names{i});
    end
    table.(col_name)(:) = value;
end

end

