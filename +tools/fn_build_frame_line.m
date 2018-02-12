function [ ] = fn_build_frame_line(input_dir, col_id,beam_id,bay_length,direction,weight)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

id = 0;

for i = 1:length(direction)
    num_bays = bay_length.(direction{i});
    % Columns
    if col_id(i) > 0
        for j = 1:(length(num_bays)+1)
            id = id + 1;
            grid_line.id(id) = id;
            grid_line.grid_line_id(id) = i;
            grid_line.element_id(id) = col_id(i);
            if strcmp(direction{i},'x')
                grid_line.orientation(id) = 1;
            elseif strcmp(direction{i},'z')
                grid_line.orientation(id) = 4;
            else
                error('Grid Frame Direction Not Recognized')
            end
            if j == 1
                grid_line.x_start(id) = 0;
            else
                grid_line.x_start(id) = sum(num_bays(1:(j-1)));
            end
            grid_line.x_end(id) = grid_line.x_start(id);
            grid_line.y_start(id) = 0;
            grid_line.y_end(id) = 1;
            grid_line.weight(id) = 0;
        end
    end
    
    % Beams
    if beam_id(i) > 0
        for j = 1:length(num_bays)
            id = id + 1;
            grid_line.id(id) = id;
            grid_line.grid_line_id(id) = i;
            grid_line.element_id(id) = beam_id(i);
            if strcmp(direction{i},'x')
                grid_line.orientation(id) = 2;
            elseif strcmp(direction{i},'z')
                grid_line.orientation(id) = 3;
            else
                error('Grid Frame Direction Not Recognized')
            end
            if j == 1
                grid_line.x_start(id) = 0;
            else
                grid_line.x_start(id) = sum(num_bays(1:(j-1)));
            end
            grid_line.x_end(id) = grid_line.x_start(id) + num_bays(j);
            grid_line.y_start(id) = 1;
            grid_line.y_end(id) = 1;
            grid_line.weight(id) = weight;
        end
    end
   

end
T = table(grid_line.id', grid_line.grid_line_id', grid_line.element_id', grid_line.orientation', grid_line.x_start', grid_line.y_start', grid_line.x_end', grid_line.y_end', grid_line.weight');
T.Properties.VariableNames = {'id' 'grid_line_id' 'element_id' 'orientation' 'x_start' 'y_start' 'x_end' 'y_end' 'weight'};
writetable(T,[input_dir filesep 'grid_line.csv'],'WriteVariableNames',true);

end

