function [ ] = fn_build_frame_line(input_dir, col_id, beam_id, wall_id, bay_length, direction, grid_lines, x_offset_start, x_offset_end, trib_wt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

id = 0;

for i = 1:length(direction)
    num_bays = bay_length.(direction{i})(grid_lines{i}.bays);
    % Columns --- COULD MAKE FUNCTION
    if col_id(i) > 0
        for j = 1:(length(num_bays)+1)
            id = id + 1;
            grid_line.id(id,1) = id;
            grid_line.grid_line_id(id,1) = i;
            grid_line.element_id(id,1) = col_id(i);
            if strcmp(direction{i},'x')
                grid_line.orientation(id,1) = 1;
            elseif strcmp(direction{i},'z')
                grid_line.orientation(id,1) = 4;
            else
                error('Grid Frame Direction Not Recognized')
            end
            if j == 1
                grid_line.x_start(id,1) = 0;
            else
                grid_line.x_start(id,1) = sum(num_bays(1:(j-1)));
            end
            grid_line.x_end(id,1) = grid_line.x_start(id);
            grid_line.y_start(id,1) = 0;
            grid_line.y_end(id,1) = 1;
            grid_line.trib_wt(id,1) = 0;
        end
    end
    
    % Beams --- COULD MAKE FUNCTION
    if beam_id(i) > 0
        for j = 1:length(num_bays)
            id = id + 1;
            grid_line.id(id,1) = id;
            grid_line.grid_line_id(id,1) = i;
            grid_line.element_id(id,1) = beam_id(i);
            if strcmp(direction{i},'x')
                grid_line.orientation(id,1) = 2;
            elseif strcmp(direction{i},'z')
                grid_line.orientation(id,1) = 3;
            else
                error('Grid Frame Direction Not Recognized')
            end
            if j == 1
                grid_line.x_start(id,1) = 0;
            else
                grid_line.x_start(id,1) = sum(num_bays(1:(j-1)));
            end
            grid_line.x_end(id,1) = grid_line.x_start(id) + num_bays(j);
            grid_line.y_start(id,1) = 1;
            grid_line.y_end(id,1) = 1;
            grid_line.trib_wt(id,1) = trib_wt(i);
        end
    end
   
    % Walls --- COULD MAKE FUNCTION
    if wall_id(i) > 0
        for j = 1:length(num_bays)
            id = id + 1;
            grid_line.id(id,1) = id;
            grid_line.grid_line_id(id,1) = i;
            grid_line.element_id(id,1) = wall_id(i);
            if strcmp(direction{i},'x')
                grid_line.orientation(id,1) = 1;
            elseif strcmp(direction{i},'z')
                grid_line.orientation(id,1) = 4;
            else
                error('Grid Frame Direction Not Recognized')
            end
            if j == 1
                grid_line.x_start(id,1) = 0 + x_offset_start(i);
            else
                grid_line.x_start(id,1) = sum(num_bays(1:(j-1))) + x_offset_start(i);
            end
            grid_line.x_end(id,1) = grid_line.x_start(id) - x_offset_start(i) + num_bays(j) - x_offset_end(i);
            grid_line.y_start(id,1) = 0;
            grid_line.y_end(id,1) = 1;
            grid_line.trib_wt(id,1) = 0;
        end
    end

end
% T = table(grid_line.id', grid_line.grid_line_id', grid_line.element_id', grid_line.orientation', grid_line.x_start', grid_line.y_start', grid_line.x_end', grid_line.y_end');
% T.Properties.VariableNames = {'id' 'grid_line_id' 'element_id' 'orientation' 'x_start' 'y_start' 'x_end' 'y_end'};
T = struct2table(grid_line);
writetable(T,[input_dir filesep 'grid_line.csv'],'WriteVariableNames',true);

end

