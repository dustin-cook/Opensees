function [ ] = fn_make_directory( dir_path )
% Given a file path, check to see if a directory exists, if it does, wipe
% it and create an clean one.

if exist(dir_path,'dir')
    [status,msg] = rmdir(dir_path, 's'); % Remove old directory
    if status == 0
        disp(msg)
        error('Failed to delete analysis directory')
    end
end
mkdir(dir_path); % Make a fresh clear directory

end

