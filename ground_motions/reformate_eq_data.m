clear
close all
clc

g_conversion = 981; %(cm/s2)

eq_dir = ['ICSB_1979'];
filename = [eq_dir filesep 'eq_vert_ground.txt'];
eq_data = load(filename);
eq_data_g = eq_data/g_conversion;

fileID = fopen([eq_dir filesep 'eq_vert_ground.tcl'],'w');
for i = 1:length(eq_data_g(:,1))
    fprintf(fileID,'%d \n',eq_data_g(i,:));
end
fclose(fileID);