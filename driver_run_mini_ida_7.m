%% Script to Run and IDA of a building with a single ground motion %%
clear all
close all
clc

%% Initial Setup
% Set Input Parameters
scale_factors = [0.7];
name = 'scale_7';

% Sa Values
sa_x_t1 = 0.35;
sa_z_t1 = 0.88;

% Load in Model Tables
node = readtable('mini_ida/model_tables/node.csv','readVariableNames',true);
story = readtable('mini_ida/model_tables/story.csv','readVariableNames',true);

% Import packages
import opensees.post_process.*
import plotting_tools.*

% make folders
mkdir('mini_ida/model/outputs')
mkdir('mini_ida/ida_results')
mkdir(['mini_ida/ida_results' filesep name])

%% Run Opensees Model
for i = 1:length(scale_factors)
    % Change EQ Scale factor
    scale_factor = scale_factors(i);
    sa_x(i,1) = sa_x_t1*scale_factor;
    sa_z(i,1) = sa_z_t1*scale_factor;
    
    % Write master run file
    file_name = ['mini_ida/model/master_run_analysis_' name '.tcl'];
    fileID = fopen(file_name,'w');
    fprintf(fileID,'wipe\n');
    fprintf(fileID,'source mini_ida/model/model.tcl\n'); 
    fprintf(fileID,'source mini_ida/model/loads.tcl\n');
    fprintf(fileID,'source mini_ida/model/recorders_%s.tcl\n',name);
    fprintf(fileID,'source mini_ida/model/run_analysis.tcl\n');
    fclose(fileID);
    
    % Write Recorders File
    file_name = ['mini_ida/model/recorders_' name '.tcl'];
    fileID = fopen(file_name,'w');
    fprintf(fileID,'puts "Defining Recorders ..."\n');
    for n = 1:height(node)
       if node.record_disp(n)
            fprintf(fileID,'recorder Node -xml mini_ida/model/outputs/%s_nodal_disp_%s.xml -time -node %i -dof 1 3 disp\n',name,num2str(node.id(n)),node.id(n));
       end
    end
    fprintf(fileID,'puts "Defining Recorders Complete"\n');
    fclose(fileID);
    
    % Write GM file
    file_name = ['mini_ida/model/set_gm_and_scale.tcl'];
    fileID = fopen(file_name,'w');
    fprintf(fileID,'puts "Defining Ground Motions and Scale ..." \n');
    fprintf(fileID,'timeSeries Path 1 -dt $dt -filePath ground_motions/ICSB_1979/chan_13_accel_15.tcl -factor %f \n',386);
    fprintf(fileID,'pattern UniformExcitation 3 1 -accel 1 -fact %f \n',scale_factor);
    fprintf(fileID,'timeSeries Path 2 -dt $dt -filePath ground_motions/ICSB_1979/chan_11_accel_15.tcl -factor %f \n',386);
    fprintf(fileID,'pattern UniformExcitation 4 3 -accel 2 -fact %f \n',scale_factor);
    fclose(fileID);

    % Call Opensees
    command = ['/projects/duco1061/software/OpenSeesSP/bin/OpenSeesSP ' 'mini_ida' filesep 'model' filesep 'master_run_analysis_' name '.tcl'];
%     command = ['openseesSP ' 'mini_ida' filesep 'model' filesep 'master_run_analysis_' name '.tcl'];
    [status,cmdout] = system(command,'-echo');
    
    % Grab displacements
    for n = 1:height(node)
           if node.record_disp(n)
               [ node_disp_raw ] = fn_xml_read(['mini_ida' filesep 'model' filesep 'outputs' filesep name '_nodal_disp_' num2str(node.id(n)) '.xml']);
               node_disp_raw = node_disp_raw'; % flip to be node per row
               node.disp_x{n} = node_disp_raw(2,:); 
               node.max_disp_x(n) = max(abs(node_disp_raw(2,:)));
               node.disp_z{n} = node_disp_raw(3,:);  
               node.max_disp_z(n) = max(abs(node_disp_raw(3,:)));
           else
               node.disp_x{n} = [];
               node.max_disp_x(n) = NaN;
               node.disp_z{n} = [];
               node.max_disp_z(n) = NaN;
           end
    end
    
    % Calc Story Drift
    [ story.max_drift_x ] = fn_drift_profile( node.disp_x, story, node );
    [ story.max_drift_z ] = fn_drift_profile( node.disp_z, story, node );
    
    max_drift_x(i,1) = max(story.max_drift_x);
    max_drift_z(i,1) = max(story.max_drift_z);
end

%% Plot and Save Results
plot(max_drift_x,sa_x)
xlabel('Max Drift')
ylabel('Sa(T_1) (g)')
fn_format_and_save_plot( ['mini_ida' filesep 'ida_results' filesep name], 'IDA Plot EW Frame Direction', 2 )

plot(max_drift_z,sa_z)
xlabel('Max Drift')
ylabel('Sa(T_1) (g)')
fn_format_and_save_plot( ['mini_ida' filesep 'ida_results' filesep name], 'IDA Plot NS Wall Direction', 2 )

% Save results as csv
ida.sa_x = sa_x;
ida.sa_z = sa_z;
ida.max_drift_x = max_drift_x;
ida.max_drift_z = max_drift_z;
ida_table = struct2table(ida);
writetable(ida_table,['mini_ida' filesep 'ida_results' filesep name filesep 'ida_results_.csv'])
