function [ ] = fn_plot_spectra(node, type, gm, plot_dir, read_dir_opensees, pile_model)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Notes
% Roof Nodes: 6170, 6171, 6176, 6177

%% Initial Setup
import plotting_tools.*

% Defin color scheme
matlab_colors =  [0, 0.4470, 0.7410;
                  0.8500, 0.3250, 0.0980;
                  0.9290, 0.6940, 0.1250;
                  0.4940, 0.1840, 0.5560;
                  0.4660, 0.6740, 0.1880;
                  0.3010, 0.7450, 0.9330;
                  0.6350, 0.0780, 0.1840];
% Define Recorder channels
if strcmp(type,'Roof')
    channel.x = 'chan_4';
    channel.z = 'chan_2';
elseif strcmp(type,'Ground')
    if pile_model
        channel.x = 'free_field_1';
        channel.z = 'free_field_3';
    else
        channel.x = 'chan_13';
        channel.z = 'chan_11';
        channel.y = 'chan_12';
    end
end

%% Run Spectra Tool to generate roof response spectra
if strcmp(type,'Roof')
    % Set Inputs
    spectra_write_dir = '..\toolbelt\Spectra Tool\inputs\ICSB_analysis';
    center_x = 671;
    center_z = 300;
    node_roof_center_id = node.id(node.x == center_x & node.z == center_z & node.story == 6);
    roof_center_id_TH = load([read_dir_opensees filesep 'node_TH_' num2str(node_roof_center_id) '.mat']);
    fileID = fopen([spectra_write_dir filesep 'roof_accel_x.txt'],'w');
    fprintf(fileID,'%d\n',gm.x.eq_dt);
    fprintf(fileID,'%d\n',roof_center_id_TH.nd_TH.accel_x_abs_TH');
    fclose(fileID);
    fileID = fopen([spectra_write_dir filesep 'roof_accel_z.txt'],'w');
    fprintf(fileID,'%d\n',gm.x.eq_dt);
    fprintf(fileID,'%d\n',roof_center_id_TH.nd_TH.accel_z_abs_TH');
    fclose(fileID);

    % Run Method
    gm_set_name = 'ICSB_analysis';
    g_factor = 1; % denominator to convert raw data to g
    dt_idx = 1; % location of the dt value
    eq_data_idx = 2; % location of the start of the ground motion
    plot_spectra = 1;
    run('..\toolbelt\Spectra Tool\main_calc_spectra(gm_set_name, g_factor, dt_idx, eq_data_idx, plot_spectra)')
end

%% Load Data
spectra_read_dir = '..\toolbelt\Spectra Tool\outputs\ICSB_analysis';
if strcmp(type,'Roof')
    spectra_analysis_x = readtable([spectra_read_dir filesep 'roof_accel_x' filesep 'spectra.csv'],'ReadVariableNames', true);
    spectra_analysis_z = readtable([spectra_read_dir filesep 'roof_accel_z' filesep 'spectra.csv'],'ReadVariableNames', true);
end
spectra_record_x = readtable(['ground_motions' filesep 'ICSB_recordings' filesep channel.x '_accel' filesep 'spectra.csv'],'ReadVariableNames', true);
spectra_record_z = readtable(['ground_motions' filesep 'ICSB_recordings' filesep channel.z '_accel' filesep 'spectra.csv'],'ReadVariableNames', true);
if isfield(gm,'y') && isfield(channel,'y')
    spectra_record_y = readtable(['ground_motions' filesep 'ICSB_recordings' filesep channel.y '_accel' filesep 'spectra.csv'],'ReadVariableNames', true);
end

%% Plot Spectra
hold on
if strcmp(type,'Roof')
    plot(spectra_analysis_z.period,spectra_analysis_z.psa_5,'color',matlab_colors(1,:),'LineWidth',1,'DisplayName','NS Analysis')
    plot(spectra_record_z.period,spectra_record_z.psa_5,'--','color',matlab_colors(1,:),'LineWidth',1,'DisplayName','NS Recording')
else
    plot(spectra_record_z.period,spectra_record_z.psa_5,'color',matlab_colors(1,:),'LineWidth',1,'DisplayName','Plan NS (Y)')
end

if strcmp(type,'Roof')
    plot(spectra_analysis_x.period,spectra_analysis_x.psa_5,'color',matlab_colors(2,:),'LineWidth',1,'DisplayName','EW Analysis')
    plot(spectra_record_x.period,spectra_record_x.psa_5,'--','color',matlab_colors(2,:),'LineWidth',1,'DisplayName','EW Recording')
else
    plot(spectra_record_x.period,spectra_record_x.psa_5,'color',matlab_colors(2,:),'LineWidth',1,'DisplayName','Plan EW (X)')
end

if isfield(gm,'y') && isfield(channel,'y')
    plot(spectra_record_y.period,spectra_record_y.psa_5,'color',matlab_colors(3,:),'LineWidth',1,'DisplayName','Plan Vertical (Y)')
end
        
xlabel('Period, T_1(s)')
ylabel('Spectral Acceleration, Sa (\xi = 5%) (g)')
box on 
legend('Location','northeast')
legend boxoff
set(gca,'FontSize',14)
xlim([0,3])
fn_format_and_save_plot( [plot_dir filesep 'Spectra Plots'], [type ' Spectra'], 0 )
hold off
close

end

