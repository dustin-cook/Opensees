clear
close
clc

%% Inpits
direction = 'x';
channel = '4';

%% Load Data
spectra_read_dir = 'C:\Users\DustinCook\Desktop\Repos\toolbelt\Spectra Tool\outputs\ICSB_analysis';
spectra_analysis_1 = readtable([spectra_read_dir filesep 'roof_accel_' direction filesep 'spectra.csv'],'ReadVariableNames', true);
% spectra_analysis_2 = readtable(['ICBS_model_3D_fixed' filesep 'roof_accel_x_2' filesep 'spectra.csv'],'ReadVariableNames', true);
% spectra_analysis_3 = readtable(['ICBS_model_3D_fixed' filesep 'roof_accel_x_3' filesep 'spectra.csv'],'ReadVariableNames', true);
% spectra_analysis_4 = readtable(['ICBS_model_3D_fixed' filesep 'roof_accel_x_4' filesep 'spectra.csv'],'ReadVariableNames', true);
spectra_record = readtable(['ICSB_recordings' filesep 'chan_' channel '_accel' filesep 'spectra.csv'],'ReadVariableNames', true);

%% Plot Spectra
hold on
plot(spectra_analysis_1.period,spectra_analysis_1.psa_3,'b','LineWidth',1.5,'DisplayName','Analysis')
% plot(spectra_analysis_2.period,spectra_analysis_2.psa_3,'b','LineWidth',1.5,'DisplayName','Analysis')
% plot(spectra_analysis_3.period,spectra_analysis_3.psa_3,'b','LineWidth',1.5,'DisplayName','Analysis')
% plot(spectra_analysis_4.period,spectra_analysis_4.psa_3,'b','LineWidth',1.5,'DisplayName','Analysis')
plot(spectra_record.period,spectra_record.psa_3,'r','LineWidth',1.5,'DisplayName','Recording')
xlabel('Period (s)')
ylabel('Sa of roof at 3% damping (g)')
grid on
box on 
legend('Location','northeast')
set(gca,'FontSize',15)
xlim([0,3])
savefig(['Roof Spectra ' direction '.fig'])
saveas(gcf,['Roof Spectra ' direction '.png'])
hold off
close

% dir1 = fileread('C:\Users\Dustin\Desktop\Repos\Opensees\ground_motions\ICSB_recordings\chan_13_accel\chan_13_accel.txt');
% data_raw_1 = str2double(strsplit(dir1,' '));
% data_1 = data_raw_1(~isnan(data_raw_1));
% plot(linspace(0,5810*0.01,5810),[0,data_1(3:end)],'LineWidth',1.5)
% pga_dir1 = max(abs(data_1(3:end)))
% 
% figure
% dir1_alt = fileread('C:\Users\Dustin\Desktop\Repos\Opensees\ground_motions\ICSB_1979\gm_ew.tcl');
% data_raw_1_alt = str2double(strsplit(dir1_alt,' '));
% data_1_alt = data_raw_1_alt(~isnan(data_raw_1_alt));
% pga_dir1_FF = max(abs(data_1_alt))
% hold on
% plot(linspace(0,5810*0.01,5810),[0,data_1(3:end)],'LineWidth',1.5)
% plot(linspace(0,5689*0.01,5689),[0,data_1_alt],'LineWidth',1.5)
% hold off
% 
% figure
% dir2 = fileread('C:\Users\Dustin\Desktop\Repos\Opensees\ground_motions\ICSB_recordings\chan_11_accel\chan_11_accel.txt');
% data_raw_2 = str2double(strsplit(dir2,' '));
% data_2 = data_raw_2(~isnan(data_raw_2));
% plot(linspace(0,5809*0.01,5809),[0,data_2(3:end)],'LineWidth',1.5)
% pga_dir1 = max(abs(data_2(3:end)))
% 
% figure
% dir3 = fileread('C:\Users\Dustin\Desktop\Repos\Opensees\ground_motions\ICSB_recordings\chan_12_accel\chan_12_accel.txt');
% data_raw_3 = str2double(strsplit(dir3,' '));
% data_3 = data_raw_3(~isnan(data_raw_3));
% plot(linspace(0,5809*0.01,5809),[0,data_3(3:end)],'LineWidth',1.5)
% pga_dir1 = max(abs(data_3(3:end)))