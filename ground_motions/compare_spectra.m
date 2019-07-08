clear
close
clc

%% Inputs
plot_name = 'NS Free field v base motion';
spectra_1_path = 'ICSB_1979\spectra_chan_11_accel.csv';
spectra_2_path = 'ICSB_1979\spectra_gm_ns_free_field.csv';
spectra_1_name = 'NS Base Recording';
spectra_2_name = 'NS Free Field';

%% Load Data
spectra_1 = readtable(spectra_1_path,'ReadVariableNames', true);
spectra_2 = readtable(spectra_2_path,'ReadVariableNames', true);

%% Plot Spectra
hold on
plot(spectra_1.period,spectra_1.psa_5,'b','LineWidth',1.5,'DisplayName',spectra_1_name)
plot(spectra_2.period,spectra_2.psa_5,'r','LineWidth',1.5,'DisplayName',spectra_2_name)
xlabel('Period (s)')
ylabel('Spectral Acceleration (\xi = 5%) (g)')
grid on
box on 
legend('Location','northeast')
set(gca,'FontSize',15)
xlim([0,3])
savefig([plot_name '.fig'])
saveas(gcf,[plot_name '.png'])
hold off
close
