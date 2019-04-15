load([read_dir filesep 'element_TH_1.mat'])
hin_bot = load([read_dir filesep 'hinge_TH_1.mat']);
hin_top = load([read_dir filesep 'hinge_TH_3.mat']);

figure
hold on
plot(0.01:0.01:49.95,ele_TH.M_TH_1/1000,'LineWidth',1.5,'DisplayName','BeamColumn Element Force Recorder')
plot(0.01:0.01:50,hin_bot.hin_TH.force_TH/1000,'LineWidth',1.5,'DisplayName','Zero Length Element Force Recorder')
% plot(hin_top.hin_TH.force_TH)

grid on
grid minor
box on
legend('location','SouthEast')
xlabel('Time (s)')
ylabel('Moment Demand (k-in)')

figure 
plot(hin_bot.hin_TH.deformation_TH,hin_bot.hin_TH.force_TH,'LineWidth',1.5)