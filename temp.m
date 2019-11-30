clear all
close all
clc

time =          [0 1 1 1.3 1.3 2.1 2.1 2.3 2.3  2.4  2.4 2.5 2.5  2.7  2.7  3    3 3.5];
functionality = [1 1 0 0   0.2 0.2 0.5 0.5 0.55 0.55 0.7 0.7 0.75 0.75 0.85 0.85 1 1];
plot(time,functionality,'k','linewidth',1.25);
hold on
plot([0,3.5],[1,1],'--k')
xlim([0.5,4])
ylim([0, 1.2])
ylabel('Building Functionality')
xlabel('Time')
box off