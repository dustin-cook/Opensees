clear all
close all
clc

% Define Plot Directory
plt_dir = pwd;

% test plot schema
figure
[X,Y,Z] = peaks;
plt = surf(X,Y,Z);
title('Surf')
xlabel('test')
fn_ATC134_plot_schema( plt, plt_dir );

figure
x = 0:0.01:2*pi;
y = sin(x);
plt = plot(x,y);
title('Line')
xlabel('test')
fn_ATC134_plot_schema( plt, plt_dir );

figure
x = rand(100,1);
y = rand(100,1);
plt = scatter(x,y);
title('Scatter')
xlabel('test')
fn_ATC134_plot_schema( plt, plt_dir );
