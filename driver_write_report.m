%% Script to Call Latex and Write Benchmarking Report
clear
close all
clc
rehash
rng('shuffle')

%% User Inputs
analysis.model_name = 'ICBS_model_3D_fixed';
analysis.proceedure = 'NDP'; % LDP or NDP or test
analysis.id = '38'; % ID of the analysis for it to create its own directory

%% Write PDF Report
% Define Directories
analysis_dir = ['outputs' '/' analysis.model_name '/' analysis.proceedure '_' analysis.id];

% Write LaTeX Variables File
fileID = fopen(['report_writers' filesep 'writeupVariables.tex'],'w');
fprintf(fileID,'\\newcommand{\\modelName}{%s}\n',strrep(analysis.model_name,'_',' '));
fprintf(fileID,'\\newcommand{\\analysisName}{%s}\n',strrep([analysis.proceedure '_' analysis.id],'_',' '));
fprintf(fileID,'\\newcommand{\\plotsDirDynamic}{%s}\n',[analysis_dir '/' 'analysis_plots']);
fprintf(fileID,'\\newcommand{\\plotsDirPushover}{%s}\n',[analysis_dir '/' 'pushover']);
fclose(fileID);
    
% Call LaTex Report Writer
command = ['pdflatex ' 'report_writers' filesep 'model_results_writeup.tex' ' -job-name=AnalysisReport' ' -output-directory=' analysis_dir ' -aux-directory=report_writers' filesep 'temp'];
[status,cmdout] = system(command,'-echo');