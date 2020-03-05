clear all
close all
clc

% Define user inputs
model_name='tayo';
analysis_name='analysis_3';
element='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,201,202,203,301,302,303,401,402,403,501,502,503,601,602,603,701,702,703,801,802,803';
node='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32';
primary_nodes='8,12,16,20,24,28,32';
story_ht='3750,3000,3000,3000,3000,3000,3000';
period='1.7';

% Call function
status = driver_run_IDA(model_name, analysis_name, element, node, primary_nodes, story_ht, period);