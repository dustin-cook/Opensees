function [ story ] = fn_accel_filter(node, story, filter_freq_range, read_dir, write_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Initial Setup
import file_exchange.fn_balaji_accel_filter
import opensees.post_process.*

% Load ground motion data
load([read_dir filesep 'gm_data.mat'])

%% Filter accels
low_freq = 0; % hardcode to no high pass filter
for n = 1:height(node)
    if node.record_disp(n)
        load([read_dir filesep 'node_TH_' num2str(node.id(n)) '.mat'],'nd_TH')
        for d = 1:2
%                 [ node_accel_filtered ] = fn_fft_accel_filter( nd_TH.accel_x_rel_TH, 0.01, 0.01:0.01:50, analysis.filter_high_freq, low_freq );            
            eq_length = ground_motion.(dirs_ran{d}).eq_length;
            eq_dt = ground_motion.(dirs_ran{d}).eq_dt;
            eq_timespace = linspace(eq_dt,eq_length*eq_dt,eq_length);
            [ node_accel_filtered_butter ] = fn_balaji_accel_filter( [eq_timespace',nd_TH.accel_x_rel_TH'], eq_dt, filter_freq_range, 5, 'bandpass');
            node_accel_filtered = node_accel_filtered_butter(:,2)';
            
            nd_TH.(['accel_' dirs_ran{d} '_rel_TH']) = node_accel_filtered; % Convert to G  
            nd_TH.(['accel_' dirs_ran{d} '_abs_TH']) = node_accel_filtered + eq.(dirs_ran{d})';
            node.(['max_accel_' dirs_ran{d} '_rel'])(n) = max(abs(node_accel_filtered));
            node.(['max_accel_' dirs_ran{d} '_abs'])(n) = max(abs(node_accel_filtered + eq.(dirs_ran{d})'));


        end
        save([write_dir filesep 'node_TH_' num2str(node.id(n)) '.mat'],'nd_TH')
    end
end

for d = 1:2
    [ story.(['max_accel_' dirs_ran{d}]) ] = fn_calc_max_repsonse_profile( node.(['max_accel_' dirs_ran{d} '_abs']), story, node, 0 );
    story.(['max_accel_center_' dirs_ran{d}]) = node.(['max_accel_' dirs_ran{d} '_abs'])(node.center == 1 & node.record_accel == 1 & node.story > 0);
end


end

