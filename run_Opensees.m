% Run Truss Tcl file in opensees
clear
close
clc

tic
%% Initial Setup 
import tools.*

%% Define Analysis parameters
analysis.id = 'test';
analysis.type = 4; % 1 = static force analysis, 2 = pushover analysis, 3 = dynamic analysis, 4 = calculate spectra
analysis.max_displ = 1; % only for pushover analsys

% Create Outputs Directory
output_dir = ['Analysis' filesep analysis.id];
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% Define EQ ground motion
analysis.eq_dir = 'ground_motions/ATC_63_far_feild';
eqs = dir([analysis.eq_dir filesep '*eq_*']);
num_eqs = length(eqs);

%% Define Model Parameters
% 2d linear model single frame, 1 bay, uniform members
num_bays = 0;
num_stories = 1;
story_ht_in = 1000;
bay_width_in = 240;
foundation_fix = [1 1 1];
story_force_k = 0;
story_weight_k = 0;
damp_ratio = 0.05;
story_mass = 1;
E = 60000;
A = 9999999999999;
I = 13824;

if analysis.type == 4
    periods = 0.05:0.05:3;
else
    periods = sqrt(story_mass/(3*E*I/story_ht_in^3));
end

%% Main Analysis
for i = 1:1
    % EQ this run
    analysis.eq_name = eqs(i).name;
    dt = 0.01;
    eq = load([analysis.eq_dir filesep analysis.eq_name]);
    
    for j = 1:length(periods)
        % Model properties for this run
        I = (story_mass * story_ht_in^3) / (3 * periods(j)^2 * E );
        period = sqrt(story_mass/(3*E*I/story_ht_in^3))
        
        %% Create Model Databases
        [ node, element ] = fn_model_table( num_stories, num_bays, story_ht_in, bay_width_in, foundation_fix, story_mass, story_weight_k, story_force_k, A, E, I );

        %% Write TCL file
        fn_build_model( analysis.id, node, element )
        fn_define_recorders( analysis, node.id, element.id )
        fn_define_loads( analysis, damp_ratio, node, dt )
        fn_define_analysis( analysis, node.id, dt, eq )

        %% Run Opensees
        command = ['opensees ' 'Analysis' filesep analysis.id filesep 'run_analysis.tcl'];
        system(command);

        %% Load outputs and plot
        % Nodal Displacement (i)
        node.disp_x = dlmread([output_dir filesep 'nodal_disp_x.txt'],' ')';
        node.disp_y = dlmread([output_dir filesep 'nodal_disp_y.txt'],' ')';

        % Nodal Reaction (k)
        node.reac_x = dlmread([output_dir filesep 'nodal_reaction_x.txt'],' ')';
        node.reac_y = dlmread([output_dir filesep 'nodal_reaction_y.txt'],' ')';
        
        % Nodal Acceleration (in/s2)
        node.accel_x = dlmread([output_dir filesep 'nodal_accel_x.txt'],' ')';
        node.accel_y = dlmread([output_dir filesep 'nodal_accel_y.txt'],' ')';

        % Element Forces
        for k = 1:length(element.id)
            ele_force = dlmread([output_dir filesep 'element_' num2str(k) '_force.txt'],' ');
            element.fx1{k} = ele_force(:,1)';
            element.fy1{k} = ele_force(:,2)';
            element.mz1{k} = ele_force(:,3)';
            element.fx2{k} = ele_force(:,4)';
            element.fy2{k} = ele_force(:,5)';
            element.mz2{k} = ele_force(:,6)';
        end

        %% Display Results
        if analysis.type == 1 % static load analysis
            disp(element)
        elseif analysis.type == 2 % pushover analysis
            if num_bays ==0
                plot(abs(node.disp_x(end,:)),abs(node.reac_x(1,:)));
            else
                plot(abs(node.disp_x(end,:)),abs(sum(node.reac_x(1:(num_bays+1),:))));
            end
            grid on
        elseif analysis.type == 3 % dynamic analysis
            plot(node.disp_x(end,:))
        elseif analysis.type == 4 % spectra analysis
            sa(i,j) = (max(abs(eq'+node.accel_x(end,:))))/386;
            sd(i,j) = max(abs(node.disp_x(end,:)));
            psa(i,j) = max(abs(node.disp_x(end,:)*((((2*pi)/periods(j))^2)/386)));
        else
            error('Unkown Analysis Type')
        end
    end
end

if analysis.type == 4 % spectra analysis
%     med_sa = median(sa);
%     med_sd = median(sd);
%     med_psa = median(psa);
    
    hold on
    grid on
    for i = 1:length(sa(:,1))
        plot(periods,sa(i,:),'b','lineWidth',0.5)
%         plot(periods,psa(i,:),'k','lineWidth',0.5)
    end
%     plot(periods,med_sa,'k','lineWidth',2)
%     plot(periods,med_sd,'r','lineWidth',2.5)
%     plot(periods,med_psa,'k','lineWidth',2)
end

toc
