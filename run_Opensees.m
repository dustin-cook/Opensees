% Run Truss Tcl file in opensees
clear
close
clc

tic
%% Initial Setup 
import tools.*

%% Define Analysis parameters
analysis.id = 'test';
analysis.type = 3; % 1 = static force analysis, 2 = pushover analysis, 3 = dynamic analysis, 4 = calculate spectra
analysis.max_displ = 1; % only for pushover analsys
analysis.time_step = 0.01;

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
num_bays = 5;
num_stories = 6;
story_ht_in = [16.5 13.5 13.5 13.5 13.5 13.2]*12;
bay_width_in = 300;
foundation_fix = [1 1 1];
story_weight_k = [644.53125	681.1523438	681.1523438	681.1523438	681.1523438	556.640625];
story_mass = story_weight_k/386;
story_force_k = [0 0 0 0 0 0];
damp_ratio = 0.02;
E_col = 63000;
A_col = 576;
I_col = 27648;
E_bm = 71000;
A_bm = 720;
I_bm = 54000;

%% Start Analysis
% Calculate Additional Building Parameters
floor_ht = 0;
for i = 1:length(story_ht_in)
    floor_ht = [floor_ht sum(story_ht_in(1:i))];
end
building_ht = sum(story_ht_in);
building_mass = sum(story_mass);

% choose periods to run based on analysis type
if analysis.type == 4
    periods = 0.001:0.05:3.001;
else
    periods = sqrt(building_mass/(3*E_col*I_col/building_ht^3))*2*pi();
end

% Main Opensees Analysis
for i = 1:1%length(eqs)
    % EQ this run
    analysis.eq_name = eqs(i).name;
    dt = 0.01;
    eq = load([analysis.eq_dir filesep analysis.eq_name]);
    
    for j = 1:length(periods)
        % Model properties for this run
        I_col = (building_mass * building_ht^3) / (3 * (periods(j)/(2*pi))^2 * E_col );
        T = sqrt(building_mass/(3*E_col*I_col/building_ht^3))*2*pi;
        
        %% Create Model Databases
        [ node, element ] = fn_model_table( num_stories, num_bays, floor_ht, bay_width_in, foundation_fix, story_mass, story_weight_k, story_force_k, A_col, E_col, I_col, A_bm, E_bm, I_bm );

        %% Write TCL file
        fn_build_model( analysis.id, node, element )
        fn_define_recorders( analysis, node.id, element.id )
        fn_define_loads( analysis, damp_ratio, node, dt )
        fn_define_analysis( analysis, node.id, eq, dt )
        fn_eigen_analysis( analysis, node.id )
        
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
            max_displ_all_nodes = max(abs(node.disp_x),[],2);
            for story = 1:num_stories+1
                max_displ_profile(story) = max(max_displ_all_nodes((story-1)*(num_bays+1)+1:story*(num_bays+1)));
            end
            max_drift_profile = (max_displ_profile(2:end) - max_displ_profile(1:end-1))./story_ht_in;
            max_accel_all_nodes = max(abs(node.accel_x),[],2);
            for story = 1:num_stories+1
                max_accel_profile(story) = max(max_accel_all_nodes((story-1)*(num_bays+1)+1:story*(num_bays+1)));
            end
        elseif analysis.type == 4 % spectra analysis
            total_eq_time = length(eq)*dt;
            old_time_vec = dt:dt:total_eq_time;
            new_time_vec = analysis.time_step:analysis.time_step:total_eq_time;
            new_eq_vec = interp1([0,old_time_vec],[0;eq],new_time_vec);
            
            sa(i,j) = max(abs(new_eq_vec+node.accel_x(end,:)/386));
            sd(i,j) = max(abs(node.disp_x(end,:)));
            psa(i,j) = ((2*pi/T)^2)*sd(i,j)/386;
        else
            error('Unkown Analysis Type')
        end
    end
end

%% Post Process Plotting
if analysis.type == 4 % spectra analysis
    med_sa = median(sa);
    med_sd = median(sd);
    med_psa = median(psa);
    
    figure
    hold on
    grid on
    xlabel('Period (s)')
    ylabel('Sa (g)')
    for i = 1:length(sa(:,1))
        plot3 = plot(periods,sa(i,:),'b','lineWidth',0.5);
        plot3 = plot(periods,psa(i,:),'k','lineWidth',0.5);
    end
%     plot(periods,med_sa,'r','lineWidth',2)
%     plot(periods,med_psa,'k','lineWidth',2)
    saveas(plot3,'DC_spectra.png')
    hold off
    close
end

toc
