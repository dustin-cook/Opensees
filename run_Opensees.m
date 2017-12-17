% Run Truss Tcl file in opensees
clear
close
clc

tic
%% Initial Setup 
import tools.*

%% Load Analysis and Model parameters
analysis = readtable(['inputs' filesep 'analysis.csv'],'ReadVariableNames',true);
model_table = readtable(['inputs' filesep 'model_2.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);

%% Start Analysis
story_force_k = [0 0 0 0 0 0];

% Create Outputs Directory
output_dir = ['outputs/' analysis.model_name{1}];
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% Define EQ ground motion
eqs = dir([analysis.eq_dir{1} filesep '*eq_*']);
num_eqs = length(eqs);

% Break up model table to be used by this method (should be improved later)
[ I_bm, A_bm, E_bm, I_col, A_col, E_col, damp_ratio, story_weight_k, foundation_fix, story_ht_in, bay_width_in, num_stories, num_bays, floor_ht, story_mass, building_ht, building_mass ] = fn_model_translation( model );

% choose periods to run based on analysis type
if analysis.type == 4
    periods = 0.001:0.05:3.001;
else
    periods = sqrt(building_mass/(3*E_col(1)*I_col(1)/building_ht^3))*2*pi();
end

% Main Opensees Analysis
for i = 1:1%length(eqs)
    % EQ this run
    eq = load([analysis.eq_dir{1} filesep analysis.eq_name{1}]);
    
    for j = 1:length(periods)
        % Model properties for this run
%         I_col = (building_mass * building_ht^3) / (3 * (periods(j)/(2*pi))^2 * E_col(1) );
%         T = sqrt(building_mass/(3*E_col(1)*I_col/building_ht^3))*2*pi;
        
        %% Create Model Databases
        [ node, element ] = fn_model_table( num_stories, num_bays, floor_ht, bay_width_in, foundation_fix, story_mass, story_weight_k, story_force_k, A_col, E_col, I_col, A_bm, E_bm, I_bm );

        %% Write TCL file
        fn_build_model( output_dir, node, element )
        fn_define_recorders( output_dir, analysis.type, node.id, element.id )
        fn_define_loads( output_dir, analysis, damp_ratio, node, analysis.eq_dt )
        fn_define_analysis( output_dir, analysis, node.id, eq, analysis.eq_dt )
        fn_eigen_analysis( output_dir, analysis.time_step, node.id )
        
        %% Run Opensees
        command = ['opensees ' output_dir filesep 'run_analysis.tcl'];
        system(command);

        %% Load outputs and plot
        % Nodal Displacement (i)
        node.disp_x = dlmread([output_dir filesep 'nodal_disp_x.txt'],' ')';
        node.disp_y = dlmread([output_dir filesep 'nodal_disp_y.txt'],' ')';

        % Nodal Reaction (k)
        node.reac_x = dlmread([output_dir filesep 'nodal_reaction_x.txt'],' ')';
        node.reac_y = dlmread([output_dir filesep 'nodal_reaction_y.txt'],' ')';
        
        % Nodal Acceleration (in/s2)
        node.accel_x_rel = dlmread([output_dir filesep 'nodal_accel_x.txt'],' ')';
        node.accel_y_rel = dlmread([output_dir filesep 'nodal_accel_y.txt'],' ')';
        node.accel_x_abs = node.accel_x_rel + ones(length(node.id),1)*eq'*386;
        node.accel_y_abs = node.accel_y_rel + ones(length(node.id),1)*eq'*386;

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
            figure
            plot((1:length(eq))*analysis.eq_dt,node.disp_x(end,:))
            xlabel('time (s)')
            ylabel('Roof Displacemnet (in)')
            grid on
            figure
            hold on
            grid on
            xlabel('time (s)')
            ylabel('Roof Acceleration (g)')
            plot((1:length(eq))*analysis.eq_dt,node.accel_x_abs(end,:)/386,'DisplayName','Roof')
            plot((1:length(eq))*analysis.eq_dt,node.accel_x_abs(1,:)/386,'DisplayName','Ground')
            legend('show')
            hold off
            max_displ_all_nodes = max(abs(node.disp_x),[],2);
            for story = 1:num_stories+1
                max_displ_profile(story) = max(max_displ_all_nodes((story-1)*(num_bays+1)+1:story*(num_bays+1)));
            end
            max_drift_profile = (max_displ_profile(2:end) - max_displ_profile(1:end-1))./story_ht_in;
            max_accel_all_nodes = max(abs(node.accel_x_rel),[],2);
            for story = 1:num_stories+1
                max_accel_profile(story) = max(max_accel_all_nodes((story-1)*(num_bays+1)+1:story*(num_bays+1)));
            end
        elseif analysis.type == 4 % spectra analysis
            total_eq_time = length(eq)*analysis.eq_dt;
            old_time_vec = analysis.eq_dt:analysis.eq_dt:total_eq_time;
            new_time_vec = analysis.time_step:analysis.time_step:total_eq_time;
            new_eq_vec = interp1([0,old_time_vec],[0;eq],new_time_vec);
            
            sa(i,j) = max(abs(new_eq_vec+node.accel_x_rel(end,:)/386));
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
