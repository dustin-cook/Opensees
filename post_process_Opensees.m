% Run Truss Tcl file in opensees
clear
close
clc

%% Load Analysis and Model parameters
analysis = readtable(['inputs' filesep 'analysis.csv'],'ReadVariableNames',true);
model_table = readtable(['inputs' filesep 'model.csv'],'ReadVariableNames',true);
model = model_table(model_table.id == analysis.model_id,:);
output_dir = ['outputs/' model.name{1}];

for i = 1:1%length(eqs)
    %% Load outputs and plot
    load([output_dir filesep 'analysis_data.mat'])
    
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
        savefig([output_dir filesep 'Roof_Displacemnet.fig'])
        figure
        hold on
        grid on
        xlabel('time (s)')
        ylabel('Roof Acceleration (g)')
        plot((1:length(eq))*analysis.eq_dt,node.accel_x_abs(end,:)/386,'DisplayName','Roof')
        plot((1:length(eq))*analysis.eq_dt,node.accel_x_abs(1,:)/386,'DisplayName','Ground')
        legend('show')
        savefig([output_dir filesep 'Roof_Acceleration.fig'])
        hold off
        max_displ_all_nodes = max(abs(node.disp_x),[],2);
        max_accel_all_nodes = max(abs(node.accel_x_rel),[],2);
        max_displ_profile = max(max_displ_all_nodes(1:(model.num_bays+1)));
        max_accel_profile = max(max_accel_all_nodes(1:(model.num_bays+1)));
        for story = 1:(model.num_stories-1)
            max_displ_profile = [max_displ_profile, max(max_displ_all_nodes(story*num_node_story+(1:(model.num_bays+1))))];
            max_accel_profile = [max_accel_profile, max(max_accel_all_nodes(story*num_node_story+(1:(model.num_bays+1))))];
        end
        max_displ_profile = [max_displ_profile, max(max_displ_all_nodes((end-(model.num_bays*3+1)):end))];
        max_accel_profile = [max_accel_profile, max(max_accel_all_nodes((end-(model.num_bays*3+1)):end))];
        max_drift_profile = (max_displ_profile(2:end) - max_displ_profile(1:end-1))./story_ht;
    else
        error('Unkown Analysis Type')
    end
    
    % Save Data
    save([output_dir filesep 'analysis_data'])

end

