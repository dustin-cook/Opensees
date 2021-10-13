%% Script to Run and IDA of a building with a single ground motion %%
clear all
close all
clc

%% Assumptions
% 1) 3D model

%% User Inputs
% Define Model
% analysis.model_id = 4;
analysis.model_type = 3; % 1 = SDOF, 2 = MDOF (default), 3 = Archetype model
analysis.proceedure = 'P58'; 
analysis.id = 'Dissertation_Study'; % ID of the analysis for it to create its own directory
analysis.gm_set = 'FEMA_far_field';
analysis.run_z_motion = 0;

% Analysis options
analysis.summit = 0;
analysis.run_parallel = 0;
analysis.run_ida = 0;
analysis.post_process_ida = 0;
analysis.create_fragilities = 1;
% analysis.plot_ida = 0;
% analysis.detialed_post_process = 0;
analysis.run_sa_stripes = 1;
% analysis.scale_increment = 0.25;
% analysis.sa_stripes = [0.2 0.4];
analysis.collapse_drift = 0.10; % Drift ratio at which we are calling the model essentially collapsed
analysis.clear_existing_data = 0;
analysis.general_ida = 0;
analysis.write_xml = 1;

% Nonlinear options
analysis.nonlinear = 1;
analysis.stories_nonlinear = inf;
analysis.stories_nonlinear_low = 0;
analysis.elastic_beams = 0;
analysis.fiber_walls = 0;
analysis.hinge_stiff_mod = 10;

% Secondary options
analysis.dead_load = 1; % Dead load factor
analysis.live_load = 1; % live load factor
analysis.opensees_SP = 0;
analysis.type = 1;
analysis.damping = 'rayleigh';
analysis.damp_ratio = 0.03;
analysis.solution_algorithm = 1;
analysis.initial_timestep_factor = 1;
analysis.suppress_outputs = 1;
analysis.play_movie = 0;
analysis.movie_scale = 10;
analysis.algorithm = 'KrylovNewton';
analysis.integrator = 'Newmark 0.5 0.25';
analysis.joint_model = 1;
analysis.joint_explicit = 0;
analysis.simple_recorders = 1;
analysis.additional_elements = 1;

% Site inputs
% Curt's thesis site
% site.lat = 33.996;
% site.lng = -118.162;
% Generic LA site
site.lat = 34.05;
site.lng = -118.25;
site.vs30 = 537;
site.sds = 1.0;
site.sd1 = 0.6;
hazard.afe = 1/475;

% Define models to run
model_data = readtable(['inputs' filesep 'archetype_models.csv'],'ReadVariableNames',true);
model_data = model_data(model_data.num_stories == 12,:);
% model_data = model_data(model_data.num_stories ~= 20,:);
% model_data = model_data(model_data.ie ~= 1.25,:);
% model_data = model_data(~contains(model_data.name,'drift15'),:);
num_models = height(model_data);

% Define remote directory
remote_dir = ['G:\My Drive\Dissertation Archetype Study\Archetypes RCMF\Archetype Model Responses'];

%% Import Packages
import ida.*
import asce_7.*
import build_model.main_build_model
import opensees.write_tcl.fn_define_model
import opensees.main_eigen_analysis
import usgs.*

%% Pull Hazard Data
rps2run = [43, 72, 108, 224, 475, 975, 2475, 4975];
% rps2run = [10, 50, 100, 150, 250, 500, 750, 1000, 1500, 2500];
[sa_spectra, sa_periods] = fn_call_USGS_hazard_API('E2014', site.lat, site.lng, site.vs30, 1./rps2run);

col_data = table;
p_exceed_ac = [];
for m = 1:num_models % run for each model     
    %% Initial Setup
    % Load basic model data
    model = model_data(m,:);
    analysis.model_id = model.id;
    fprintf('Running Model %i of %i: %s\n', m, num_models, model.name{1})

    % Define read and write directories
    main_dir = ['outputs' '/' model.name{1} '/' analysis.proceedure '_' analysis.id];
    model_dir = [main_dir '/' 'model_data'];
    if ~exist(model_dir,'dir')
        mkdir(model_dir)
    end
    tcl_dir = [main_dir '/' 'opensees_data'];
    if ~exist(tcl_dir,'dir')
        mkdir(tcl_dir)
    end
    model_remote_dir = [remote_dir filesep model.name{1}];
    if ~exist(model_remote_dir,'dir')
        mkdir(model_remote_dir)
    end
%     ELFP_model_dir = ['outputs' '/' model.name{1} '/' 'ELFP' '_' analysis.id '/' 'model_data'];
    % tcl_dir = ['outputs' '/' model.name{1} '/' 'ELFP' '_' analysis.id '/' 'eigen_analysis'];
    % asce41_dir = [main_dir '/' 'asce_41_data'];

    % Load ground motion data
    gm_set_table = readtable(['ground_motions' filesep analysis.gm_set filesep 'ground_motion_set.csv'],'ReadVariableNames',true);
    gm_median_pga = median(gm_set_table.pga);

    %% Build Model
    analysis.out_dir = main_dir;
    disp('Building Model ...')

    % Create model tables
    main_build_model( model, analysis, [] )

    % Load Model Data
    node = readtable([model_dir filesep 'node.csv'],'readVariableNames',true);
    story = readtable([model_dir filesep 'story.csv'],'readVariableNames',true);
    hinge = readtable([model_dir filesep 'hinge.csv'],'readVariableNames',true);
    element = readtable([model_dir filesep 'element.csv'],'readVariableNames',true);
    joint = readtable([model_dir filesep 'joint.csv'],'readVariableNames',true);

    % Write model tcl file
    [ ~ ] = fn_define_model( tcl_dir, node, element, joint, hinge, analysis, model.dimension, story, [], model );

    %% Run Eigen Analysis
%     analysis.nonlinear = 0;
    [ model ] = main_eigen_analysis( model, analysis );

    %% Modify Model Data
    % Set period variable
    model.T1_x
    ida_results.period = model.T1_x;

    % Factor Loads
    [ element ] = fn_factor_loads( analysis, element, [] );

    %% Interpolate Site Hazard for building period
%     [design_values] = fn_call_USGS_design_value_API(1, 'asce7-16', site.lat, site.lng, 'II', model.site_class);
    period_cropped = min(max(model.T1_x,min(sa_periods)),max(sa_periods)); % limit building period to the available periods
%     sa_475 = interp1(sa_periods,sa_spectra,period_cropped);
    sa_levels = interp1(sa_periods,sa_spectra,period_cropped);
    
    % design level
    sa_dbe = min(site.sds,site.sd1/model.T1_x);

    % MCE level
    ida_results.mce = sa_dbe*1.5;
    
%     analysis.sa_stripes = [sa_475, sa_dbe];
    analysis.sa_stripes = sa_levels;

    %% Run Opensees Models
    if analysis.run_ida || analysis.post_process_ida
        fn_master_IDA(analysis, model, story, element, node, hinge, joint, gm_set_table, ida_results, tcl_dir, main_dir)
    end
    
    %% Create Response and Consequence Fragilities
    if analysis.create_fragilities
        write_dir = [main_dir '/' 'IDA' '/' 'Fragility Data'];
        if ~exist(write_dir,'dir')
            mkdir(write_dir)
        end
        [col_med, col_beta, p_col_mce, p_col_dbe] = fn_create_collapse_fragility(analysis, gm_set_table, ida_results, main_dir, write_dir);
        col_data.id(m,1) =  model.id;
        col_data.model_name(m,1) =  model.name;
        col_data.med_sa(m,1) = col_med;
        col_data.beta(m,1) = col_beta;
        col_data.p_col_mce(m,1) = p_col_mce;
        col_data.p_col_dbe(m,1) = p_col_dbe;
    end

    %% Post Process EDP data
    disp('Collecting EDP data for FEMA P-58 Assessment ...')
    id = 0;
    idr = [];
    pfa = [];
    max_a_ratio = [];
    for gm = 1:height(gm_set_table)
        gm_dir = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_table.set_id(gm)) '_' num2str(gm_set_table.pair(gm))];
        sa_folders = dir([gm_dir filesep 'Sa_*']);
        for s = 1:length(sa_folders)
            % Load data
            outputs_dir = [main_dir '/' 'IDA' '/' 'Summary Data' '/' 'GM_' num2str(gm_set_table.set_id(gm)) '_' num2str(gm_set_table.pair(gm)) '/' sa_folders(s).name];
            outputs_file = [outputs_dir filesep 'summary_results.mat'];
            story_file = [outputs_dir filesep 'story_analysis.mat'];
            hinge_file = [outputs_dir filesep 'hinge_analysis.mat'];
            
            % Set values
            max_a_ratio.gm(gm,1) = gm;
            sa_val = round(str2double(strrep(strrep(sa_folders(s).name,'Sa_',''),'_','.')),3);
            
            % Load Data
            if exist(outputs_file,'file')
                load(outputs_file)
            else
                error('NEED TO EXPLICITLY HANDLE MISSING SUMMARY TABLE')
            end
            
            % Component response
            if exist(hinge_file,'file')
                load(hinge_file)
                if summary.collapse > 0 % set collapse values = NaN
                    max_a_ratio.(sa_folders(s).name)(gm,1) = NaN;
                else
                    max_a_ratio.(sa_folders(s).name)(gm,1) = max(hinge.a_ratio);
                end
            else
                max_a_ratio.(sa_folders(s).name)(gm,1) = NaN; % didnt collect hinge properties for this ground motion
            end
            
            % EDPs
            if summary.collapse == 0
                if exist(story_file,'file')
                    load(story_file)

                    id = id + 1;
                    % Formulate the IDR table for the P-58 assessment - X direction
                    idr.sa(id,1) = sa_val;
                    idr.direction(id,1) = 1;
                    idr.gm(id,1) = gm;%gm_set_table.eq_name{gm};
                    for n = 1:height(story)
                        if summary.collapse > 0 % set collapse values = NaN
                            idr.(['story_' num2str(n)])(id,1) = NaN;
                        else
                            idr.(['story_' num2str(n)])(id,1) = story.max_drift_x(n);
                        end
                    end

                    % Formulate the PFA table for the P-58 assessment - X direction
                    pfa.sa(id,1) = sa_val;
                    pfa.direction(id,1) = 1;
                    pfa.gm(id,1) = gm;%gm_set_table.eq_name{gm};
                    pfa.floor_1(id,1) = summary.pga_x;
                    for n = 1:height(story)
                        if summary.collapse > 0 % set collapse values = NaN
                            pfa.(['floor_' num2str(n+1)])(id,1) = NaN;
                        else
                            pfa.(['floor_' num2str(n+1)])(id,1) = story.max_accel_x(n);
                        end
                    end

                    id = id + 1;

                    % Formulate the IDR table for the P-58 assessment - Z direction 
                    if analysis.run_z_motion
                        idr.sa(id,1) = sa_val;
                        idr.direction(id,1) = 2;
                        idr.gm(id,1) = gm;%gm_set_table.eq_name{gm};
                        for n = 1:height(story)
                            if summary.collapse > 0 % set collapse values = NaN
                                idr.(['story_' num2str(n)])(id,1) = NaN;
                            else
                                idr.(['story_' num2str(n)])(id,1) = story.max_drift_z(n);
                            end
                        end
                    else
                        idr.sa(id,1) = sa_val;
                        idr.direction(id,1) = 2;
                        idr.gm(id,1) = gm;%gm_set_table.eq_name{gm};
                        for n = 1:height(story)
                            if summary.collapse > 0 % set collapse values = NaN
                                idr.(['story_' num2str(n)])(id,1) = NaN;
                            else
                                idr.(['story_' num2str(n)])(id,1) = story.max_drift_x(n);
                            end
                        end
                    end

                    % Formulate the PFA table for the P-58 assessment - Z direction 
                    if analysis.run_z_motion
                        pfa.sa(id,1) = sa_val;
                        pfa.direction(id,1) = 2;
                        pfa.gm(id,1) = gm;%gm_set_table.eq_name{gm};
                        pfa.floor_1(id,1) = summary.pga_z;
                        for n = 1:height(story)
                            if summary.collapse > 0 % set collapse values = NaN
                                pfa.(['floor_' num2str(n+1)])(id,1) = NaN;
                            else
                                pfa.(['floor_' num2str(n+1)])(id,1) = story.max_accel_z(n);
                            end
                        end
                    else
                        pfa.sa(id,1) = sa_val;
                        pfa.direction(id,1) = 2;
                        pfa.gm(id,1) = gm;%gm_set_table.eq_name{gm};
                        pfa.floor_1(id,1) = summary.pga_x;
                        for n = 1:height(story)
                            if summary.collapse > 0 % set collapse values = NaN
                                pfa.(['floor_' num2str(n+1)])(id,1) = NaN;
                            else
                                pfa.(['floor_' num2str(n+1)])(id,1) = story.max_accel_x(n);
                            end
                        end
                    end
                else
                    error('NEED TO EXPLICITLY HANDLE MISSING STORY TABLE')
                end
            end
        end
    end

    % Convert to tables and sort
    idr_table = struct2table(idr);
    pfa_table = struct2table(pfa);
    idr_table_sort = [];
    pfa_table_sort = [];
    sa_vals = unique(idr_table.sa);
    for i = 1:length(sa_vals)
        for d = 1:2
            filt = idr_table.sa == sa_vals(i) & idr_table.direction == d;
            idr_table_sect = idr_table(filt,:);
            idr_table_sort = [idr_table_sort; idr_table_sect];
            pfa_table_sect = pfa_table(filt,:);
            pfa_table_sort = [pfa_table_sort; pfa_table_sect];
        end
    end
    
    % Remove ground motions that collapse (make sure there are not too many
    % cases)

    % Write tables to save directory
    write_dir = [main_dir '/' 'IDA' '/' 'EDPs'];
    mkdir(write_dir)
    writetable(idr_table_sort,[write_dir filesep 'idr.csv'])
    writetable(pfa_table_sort,[write_dir filesep 'pfa.csv'])
    max_a_ratio_table = struct2table(max_a_ratio);
    writetable(max_a_ratio_table,[write_dir filesep 'a_ratio.csv'])
    
    % Save P-58 input run data to remote location
    writetable(idr_table_sort,[model_remote_dir filesep 'idr.csv'])
    writetable(pfa_table_sort,[model_remote_dir filesep 'pfa.csv'])
    writetable(model,[model_remote_dir filesep 'model.csv'])
    writetable(story,[model_remote_dir filesep 'story.csv'])
    
    % Collect a_ratio data for ACI work
    num_gms_exceed_ac = [];
    num_gms = height(max_a_ratio_table);
    for s = 1:length(analysis.sa_stripes)
        gm_data = max_a_ratio_table{:,s+1};
        num_gms_exceed_ac(s) = sum(gm_data > 0.5 | isnan(gm_data));
    end
    p_exceed_ac(m,:) = num_gms_exceed_ac / num_gms;
end

% write collapse data for all models
writetable(col_data,[remote_dir filesep 'col_data.csv'])

% write ACI data
csvwrite([remote_dir filesep 'ACI_prob_exceedance.csv'],p_exceed_ac)

% plot ACI data
hold on
plt = plot(analysis.sa_stripes, p_exceed_ac,'-o');
legend(plt,strrep(model_data.name,'_',' '));
legend('location', 'northwest')
plot((2/3)*[ida_results.mce, ida_results.mce], [0, 1],'--k','HandleVisibility','off')
text(1.05*(2/3)*ida_results.mce, 0.03, 'DE','FontName', 'times', 'FontSize', 8)
plot([ida_results.mce, ida_results.mce], [0, 1],'--k','HandleVisibility','off')
text(1.05*ida_results.mce, 0.03, 'MCE_r','FontName', 'times', 'FontSize', 8)
upper_y = max(max(p_exceed_ac));
ylim([0 upper_y])
xlabel('Sa T_1 (g)')
ylabel('Probability of Exceeding 0.5*a')
box on
grid on
set(gcf,'position',[10,10,500,300])
set(gca,'fontname','times')
plt_name = 'ACI_plt';
savefig([remote_dir filesep plt_name '.fig'])
saveas(gcf,[remote_dir filesep plt_name],'png')
close