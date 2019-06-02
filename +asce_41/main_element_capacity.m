function [ element, joint ] = main_element_capacity( story, ele_prop_table, element, analysis, joint, read_dir, write_dir )
% Description: Main script that calculates the strength of each element in 
% the model according to ASCE 41 

% Created By: Dustin Cook
% Date Created: 1/3/2019

% Inputs:

% Outputs:

% Assumptions:


%% Import Packages
import asce_41.*

%% Begin Method
for i = 1:length(element.id)
    ele = element(i,:);
    ele_id = ele.ele_id;
    ele_prop = ele_prop_table(ele_prop_table.id == ele_id,:);
    TH_file = [read_dir filesep 'element_TH_' num2str(ele.id) '.mat'];
    if ~exist(TH_file,'file')
        ele_TH = [];
    else
        load(TH_file)
    end
    
    %% Calculate Element Capacties
    if isempty(ele_TH)
        [ ele, ~, ~ ] = fn_element_capacity( story, ele, ele_prop, ele_TH, analysis.nonlinear );
        ele_PM = [];
    else
        [ ele, ele_TH, ele_PM ] = fn_element_capacity( story, ele, ele_prop, ele_TH, analysis.nonlinear );
        % Save time histories and clear data
        save([write_dir filesep 'element_TH_' num2str(ele.id) '.mat'],'ele_TH')
        save([write_dir filesep 'element_PM_' num2str(ele.id) '.mat'],'ele_PM')
    end
    
    % Display Progress
    if ~analysis.suppress_outputs
        disp([num2str(i), ' out of ', num2str(length(element.id)) ' elements complete' ])
    end
    
    %% Caculate required development length and make sure there is enough
    if ele.id == 17
        test = 5;
    end
    [ ele.pass_aci_dev_length, ele.ld_avail, ele.ld_req ] = fn_development_check( ele, ele_prop );
    
    % Save capcity for convergence check (just use onesided capacity for
    % now)
    if strcmp(ele.type,'wall')
        ele.capacity = ele.Vn_1;
    else
        ele.capacity = ele.Mn_pos_1;
    end
    
    % Save data to main table
    ele_to_save(i,:) = ele;
end

element = ele_to_save;

%% Calculate Joint Properties;
if ~isempty(joint)
    [ joint ] = fn_joint_capacity( joint, element, ele_prop_table );
end

end

