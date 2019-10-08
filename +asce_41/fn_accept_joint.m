function [joint] = fn_accept_joint(joint,joint_explicit,read_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for i = 1:height(joint)
    if joint_explicit
        TH_file = [read_dir filesep 'joint_TH_' num2str(joint.id(i)) '.mat'];
        if exist(TH_file,'file')
            load(TH_file)
            jnt_yeild_rot = joint.Mn(i)/(101*6*joint.e(i)*joint.iz(i)/joint.h(i));
            joint.max_rotation(i) = max(abs(jnt_TH.deformation_TH));
            if joint.max_rotation(i) <= joint.io(i) + jnt_yeild_rot
                joint.accept(i) = 1; % Passes IO
            elseif joint.max_rotation(i) <= joint.ls(i) + jnt_yeild_rot
                joint.accept(i) = 2; % Passes LS
            elseif joint.max_rotation(i) <= joint.cp(i) + jnt_yeild_rot
                joint.accept(i) = 3; % Passes CP
            else
                joint.accept(i) = 4; % Fails all Acceptance Criteria
            end
        end
    else
        % Check force controlled acceptance of joints (assumes joints are
        % critical elements)
        force_control_demands_CP = 1*1.3*(joint.Vmax(i) - joint.V_grav(i)) + joint.V_grav(i);
        force_control_demands_LS = 1.3*1.3*(joint.Vmax(i) - joint.V_grav(i)) + joint.V_grav(i);
        force_control_demands_IO = force_control_demands_LS;

        if force_control_demands_IO/joint.Vj(i) < 1.0
            joint.accept(i) = 1; % Passes IO
        elseif force_control_demands_LS/joint.Vj(i) < 1.0
            joint.accept(i) = 2; % Passes LS
        elseif force_control_demands_CP/joint.Vj(i) < 1.0
            joint.accept(i) = 3; % Passes CP
        else
            joint.accept(i) = 4; % Fails all Acceptance Criteria
        end
    end
end

end

