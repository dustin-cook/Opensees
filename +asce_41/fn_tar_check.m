function [ TAR ] = fn_tar_check( max_disp_profile, ave_disp_profile )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

TAR = max_disp_profile(2:end) ./ ave_disp_profile(2:end);
if sum(TAR > 1.2) > 0
    warning('Torsional Amplification Exceeds Limit. Forces and displacements caused by accidental torsion shall be amplified by a factor Ax')
end

end

