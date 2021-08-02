function [ element ] = fn_factor_loads( analysis, element, site )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if analysis.type == 4
    % Factor for design (ie static loading)
    element.gravity_load =   (analysis.eq_vert_load_factor*site.Sds + analysis.dead_load)*element.dead_load + ...
                              analysis.live_out_load*element.live_load.*~element.inner_bay + analysis.live_in_load*element.live_load.*element.inner_bay;
    element.gravity_load_1 = (analysis.eq_vert_load_factor*site.Sds + analysis.dead_load)*element.dead_load_1 + ...
                              analysis.live_out_load*element.live_load_1.*~element.inner_bay + analysis.live_in_load*element.live_load_1.*element.inner_bay;
    element.gravity_load_2 = (analysis.eq_vert_load_factor*site.Sds + analysis.dead_load)*element.dead_load_2 + ...
                              analysis.live_out_load*element.live_load_2.*~element.inner_bay + analysis.live_in_load*element.live_load_2.*element.inner_bay;
else
    % No factors for all else
    element.gravity_load =   analysis.dead_load*element.dead_load   + analysis.live_load*element.live_load;
    element.gravity_load_1 = analysis.dead_load*element.dead_load_1 + analysis.live_load*element.live_load_1;
    element.gravity_load_2 = analysis.dead_load*element.dead_load_2 + analysis.live_load*element.live_load_2;
end

end

