function [design_values] = fn_call_USGS_design_value_API(id, reference_doc, lat, lng, risk_cat, site_class)
% Description: Function to Call the USGG Design Values API and return,
% code short and 1 second spectra values for the base line code, risk
% targeted approach, and uniform hazard approach

% Created by: Dustin Cook
% Date Created: 2/15/2019

% Inputs:
%   id - numeric or integer unique identifier for the site
%   reference_doc - string represention of the code edition to be called
%                   options include [ asce7-05, asce7-10, asce7-16, 
%                                     nehrp-2009, nehrp-2015, nehrp-2020
%                                     asce41-13, asce41-17, and others]
%   lat - numeric lattitide of the site
%   lng - numeric longitude of the site
%   risk_cat - string representation of the risk category (I, II, III or IV)
%   site_class - sting representation of the site class (A, B, C, D, E, or F)

% Ouputs:
%   design_values.ss = short period code spectral acceleration
%   design_values.s1 = s1 second period code spectral acceleration
%   design_values.ssrt = short period risk targeted spectral acceleration
%   design_values.s1rt = s1 second period risk targeted spectral acceleration
%   design_values.ssuh = short period uniform hazard spectral acceleration
%   design_values.s1uh = s1 second period uniform hazard spectral acceleration

try
    %% Call USGS API
    options = weboptions('Timeout', 30);
    DATA = webread( [ 'https://earthquake.usgs.gov/ws/designmaps/' ...
                      reference_doc '.json?' ], ...
                      'latitude', lat, ...
                      'longitude', lng, ...
                      'riskCategory', risk_cat, ...
                      'siteClass', site_class, ...
                      'title', num2str(id), ...
                      options);

    %% Post Process Return
    design_values.ss = DATA.response.data.ss;
    design_values.s1 = DATA.response.data.s1;

    % Check if Fa is empty and why
    if isempty(DATA.response.data.fa)
        design_values.fa = NaN;
        design_values.sms = NaN;
        design_values.sds = NaN;
    else
        design_values.fa = DATA.response.data.fa;
        design_values.sms = DATA.response.data.sms;
        design_values.sds = DATA.response.data.sds;
    end

    % Check if Fv is empty and why
    if isempty(DATA.response.data.fv)
        design_values.fv = NaN;
        design_values.sm1 = NaN;
        design_values.sd1 = NaN;
    else
        design_values.fv = DATA.response.data.fv;
        design_values.sm1 = DATA.response.data.sm1;
        design_values.sd1 = DATA.response.data.sd1;
    end

    % Check if SDC is empty
    if isempty(DATA.response.data.sdc)
        design_values.sdc = NaN;
    else
        design_values.sdc = DATA.response.data.sdc;
    end


    % Check to see if risk targeted values exist
    if isfield(DATA.response.data,'ssrt') && isfield(DATA.response.data,'s1rt')
        design_values.ssrt = DATA.response.data.ssrt;
        design_values.s1rt = DATA.response.data.s1rt;
    else
        design_values.ssrt = NaN;
        design_values.s1rt = NaN;
    end

    % Check to see if Uniform Hazard values exist
    if isfield(DATA.response.data,'ssuh') && isfield(DATA.response.data,'s1uh')
        design_values.ssuh = DATA.response.data.ssuh;
        design_values.s1uh = DATA.response.data.s1uh;
    else
        design_values.ssuh = NaN;
        design_values.s1uh = NaN;
    end
    
catch
    % If the call fails
    design_values.ss = NaN;
    design_values.s1 = NaN;
    design_values.fa = NaN;
    design_values.sms = NaN;
    design_values.sds = NaN;
    design_values.fv = NaN;
    design_values.sm1 = NaN;
    design_values.sd1 = NaN;
    design_values.sdc = NaN;
    design_values.ssrt = NaN;
    design_values.s1rt = NaN;
    design_values.ssuh = NaN; 
    design_values.s1uh = NaN;
end

end

