function [spectra, periods] = fn_call_USGS_hazard_API(edition, lat, lng, vs30, afe)
% Description: Function to Call the USGG Design Values API and return,
% code short and 1 second spectra values for the base line code, risk
% targeted approach, and uniform hazard approach

% Created by: Dustin Cook
% Date Created: 7/26/2021

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

% Set parameters
options = weboptions('Timeout', 30);         
periods = [0.1; 0.2; 0.3; 0.5; 0.75; 1; 2; 3]; % only these periods return values from the API call
periods_str = {'0P1' '0P2' '0P3' '0P5' '0P75' '1P0' '2P0' '3P0'}; % only these periods return values from the API call


try
    %% Call USGS Hazards API
    for i = 1:length(periods)
        DATA = webread(['https://earthquake.usgs.gov/nshmp-haz-ws/hazard/' edition ...
            '/WUS/' num2str(lng) '/' num2str(lat) '/SA'  periods_str{i}  '/'...
            num2str(vs30)], options);
        
        sa_vals = DATA.response.metadata.xvalues;
        afe_vals = DATA.response.data.yvalues;
        crop_filt = afe_vals > 0;
        spectra(i,:) = interp1(afe_vals(crop_filt),sa_vals(crop_filt),afe);
    end
catch
    error('Problem grabbing site hazards')
end

end

