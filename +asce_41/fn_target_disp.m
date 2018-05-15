function [ target_disp_in ] = fn_target_disp( strength_ratio, site_class, mode_shape, num_stories, period, sa, asce_41_13_cm )
% FUNCTION DESCRIPTION
% The funciton compute_target_roof_displacement follows section 3 in the
% document "SP3 Implementation of the P-58 Simplified Method" and uses 
% values input by the user in order to calculate the fisrt mode elastic 
% target roof displacement according to ASCE 41-06.

%   Created By: Dustin Cook
%   Date Created: 10/22/15
%   Date Last Modified: 1/7/16

% FUNCTION INPUTS 
%   S: Strength Ratio, as descibed by the P-58 simplified method, for one 
%      level of Sa. Float value. Unitless ratio.
%   siteClass: String that descibes the building site class. String is one
%              capitol letter A through F.
%   modeShape: First deflected mode shape of the structure. Single        
%              dimensional array that has the length of number of stories
%              plus one. Mode shape should be normalized at the roof and
%              the first value in the array should be zero. Unitless.
%   numStories: Total number of stories in the structure. Float value.
%               Units of number of stories.
%   T: First mode building period. Float value with units os seconds.
%   SA: Spectral Accelation value for one level of hazard. Float value 
%       with units of g. 
%   buildingType: String that describes the structural system.
%                 'mf' = Moment Frame
%                 'lf' = Light Frame
%                 'sw' = Shear Wall
%                 'mf_sw" = Dual System: Moment Frame and Shear Wall
%                 'cbf' = Concentrically Braced Frame
%                 'ebf' = Eccentrically Braced Frame
%                 'brb' = Buckling Restrained Braced Frame
%
% FUNCTION OUTPUTS
%   TD: Target Displacement, as descibed by equation 3-14 of ASCE 41-06, 
%       for each level of Sa. Float value with units of inches.

%% Compute the Modal Participation Factor 
G1 = sum(mode_shape)/sum(mode_shape.^2);

%% Determine the Modal Mass Factor (Cm)
if period<=1 && num_stories >=3
    CM = asce_41_13_cm;
else
    CM = 1;
end

%% Compute Inelastic Modification Factor (C1)
% Determine a factor relating to soil conditions
switch site_class 
    case 'A' 
        a = 130;
    case 'B'
        a = 130;
    case 'C' 
        a = 90;
    case 'D' 
        a = 60;
    case 'E' 
        a = 60;
    case 'F'
        a = 60;
    otherwise 
        error('ASCE41:InvalidCase','Unexpected Site Class. Check Inputs')
end

% Calculte R factor
R = strength_ratio*CM;

% Calculate C1 based on period (T) 
if period < 0.2 
    c1 = 1+(R-1)/(0.04*a);
elseif period <= 1
    c1 = 1+(R-1)/(a*period^2);
else
    c1 = 1;
end

%% Compute hysteretic Modification Factor (C2)
if period < 0.2  
    c2 = 1+((R-1).^2)/32;
elseif period <= 0.7
    c2 = 1+(1/800)*((R-1).^2)/(period^2);
else
    c2 = 1;
end

%% Compute Tareget Displacement using equation 3-14 from ASCE 41-06
g = 386.4; % acceleration of gravity in in/s^2
target_disp_in = G1*c1*c2*sa*g*period^2/(4*pi()^2);

end

