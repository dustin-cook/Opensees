function [ data ] = fn_xml_read( file_name )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import file_exchange.*

%% XML to Struct File Exchange Method
xml_read = xmlread(file_name);
[ xml_struct ] = xml2struct( xml_read );
data = str2num(xml_struct.OpenSees.Data.Text);

%% Reg Exp Method
% text = fileread(file_name);
% expr = '(?<=<Data>).*';
% text = regexp(text,expr,'match');
% expr = '.*(?=</Data>)';
% text = regexp(text,expr,'match');
% expr = '\s*\S*';
% matches = regexp(text{1},expr,'match');
% data = str2double(matches{1});
% data = data(1:(end-1)); % remove last NaN

end

