function [ data ] = fn_xml_read( file_name )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import file_exchange.*

%% Begin Method
xml_read = xmlread(file_name);
[ xml_struct ] = xml2struct( xml_read );
data = str2num(xml_struct.OpenSees.Data.Text);

end

