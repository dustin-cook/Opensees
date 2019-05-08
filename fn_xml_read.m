function [ data ] = fn_xml_read( file_name )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Import Packages
import file_exchange.*

%% Begin Method
xml_read = xmlread(file_name);
[ xml_struct ] = xml2struct( xml_read );
clear xml_read
% filtered_txt = strrep(xml_struct.OpenSees.Data.Text,'-1.#QNAN','0');
% filtered_txt = strrep(filtered_txt,'-1.#IND','0');
data = str2num(xml_struct.OpenSees.Data.Text);

end

