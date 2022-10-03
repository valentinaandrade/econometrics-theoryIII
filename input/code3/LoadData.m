function [data info] = LoadData(info)

% LoadData
%
% Usage:
%   [data info] = LoadData(info);
%
% Purpose:
%   Loads data for a reduced-form VAR of order p.
%
% Input:
%   info    structure including
%   .file       string, source file with data (including path)
%   .first_obs  double, first observation in sample
%   .last_obs   double, last observation in sample
% 
% Output:
%   data    T+p x N matrix, data with obs./variables in rows/columns
%   info    structure;
%   .names      1 x N vector of strings, variable names   
%   .dates      T+p x 1 vector, dates corresponding to sample
%
% Author:
%   Markus Kirchner, April 2012

% Get input
file       = info.file;
first_obs  = round2(info.first_obs,4);
last_obs   = round2(info.last_obs,4);

% Load data from file
[raw_data txt] = xlsread(file);
var_names      = txt(2:end);

% Get the sample to be used
dates       = round2(raw_data(:,1),4); 
obs_id      = 1:size(raw_data,1);
sample      = obs_id(dates==first_obs):obs_id(dates==last_obs);
select_data = raw_data(sample,2:end);
dates       = dates(sample);

% Provide output
data       = select_data;
info.dates = dates;
info.names = var_names;