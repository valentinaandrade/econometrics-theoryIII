function [data info] = PrepareRegressionData(raw_data,info)

% LoadData
%
% Usage:
%   [data info] = PrepareRegressionData(raw_data,info);
%
% Purpose:
%   Prepares data for a reduced-form vector autoregression that allows
%   for general exclusion restrictions. The data is put into SUR form and
%   OLS form to be used later on according to whether restrictions are
%   imposed or not.
%
% Input:
%   raw_data    T+p x N matrix, data with obs./variables in rows/columns
%   info    structure including
%   .names      1 x N vector of strings, variable names 
%   .lhs_names  1 x n vector of strings, left-hand side variables 
%   .rhs_names  n x max(k_i) matrix of strings, rhs variables per equation
%   .exo_names  1 x m vector of strings, exogenous variables
%   .do_cons    boolean, if 1 add constant, else 0
%   .do_lin     boolean, if 1 add linear trend, else 0
%   .do_quad    boolean, if 1 add linear-quadratic trend, else 0
%   .do_cub     boolean, if 1 add linear-quadratic trend, else 0
%   .sea_freq   integer, seasonal frequency of data set 
%               (4 = quarterly, 12 = monthly, 0 = no seasonal dummies)
%   .lag_order  integer, VAR lag order p
% 
% Output:
%   data    structure;
%   .all        T+p x n matrix, data with obs./variables in rows/columns
%   .exo        T+p x m matrix, exogenous variables and deterministic terms
%   .y_sur      Tn x 1 vector, SUR endogenous variables
%   .X_sur      Tn x sum(k_i) array, SUR exogenous variables
%   .Y_ols      T x n matrix, OLS endogenous variables
%   .X_ols      T x k_i matrix, OLS exog. variables (k_i = k_j for all i,j)
%   info    structure including
%   .nexo       integer, number of exogenous terms M
%   .nobs       integer, number of observations T
%   .nvars      integer, number of variables n
%   .names      1 x n vector of strings, variable names
%
% Author:
%   Markus Kirchner, April 2012

% Get input
var_names = info.names;
exo_names = info.exo_names;
lhs_names = info.lhs_names;
rhs_names = info.rhs_names;
n_vars    = numel(lhs_names);
d         = info.sea_freq;
do_cons   = info.do_cons;
do_lin    = info.do_lin;
do_quad   = info.do_quad;
do_cub    = info.do_cub;
p         = info.lag_order;

% Check inputs
if (~isempty(char(info.exo_names{:})))
    n_exo  = numel(exo_names);
else
    n_exo  = 0;
end

% Get data information
n_obs = size(raw_data,1);

% Allocate memory
endo_data = zeros(n_obs,n_vars);
exo_data  = zeros(n_obs,n_exo);
pos_exo   = zeros(1,n_exo);
pos_endo  = zeros(1,n_vars);
cons      = []; 
trend     = []; 
dums      = [];

% Locate positions of variables in Y_t
for i=1:n_vars
    pos_endo(i) = loc(char(var_names),char(lhs_names(i))); 
    names(i)    = var_names(pos_endo(i));
end

% Locate positions of variables in X_t
for i=1:n_exo
    pos_exo(i) = loc(char(var_names),char(exo_names(i)));
end

% Prepare left-hand side variables
for i=1:n_vars
    endo_data(:,i) = raw_data(:,pos_endo(i));
end

% Stack data into a single vector for SUR form
y_sur = reshape(endo_data(1+p:end,:),(n_obs-p)*n_vars,1); 

% Put data into a matrix for OLS form
Y_ols = endo_data(1+p:end,:);

% Get exogenous variables
for i=1:n_exo
    exo_data(:,i) = raw_data(:,pos_exo(i));
end

% Check deterministic terms
if do_cons
    cons = ones(n_obs,1);
end
if do_lin
    trend = (1:n_obs)';
end
if do_quad
    trend = (1:n_obs)';
    trend = [trend trend.^2];
end
if do_cub
    trend = (1:n_obs)';
    trend = [trend trend.^2 trend.^3];
end
if d~=0
   dums        = repmat(eye(d),round(n_obs/d+1),1);   
   dums(:,end) = [];                               
   while size(dums,1)>n_obs
       dums(end,:) = [];
   end
end

% Add deterministic terms
exo_data = [cons trend dums exo_data];

% Number of exogenous terms
n_exo = size(exo_data,2);

% Prepare right-hand side variables for SUR form
X_sur = [];
for i=1:size(rhs_names,1)
    X_i = exo_data(1+p:end,:);  % exogenous terms
    k=1;
    while k<=p      % go through lags
        for j=1:size(rhs_names,2)   % check rhs variables
            if (~isempty(char(rhs_names{i,j})))
                pos_x = loc(char(var_names),char(rhs_names(i,j)));
                X_i = [X_i raw_data(1+p-k:end-k,pos_x)];
            end
        end
        k = k+1;
    end
    X_sur = blkdiag(X_sur,X_i);
end

% Prepare right-hand side variables for OLS form
X_ols = exo_data(1+p:end,:);  % exogenous terms
k=1;
while k<=p      % go through lags
    for i=1:n_vars  
        pos_x = loc(char(var_names),char(lhs_names(i)));
        X_ols = [X_ols raw_data(1+p-k:end-k,pos_x)];
    end
    k = k+1;
end

% Provide output
data.all   = endo_data;
data.exo   = exo_data;
data.y_sur = y_sur;
data.X_sur = X_sur;
data.Y_ols = Y_ols;
data.X_ols = X_ols;
info.nobs  = n_obs-p;
info.nvars = n_vars;
info.nexo  = n_exo;
info.names = names;