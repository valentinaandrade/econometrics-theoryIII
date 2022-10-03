function output = ComputeHistoricalDecomposition(input,info,data)

% ComputeHistoricalDecomposition
%
% Usage:
%   output = ComputeHistoricalDecomposition(input,info,data);
%
% Purpose:
%   Provides historical decomposition for a vector autoregression with
%   fixed parameters based on a Cholesky decomposition.
%
% Input:
%   input   structure including
%   .A          n x n x p matrix, autoregr. coef. lag matrices
%   .C          m x n matrix, coefficients on exogenous variables
%   .Omega      n x n matrix, residual variance-covariance matrix
%   info    structure including
%   .lag_order  integer, VAR lag order p 
%   .nexo       integer, number of exogenous terms m
%   data    structure including
%   .all        T+p x n matrix, data with obs./variables in rows/columns
%   .exo        T+p x M matrix, deterministic terms

%
% Output:
%   output  structure "input" including
%   .hist_dec   T+p x n x n+2 array, contribution of shocks and exog. terms
%               rows: periods; columns: variables; 3rd: n shocks, then
%               init. values, then exog. 
%
% Author:
%   Markus Kirchner, July 2012

% Get inputs
A         = input.A;
C         = input.C;
Omega     = input.Omega;
shocks    = input.shocks;
[T n]     = size(shocks);
p         = info.lag_order;
m         = info.nexo;
all_data  = data.all;
exo_data  = data.exo;

% Allocate memory
decomp = zeros(T+p,n,n+2);

% Contribution of pre-sample values
for i=1:T+p
    if i<=p
        % Pre-sample values
        decomp(i,:,n+1) = all_data(i,:);
    else
        % Autoregressive parts
        for j=1:p
            decomp(i,:,n+1) = decomp(i,:,n+1) + decomp(i-j,:,n+1) * A(:,:,j)';
        end
    end
end

% Contribution of exog. terms
for i=1:T
    % Exogenous terms
    if m>0
        decomp(i+p,:,n+2) = exo_data(i+p,:) * C;
    end
    % Autoregressive parts
    for j=1:p
        decomp(i+p,:,n+2) = decomp(i+p,:,n+2) + decomp(i+p-j,:,n+2) * A(:,:,j)';
    end
end

% Cholesky factor
cho_fac = chol(Omega)';

% Contribution of shocks
for k=1:n
    % Current shock
    this_shock = zeros(1,n);
    for i=1:T
        % Autoregressive parts
        for j=1:p
            decomp(i+p,:,k) = decomp(i+p,:,k) + decomp(i+p-j,:,k) * A(:,:,j)';
        end
        % Current shock, setting others to zero
        this_shock(k) = shocks(i,k);
        % Shocks
        decomp(i+p,:,k) = decomp(i+p,:,k) + this_shock * cho_fac';
    end
end

% Prepare output
output = input;

% Provide output
output.hist_dec = decomp;