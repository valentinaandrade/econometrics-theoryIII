function output = ComputeShockDecomposition(input)

% ComputeShockDecomposition
%
% Usage:
%   output = ComputeShockDecomposition(input);
%
% Purpose:
%   Provides structural shocks for a vector autoregression with fixed
%   parameters based on a Cholesky decomposition.
%
% Input:
%   input   structure including
%   .resid      T x n matrix, regression residuals
%   .Omega      n x n matrix, residual variance-covariance matrix
%
% Output:
%   output  structure "input" including
%   .shocks     T x n matrix, structural shocks
%
% Author:
%   Markus Kirchner, June 2012

% Get inputs
resid = input.resid;
Omega = input.Omega;
[T n] = size(resid);

% Allocate memory
shocks = zeros(n,T);

% Contemp. relations
A0 = eye(n) / chol(Omega)';

% Compute structural shocks
for i=1:T
    shocks(:,i) = A0 * resid(i,:)';
end

% Prepare output
output = input;

% Provide output
output.shocks = shocks';