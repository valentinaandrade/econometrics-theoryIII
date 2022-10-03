function stable = CheckStability(vecb,info)

% CheckStability
%
% Usage:
%   stable = CheckStability(vecb,info);
%
% Purpose:
%   Checks stability of VAR coefficients.
%
% Input:
%   vecb    k x 1 vector, draw of coefficient states for period t
%   info    structure including
%   .nexo       integer, number of exogenous terms in the VAR
%   .nvars      integer, number of variables
%   .lag_order  integer, VAR lag order p
%   .rhs_names  n x max(k_i) matrix of strings, rhs variables per equation
%  
% Output:
%   stable  boolean, 1 if coefficients are stable, else 0
%
% Author:
%   Markus Kirchner, August 2012

% Bring coefficients into simult. eqs. form
B = GenerateSimultEqsMatrices(vecb,info);
        
% Generate VAR(p) lag matrices
A = GenerateLagMatrices(B);

% Generate VAR(1) companion matrix
F = GenerateCompanionMatrix(A);

% Eigenvalue check
stable = 1;
if max(abs(eig(F))>=1)
    stable = 0; 
end