function results = TestCoefficientRestrictions(data,info)

% TestBlockExogeneity
%
% Usage:
%   results = TestCoefficientRestrictions(data,info);
%
% Purpose:
%   Tests zero coefficient restrictions in a vector autoregression using a 
%   likelihood ratio test. The null hypothesis is that a selected set of 
%   autoregressive coefficients is jointly zero.
%
% Input:
%   data    structure including
%   .y_sur      Tn x 1 vector, SUR endogenous variables
%   .X_sur      Tn x sum(k_i) array, SUR exogenous variables
%   .Y_ols      T x n matrix, OLS endogenous variables
%   .X_ols      T x k_i matrix, OLS exog. variables (k_i = k_j for all i,j)
%   info    structure including
%   .conv_crit  double, convergence criterion for SUR estimation
%   .max_it     integer, maximum number of iterations for SUR estimation
%   .rhs_names  n x max(k_i) matrix of strings, rhs variables per equation
%   .nexo       integer, number of exogenous terms m+x
%   .lag_order  integer, VAR lag order p
%   .nobs       integer, number of observations T
%   .nvars      integer, number of variables n
%   .names      1 x n vector of strings, variable names
% 
% Output:
%   results structure;
%   .nobs       integer, number of observations T-p
%   .alpha      double, significance level
%   .c_val      double, critical value of likelihood ratio test
%   .lik_ratio  double, likelihood ratio test statistic
%   .P_val      double, P-value for likelihood ratio test
%   .test_type  string, type of test ("coefficient_restrictions")
%
% Author:
%   Markus Kirchner, June 2012

% Get input
p         = info.lag_order;
nobs      = size(data.all,1)-p;
alpha     = info.alpha;
names     = info.names;
rhs_names = info.rhs_names;
nvars     = info.nvars;

% Impose the same rhs variables in each equation [no restrictions]
for i=1:nvars
    for j=1:nvars
        info.rhs_names{i,j} = names{i};
    end
end

% Run unrestricted VAR regression
results_u = EstimateVAR(data,info);

% Use different right-hand side variables in each equation [restrictions]
info.rhs_names = rhs_names;

% Run restricted VAR regression
results_r = EstimateVAR(data,info);

% Calculate degrees of freedom = number of restrictions
df = sum(sum(results_r.B==0));

% Number of parameters in each eq. of unrestr. system (see Enders, 1995)
% Note: set npar = 0 for no df adjustment (see Hamilton, 1994)
Bu   = results_u.B;
Cu   = results_u.C;
npar = size(Bu,1)+size(Cu,1);

% Calculate test statistic 
Omega_r   = results_r.Omega;
Omega_u   = results_u.Omega;
lik_ratio = (nobs - npar) * ( log(det(Omega_r)) - log(det(Omega_u)) );

% Obtain critical value 
c_val = chi2inv(1-alpha,df);

% Calculate P-value
P_val = 1 - chi2cdf(lik_ratio,df);

% Provide results
results.nobs      = nobs;
results.alpha     = alpha;
results.c_val     = c_val;
results.lik_ratio = lik_ratio;
results.P_val     = P_val;
results.df        = df;
results.test_type = 'coefficient_restrictions';