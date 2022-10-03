function results = TestLagLength(data,info)

% TestLagLength
%
% Usage:
%   results = TestLagLength(data,info);
%
% Purpose:
%   Tests for appropriate lag length in a VAR, allowing for general 
%   exclusion restrictions. For a given maximum lag order z, z+1 
%   likelihood ratio tests are implemented, each  testing for inclusion 
%   of lag order j-1 vs. j. The null hypothesis is that j should not be 
%   included as additional lag. In addition, information criteria (IC)
%   are reported for each individual lag length. The lowest value 
%   of the IC indicates the preferred lag length.
%
% Input:
%   data    structure;
%   .all        T+p x n matrix, data with obs./variables in rows/columns
%   .exo        m x n matrix, exogenous variables and deterministic terms
%   info    structure including
%   .lag_order  integer, VAR lag order p
%   .max_lag    integer, maximum lag order z to be tested
%   .alpha      double, significance level
%   .rhs_names  n x max(k_i) matrix of strings, rhs variables per equation
%   .lhs_names  1 x n vector of strings, left-hand side variables 
%
% Output:
%   results structure;
%   .nobs       integer, number of observations T-max_lag
%   .alpha      double, significance level
%   .c_val      double, critical value of likelihood ratio test
%   .lik_ratio  1 x z+1 vector, likelihood ratio
%   .P_val      1 x z+1 vector, P-values for likelihood ratio tests
%   .aic        1 x z+1 vector, Akaike information criterion
%   .hqc        1 x z+1 vector, Hannan-Quinn information criterion
%   .sic        1 x z+1 vector, Schwarz information criterion
%   .lags       1 x z+1 vector, lags tested
%   .test_type  string, type of test ("lag_length")
%
% Author:
%   Markus Kirchner, May 2012

% Ignore warnings, similarly as in EViews
warning off;

% Get input
alpha     = info.alpha;
max_lag   = info.max_lag;
[T n]     = size(data.all);
rhs_names = info.rhs_names;
lhs_names = info.lhs_names;
nexo      = info.nexo;

% Prepare check of estimation method
ols_check = zeros(size(rhs_names));
for i=1:size(rhs_names,1)
   for j=1:size(rhs_names,2)
       ols_check(i,j) = ~isempty(char(rhs_names{i,j}));
   end
end

% Check whether OLS or SUR should be used
if all(all(ols_check))
    method = 'OLS';  % use OLS when there are no zero restrictions
else
    method = 'FGLS'; % use FGLS when there are zero restrictions
end

% Allocate memory
Omega     = zeros(n,n,max_lag+1);
lik_ratio = NaN(1,max_lag+1);
P_val     = NaN(1,max_lag+1);
LogL      = zeros(1,max_lag+1);
aic       = zeros(1,max_lag+1);
hqc       = zeros(1,max_lag+1);
sic       = zeros(1,max_lag+1);
   
% Estimate various VARs
for p=0:max_lag

    % Current lag order
    info.lag_order = p;

    % Adjust sample to use same data in each run
    begin = 1 + max_lag - p;
    
    % Adjust data to this sample
    estim_data.all = data.all(begin:end,:);
    estim_data.exo = data.exo(begin:end,:);
    
    % Switch over estimation method
    switch method

    case 'FGLS' % <------------- Prepare data for FGLS estimation
        
    % Left-hand side variables
    this_data.y_sur = vec(estim_data.all(1+p:end,:));

    % Right-hand side variables
    this_data.X_sur = [];
    for i=1:size(rhs_names,1)
        X_i = estim_data.exo(1+p:end,:);  % exogenous terms
        q = 1;
        while q<=p      % go through lags
            for j=1:size(rhs_names,2)   % check rhs variables
                if (~isempty(char(rhs_names{i,j})))
                    pos_x = loc(char(lhs_names),char(rhs_names(i,j)));
                    X_i = [X_i estim_data.all(1+p-q:end-q,pos_x)];
                end
            end
            q = q+1;
        end
        this_data.X_sur = blkdiag(this_data.X_sur,X_i);
    end
    
    case 'OLS' % <------------- Use OLS to estimate the VAR parameters
    
    % Left-hand side variables
    this_data.Y_ols = estim_data.all(1+p:end,:); 
   
    % Right-hand side variables
    this_data.X_ols = estim_data.exo(1+p:end,:);  % exogenous terms
    q = 1;
    while q<=p      % go through lags
        this_data.X_ols = [this_data.X_ols estim_data.all(1+p-q:end-q,:)];
        q = q+1;
    end

    end
    
    % Estimate VAR
    output = EstimateVAR(this_data,info);

    % Store relevant information
    Omega(:,:,p+1) = output.Omega;
    LogL(p+1)      = output.LogL;
    aic(p+1)       = output.aic;
    sic(p+1)       = output.sic;
    hqc(p+1)       = output.hqc;
    
end 

% Degrees of freedom = number of restrictions
df = sum( sum( (output.A(:,:,1)~=0) ) );

% Critical value LR test
c_val = chi2inv(1-alpha,df);

% LR test statistic and P-value
for i=2:max_lag+1
    
    % Current lag order unrestricted system
    p = i-1;
        
    % Number of parameters in longest unres. eq. (see Canova, 2007)
    % Note: set npar = 0 for no df adjustment
    npar = p*n+nexo;

    % LR test statistic
    lik_ratio(i) = (T - max_lag - npar) * ...
        ( log(det(Omega(:,:,i-1))) - log(det(Omega(:,:,i))) );

    % P-value LR test
    P_val(i) = 1 - chi2cdf(lik_ratio(i),df);
    
end

% Provide results
results.nobs      = T-max_lag;
results.alpha     = alpha;
results.c_val     = c_val;
results.LogL      = LogL;
results.lik_ratio = lik_ratio;
results.P_val     = P_val;
results.aic       = aic;
results.sic       = sic;
results.hqc       = hqc;
results.lags      = 0:max_lag;
results.test_type = 'lag_length';