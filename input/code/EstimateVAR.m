function results = EstimateVAR(data,info)

% EstimateVAR
%
% Usage:
%   results = EstimateVAR(data,info);
%
% Purpose:
%   Reduced-form vector autoregression estimated by iterated feasible 
%   generalized least squares [FGLS]) or OLS. The regression model in SUR
%   form is given by
%       y_it = x_it' b_i + e_it     i = 1,...,n     t = 1,...,T 
%   Equivalently:
%       y_i = X_i*b_i + e_i
%   where y_i and e_i are T x 1 vectors, X_i is a T x k_i matrix, and 
%   b_i is a k_i x 1 vector. Equivalently:
%       y = X*b + e
%   where X = diag(X_1,X_2,...,X_n), b = vec(b_1,...,b_n), etc. To estimate
%   the model OLS is used when X_i = X_j for all i,j, otherwise FGLS is 
%   used because there are exclusion restrictions imposed on the right-hand
%   side variables. The OLS form of the model in matrix notation is
%       Y = X*B + E
%   where Y, X, B, and E are of size T x n, T x k, k x n, and T x n.
%
% Input:
%   data    structure including (accepts both SUR and OLS forms)
%   .y_sur      Tn x 1 vector, SUR endogenous variables
%   .X_sur      Tn x sum(k_i) array, SUR exogenous variables
%   .Y_ols      T x n matrix, OLS endogenous variables
%   .X_ols      T x k_i matrix, OLS exog. variables (k_i = k_j for all i,j)
%   info    structure including
%   .conv_crit  double, convergence criterion for SUR estimation
%   .max_it     integer, maximum number of iterations for SUR estimation
%   .rhs_names  n x max(k_i) matrix of strings, rhs variables per equation
%   .nexo       integer, number of exogenous terms m
%   .lag_order  integer, VAR lag order p
%   .nobs       integer, number of observations T
%   .nvars      integer, number of variables n
% 
% Output:
%   results structure;
%   .vecb       n*m+sum(k)*p x 1 vector, coefficients in vector form
%   .Sigb       n*m+sum(k)*p x n*m+sum(k)*p matrix, coef. var.-cov. matrix
%   .B          max(k_i)*p x n matrix, autoregr. coef. in simult. eqs. form
%   .A          n x n x p matrix, autoregr. coef. lag matrices
%   .C          m x n matrix, coefficients on exogenous variables
%   .Omega      n x n matrix, residual variance-covariance matrix
%   .fitted     T x n matrix, fitted values
%   .resid      T x n matrix, regression residuals
%   .LogL       double, system log likelihood
%   .aic        double, system Akaike information criterion
%   .sic        double, system Schwarz information criterion
%   .hqc        double, system Hannan-Quinn information criterion
%   .method     string, estimation method used
%
% Author:
%   Markus Kirchner, June 2012

% Get inputs
n         = info.nvars;
m         = info.nexo;
p         = info.lag_order;
conv_crit = info.conv_crit;
max_it    = info.max_it;
rhs_names = info.rhs_names;

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
B = zeros(n*p,n);
A = zeros(n,n,p);
C = zeros(m,n);

% Switch over estimation method
switch method

case 'OLS' % <------------- Use OLS to estimate the VAR parameters

% Get remaining inputs
Y_ols = data.Y_ols;
X_ols = data.X_ols;
T     = size(Y_ols,1);

% OLS regression equation by equation
Bhat = (X_ols' * X_ols) \ X_ols' * Y_ols;

% Put coefficients into vector form
vecb = reshape(Bhat,size(Bhat,1)*size(Bhat,2),1);

% Fitted values and residuals
fitted = X_ols * Bhat;
resid  = Y_ols - fitted;

% Residual variance-covariance matrix
mu    = repmat(mean(resid,1)',1,T);
Omega = (resid' - mu)  * (resid' - mu)' / T;

% Coefficients variance-covariance matrix
Sigb = kron( Omega, eye(length(vecb)/n)/(X_ols' * X_ols) );

case 'FGLS' % <------------- Use FGLS to estimate the VAR parameters

% Get remaining inputs
y_sur = data.y_sur;
X_sur = data.X_sur;
T     = length(y_sur)/n;

% Estimation step 1: OLS regression (equation by equation)
bhat(:,1) = (X_sur' * X_sur) \ X_sur' * y_sur;

% GLS iterations
i = 2; stop = 0; 
while ~stop

    % Compute residuals
    ehat = y_sur - X_sur * bhat(:,i-1);
    Ehat = reshape(ehat,T,n);
    
    % Compute residual variance-covariance matrix
    Omega = Ehat' * Ehat / T;
    
    % Estimation step 2: Generalized least squares regression
    tmp       = kron( eye(n) / Omega, eye(T) );
    bhat(:,i) = (X_sur' * tmp * X_sur) \ X_sur' * tmp * y_sur;
    
    % Convergence statistic
    stat = norm( bhat(:,i) - bhat(:,i-1) ) / norm( bhat(:,i-1) );
        
    % Check convergence
    if stat<conv_crit || i>=max_it
        stop = 1;
        if i>=max_it
            disp('EstimateVAR: FGLS did not converge.')
        end
    end
    
    % Go to next iteration
    i = i+1;
    
end

% Final coefficients
vecb = bhat(:,end);

% Fitted values and residuals
fitted = X_sur * vecb;
resid  = y_sur - fitted;

% Bring those into matrix form
fitted = reshape(fitted,T,n);
resid  = reshape(resid,T,n);

% Residual variance-covariance matrix
Omega = resid' * resid / T;

% Coefficients variance-covariance matrix
tmp  = kron(eye(n) / Omega, eye(T));
Sigb = eye(length(vecb)) / (X_sur' * tmp * X_sur);

end

% Bring coefficients into simult. eqs. form
ib = 1;
for ieq=1:n
    C(:,ieq) = vecb(ib:ib+m-1); % constant terms
    ib = ib + m;
    iB = 1;
    for ilag=1:p
        for jeq=1:n
            if ~isempty(char(rhs_names{ieq,jeq}))
                B(iB,ieq) = vecb(ib); % autoregressive coefficients
                ib = ib + 1;
            end
            iB = iB + 1;
        end
    end
end

% Generate lag matrices
Bprime = B';
j      = 1;
for i=1:n:1+n*(p-1)
    if j>p
        break
    end
    A(:,:,j) = Bprime(:,i:i+n-1);
    j = j+1;
end

% Log likelihood
LogL = -.5 * T * n * ( 1 + log( 2 * pi ) ) -.5 * T * log(real(det(Omega)));

% Total number of parameters
k = sum(sum(B~=0)) + size(C,1) * size(C,2);

% Information criteria
aic = -2 * LogL / T + 2 * k / T;
sic = -2 * LogL / T + k * log(T) / T;
hqc = -2 * LogL / T + 2 * k * log(log(T)) / T;

% Provide output
results.B      = B;
results.A      = A;
results.C      = C;
results.vecb   = vecb;
results.Sigb   = Sigb;
results.Omega  = Omega;
results.fitted = fitted;
results.resid  = resid;
results.LogL   = LogL;
results.aic    = aic;
results.sic    = sic;
results.hqc    = hqc;
results.method = method;