function [bands sig] = ComputeConfidenceBands(input,info,data)

% ComputeConfidenceBands
%
% Usage:
%   [bands sig] = ComputeConfidenceBands(input,info,data);
%
% Purpose:
%   Computes bootstrapped confidence bands for impulse responses of a 
%   VAR with general exclusion restrictions.
%
% Input:
%   input   structure including
%   .A          n x n x p matrix, autoregr. coef. lag matrices
%   .C          m x n matrix, coefficients on exogenous variables
%   .resid      (T-p) x n matrix, regression residuals
%   .method     string, estimation method used
%   info    structure including
%   .lag_order  integer, VAR lag order p
%   .horizon    scalar, horizon h for impulse response functions
%   .rep        scalar, number of bootstrap replications
%   .conf       1 x c vector, confidence levels for error bands
%   .nobs       integer, number of observations of data set
%   .rhs_names  n x max(k_i) matrix of strings, rhs variables per equation
%   .lhs_names  1 x n vector of strings, left-hand side variables 
%   data    structure;
%   .all        T+p x n matrix, data with obs./variables in rows/columns
%   .exo        m x n matrix, exogenous variables and deterministic terms
%
% Output:
%   bands   n x h+1 x 2 x c array, lower (1st) and upper (2nd) bands
%   sig     n x h+1 matrix, IRF standard deviations (optional)
%
% Author:
%   Markus Kirchner, May 2012

% Provide inputs
A         = input.A;
C         = input.C;
p         = info.lag_order;
T         = size(data.all,1);
resid     = input.resid;
rep       = info.rep;
h         = info.horizon;
m         = info.nexo;
conf      = info.conf;
c         = numel(conf);
exo_data  = data.exo;
rhs_names = info.rhs_names;
lhs_names = info.lhs_names;
n         = length(lhs_names);
method    = input.method;

% Allocate memory
boot_irf  = zeros(n,h+1,rep);
bands     = zeros(n,h+1,2);
boot_data = zeros(T,n,rep);

% Generate bootstrap indices
[~,index] = bootstrp(rep,[],resid);    

% Put waitbar
h_wait = waitbar(0,'Generating bootstrap samples, please wait...');

% Generate bootstrap samples
for k=1:rep
    waitbar(k/rep,h_wait)   % update waitbar
    for i=1:T
        if i<=p
            boot_data(i,:,k) = data.all(i,:);  % use pre-sample values here
        else
            if m>0
                boot_data(i,:,k) = data.exo(i,:)*C; % exogenous variables
            end
            for j=1:p   
                boot_data(i,:,k) = boot_data(i,:,k)...
                    + boot_data(i-j,:,k)*A(:,:,j)';  % autoregressive parts          
            end
            boot_data(i,:,k) = boot_data(i,:,k)...  % stochastics
                + resid(index(i-p,k),:);    
        end
    end
end

% Close waitbar
close(h_wait)

% Put waitbar
h_wait = waitbar(0,'Bootstrapping impulse responses, please wait...');

% Generate bootstrap impulse responses
for k=1:rep
   waitbar(k/rep,h_wait)    
   
   % Switch over estimation method [OLS more accurate without restrictions]
   switch method
       
   case 'FGLS' % <------------- Prepare data for FGLS estimation
       
   % Left-hand side variables
   y_data = reshape(boot_data(1+p:end,:,k),(T-p)*n,1); 
   
   % Right-hand side variables
   X_data = [];
   for i=1:size(rhs_names,1)
       X_i = exo_data(1+p:end,:);  % exogenous terms
       q = 1;
       while q<=p      % go through lags
           for j=1:size(rhs_names,2)   % check rhs variables
               if (~isempty(char(rhs_names{i,j})))
                   pos_x = loc(char(lhs_names),char(rhs_names(i,j)));
                   X_i = [X_i boot_data(1+p-q:end-q,pos_x,k)];
               end
           end
           q = q+1;
       end
       X_data = blkdiag(X_data,X_i);
   end
   
   % Data to be used for current replication, SUR form
   this_data.y_sur = y_data;
   this_data.X_sur = X_data;
   
   case 'OLS' % <------------- Prepare data for OLS estimation
   
   % Left-hand side variables
   Y_data = boot_data(1+p:end,:,k); 
   
   % Right-hand side variables
   X_data = exo_data(1+p:end,:);  % exogenous terms
   q = 1;
   while q<=p      % go through lags
       X_data = [X_data boot_data(1+p-q:end-q,:,k)];
       q = q+1;
   end
   
   % Data to be used for current replication, OLS form
   this_data.Y_ols = Y_data;
   this_data.X_ols = X_data;
   
   end
      
   % Estimation by FGLS or OLS, is determined in the function
   results = EstimateVAR(this_data,info);
   
   % Calculate impulse responses
   boot_irf(:,:,k) = ComputeImpulseResponses(results,info);  
   
end

% Close waitbar
close(h_wait)

% Sort impulse responses
sort_irf = sort(boot_irf,3);      % sort at each t

% Compute lower and upper bands
for i=1:c
    bands(:,:,1,i) = sort_irf(:,:,ceil((1-(1+conf(i))/2)*rep));
    bands(:,:,2,i) = sort_irf(:,:,ceil((1+conf(i))/2*rep));
end

% Provide standard deviations of the bootstrapped impulse responses
if nargout==2
    sig=std(boot_irf,0,3);
end