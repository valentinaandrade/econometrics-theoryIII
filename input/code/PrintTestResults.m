function PrintTestResults(results)

% PrintTestResults(results)
%
% Usage:
%   PrintTestResults(results);
%
% Purpose:
%   Prints results of various tests conducted on a vector autoregression.
%
% Input:
%   results structure including different variables depending on test and
%   .test_type  string, type of test
%
% Output:
%   none
%
% Author:
%   Markus Kirchner, May 2012

% Get input
test_type = results.test_type;
nobs      = results.nobs;

% Switch over test types and print output
switch test_type

case 'granger_causality' % <-------------- Granger causality test

% Get remaining input  
F_stat    = results.F_stat;  
P_val     = results.P_val; 
lag_order = results.lag_order;
y_name    = results.y_name;
x_name    = results.x_name;

fid = 1;
fprintf(fid,'\n');
fprintf(fid,'********************************************************\n');
fprintf(fid,'Granger causality test (see Hamilton, 1994) \n');
fprintf(fid,'Included observations: %1.0f \n',nobs);
fprintf(fid,'Number of lags: %1.0f \n',lag_order);
fprintf(fid,'--------------------------------------------------------\n');
fprintf(fid,'Null hypothesis: %s',x_name);
fprintf(fid,' does not Granger-cause %s \n',y_name);
fprintf(fid,'--------------------------------------------------------\n');
fprintf(fid,'Test statistic (F)                            = %0.4f \n',F_stat(1));
fprintf(fid,'P-value                                       = %0.4f \n',P_val(1));
fprintf(fid,'--------------------------------------------------------\n');
fprintf(fid,'Null hypothesis: %s',y_name);
fprintf(fid,' does not Granger-cause %s \n',x_name);
fprintf(fid,'--------------------------------------------------------\n');
fprintf(fid,'Test statistic (F)                            = %0.4f \n',F_stat(2));
fprintf(fid,'P-value                                       = %0.4f \n',P_val(2));
fprintf(fid,'********************************************************\n');

case 'coefficient_restrictions' % <------- coefficient restrictions LR test

% Get remaining input
alpha     = results.alpha;
c_val     = results.c_val;    
lik_ratio = results.lik_ratio;  
P_val     = results.P_val; 
df        = results.df;

fid = 1;
fprintf(fid,'\n');
fprintf(fid,'****************************************************************\n');
fprintf(fid,'Likelihood ratio test of coef. restrictions (see Hamilton, 1994) \n');
fprintf(fid,'Included observations: %1.0f \n',nobs);
fprintf(fid,'----------------------------------------------------------------\n');
fprintf(fid,'Null hypothesis: coefficients jointly zero \n');
fprintf(fid,'Number of zero restrictions: %1.0f \n',df);
fprintf(fid,'----------------------------------------------------------------\n');
fprintf(fid,'Significance level (in percent)                       = %0.0f \n',alpha*100);
fprintf(fid,'Critical value Chi-square distribution                = %0.4f \n',c_val);
fprintf(fid,'Test statistic (likelihood ratio)                     = %0.4f \n',lik_ratio);
fprintf(fid,'P-value                                               = %0.4f \n',P_val);
fprintf(fid,'****************************************************************\n');

case 'lag_length' % <-------------------- lag length test

% Get remaining input
LogL      = results.LogL;
lik_ratio = results.lik_ratio;  
P_val     = results.P_val; 
aic       = results.aic;
sic       = results.sic;
hqc       = results.hqc;
lags      = results.lags;

% Prepare matrix to be printed
input.cnames = strvcat('LogL','LR','P-val','AIC','SIC','HQ');
input.rnames = strvcat('Lag',num2str(lags'));
fid = 1;
fprintf(fid,'\n');
fprintf(fid,'*********************************************************************\n');
fprintf(fid,'VAR lag order selection criteria (see Canova, 2007) \n');
fprintf(fid,'Included observations: %1.0f \n',nobs);
fprintf(fid,'---------------------------------------------------------------------\n');
mprint([LogL' lik_ratio' P_val' aic' sic' hqc'],input)
fprintf(fid,'---------------------------------------------------------------------\n');
fprintf(fid,'LogL  : log likelihood \n');
fprintf(fid,'LR    : sequential modified likelihood ratio test statistic \n');		
fprintf(fid,'P-val : P-value likelihood ratio test \n');
fprintf(fid,'AIC   : Akaike information criterion \n');		
fprintf(fid,'SC    : Schwarz information criterion \n');		
fprintf(fid,'HQ    : Hannan-Quinn information criterion \n');		
fprintf(fid,'*********************************************************************\n');
    
end