% Usage:
% 	ModeloBASE
% 
% Purpose:
%   This code estimates a VAR(p) allowing for exlusion restrictions on the 
%   autoregressive lag structure through iterated FGLS on SUR form (= ML).
%   The traditional OLS equation by equation (= ML) for the case without
%   exclusion restrictions is nested. The regression model in SUR form is
%       y_it = x_it' b_i + e_it     i = 1,...,n     t = 1,...,T 
%   To estimate the model OLS is used when x_it = x_jt for all i,j, 
%   otherwise FGLS is used when the RHS variables are different.
%
%   This code also calculates impulse responses to structural shocks by a 
%   Cholesky  decomposition of the estimated variance-covariance matrix of 
%   the reduced-form residuals, and computes the historical decomposition.

clear all;
close all;
clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fix the seed for the random number generators
% Note: results were generated with seed "20071984"
rng(20071984);

% Provide data file
% Format: [current path '\sub-directories\file']

info.file = [pwd '\data\data78.xlsx'];
% info.file = [pwd '\data\data85.xlsx'];

% Define left-hand side variables and their ordering 
% Format: the names need to be the same as in the data file

model = 4; 

switch model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           

    case 1

info.lhs_names  = {'rext_hp','pib','inf','tpm_hp','dtcr'};  % primera variable no responde de forma contemp. / segunda responde contemp. a la primera pero no al resto / última responde de forma contemp. a todas

% la siguiente matriz es simétrica si agrego una variable mas se debe
% agregar una fila adicional a la sgte matriz:

% Define right-hand side endogenous variables for each equation        
% Note: to impose zero restrictions use [] for the respective variable
info.rhs_names  = {'rext_hp','pib','inf','tpm_hp','dtcr';...
                   'rext_hp','pib','inf','tpm_hp','dtcr';...
                   'rext_hp','pib','inf','tpm_hp','dtcr';...
                   'rext_hp','pib','inf','tpm_hp','dtcr';...   
                   'rext_hp','pib','inf','tpm_hp','dtcr';...
                   };      

    case 2

info.lhs_names  = {'rext_hp','pib','inf_hp','tpm_hp','dtcr'};  % primera variable no responde de forma contemp. / segunda responde contemp. a la primera pero no al resto / última responde de forma contemp. a todas

% la siguiente matriz es simétrica si agrego una variable mas se debe
% agregar una fila adicional a la sgte matriz:

% Define right-hand side endogenous variables for each equation        
% Note: to impose zero restrictions use [] for the respective variable
info.rhs_names  = {'rext_hp','pib','inf_hp','tpm_hp','dtcr';...
                   'rext_hp','pib','inf_hp','tpm_hp','dtcr';...
                   'rext_hp','pib','inf_hp','tpm_hp','dtcr';...
                   'rext_hp','pib','inf_hp','tpm_hp','dtcr';...   
                   'rext_hp','pib','inf_hp','tpm_hp','dtcr';...
                   };                 
    case 3

info.lhs_names  = {'rext_hp','pib','defl','tpm_hp','dtcr'};  % primera variable no responde de forma contemp. / segunda responde contemp. a la primera pero no al resto / última responde de forma contemp. a todas

% la siguiente matriz es simétrica si agrego una variable mas se debe
% agregar una fila adicional a la sgte matriz:

% Define right-hand side endogenous variables for each equation        
% Note: to impose zero restrictions use [] for the respective variable
info.rhs_names  = {'rext_hp','pib','defl','tpm_hp','dtcr';...
                   'rext_hp','pib','defl','tpm_hp','dtcr';...
                   'rext_hp','pib','defl','tpm_hp','dtcr';...
                   'rext_hp','pib','defl','tpm_hp','dtcr';...   
                   'rext_hp','pib','defl','tpm_hp','dtcr';...
                   };      

    case 4

info.lhs_names  = {'rext_hp','pib','defl_hp','tpm_hp','dtcr'};  % primera variable no responde de forma contemp. / segunda responde contemp. a la primera pero no al resto / última responde de forma contemp. a todas

% la siguiente matriz es simétrica si agrego una variable mas se debe
% agregar una fila adicional a la sgte matriz:

% Define right-hand side endogenous variables for each equation        
% Note: to impose zero restrictions use [] for the respective variable
info.rhs_names  = {'rext_hp','pib','defl_hp','tpm_hp','dtcr';...
                   'rext_hp','pib','defl_hp','tpm_hp','dtcr';...
                   'rext_hp','pib','defl_hp','tpm_hp','dtcr';...
                   'rext_hp','pib','defl_hp','tpm_hp','dtcr';...   
                   'rext_hp','pib','defl_hp','tpm_hp','dtcr';...
                   }; 
               
    end; 
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     


% Define (optional) exogenous variables, the same for each equation   
info.exo_names  = {[]};

% Define the shocked variable
switch model
    case 1
info.shock_name = 'tpm_hp'; 
    case 2
info.shock_name = 'tpm_hp';
    case 3
info.shock_name = 'tpm_hp';
    case 4
info.shock_name = 'tpm_hp';
end;

% Set lag order p of the VAR; see "TestLagLength"
info.lag_order  = 3;   

% Set deterministic exogenous variables
% Note: set the do-variables on/off with 1/0
info.do_cons    = 0;                % add constant 
info.do_lin     = 0;                % add linear trend
info.do_quad    = 0;                % add linear-quadratic trend
info.do_cub     = 0;                % add linear-quadratic-cubic trend

% Add seasonal dummies according to seasonal frequency set here
% Note: 0 for no dummies, 4 for quarterly dummies, 12 for monthly
info.sea_freq   = 0;   

% Provide first and last observation of the estimation sample
% Format quarterly data: YEAR.QUARTER where 00=Q1, 25=Q2, 50=Q3 and 75=Q4
% Format monthly data: YEAR + MONTH/12 where MONTH = 0,...,11

info.first_obs  = 1978.00; 
% info.first_obs  = 1985.50; 
info.last_obs   = 2022.25;

% Settings to control impulse responses, bootstrap and probability bands
info.do_norm    = 1;                % 1 for normalized shocks (else 0)
info.norm_fac   = 0.01;             % normalization factor (e.g. 10 = 10%)
info.horizon    = 20;               % horizon of impulse responses
info.rep        = 10^3;             % number of bootstrap replications
info.conf       = [.68 .90 .95];    % significance levels for error bands

% Options for plots
info.widths     = [1 0.7];          % line widths (zero line, point estim.)
info.fsizes     = [8 8];            % font sizes (titles, axes)
info.area_color = [0.9 0.9 0.9];    % color of area for highest sign. level
info.area_step  = 0.1;              % reduction of color until lowest level

% Options for additional tests (e.g. lag length)
info.alpha      = 0.05;             % significance level
info.max_lag    = 8;                % maximum lag order to be tested

% Settings to control FGLS estimation when restrictions are imposed
info.conv_crit  = 1e-6;             % convergence criterion
info.max_it     = 100;              % maximum number of iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OLS or GLS regression

% Load data set
[raw_data,info] = LoadData(info);
% Prepare regression data
[data,info] = PrepareRegressionData(raw_data,info);
% Estimate the VAR
estim_results = EstimateVAR(data,info);
% Stability check
stable = CheckStability(estim_results.vecb,info);
% Warning message if coefficients unstable
if ~stable, disp('Coefficients unstable.'), end; clear stable;
% Plot regression results
PlotRegressionResults(data.all,estim_results,info);

%% Impulse responses

% Locate shocked variable
info.pos_shock = loc(char(info.names'),info.shock_name);
% Compute impulse responses
irf_point = ComputeImpulseResponses(estim_results,info);
% Compute confidence bands for impulse responses
[~,irf_sig] = ComputeConfidenceBands(estim_results,info,data);
% Plot impulse responses
PlotImpulseResponses(irf_point,irf_sig,info);

%% Additional tests

% Lag length tests
lag_test_results = TestLagLength(data,info);
% Print test results to screen
PrintTestResults(lag_test_results);
% Check whether coefficient restrictions are used
is_r = strcmp(estim_results.method,'FGLS');
% Tests of coefficient restrictions
if is_r, coef_test_results = TestCoefficientRestrictions(data,info); end
% Print test results to screen
if is_r, PrintTestResults(coef_test_results); end; clear is_r;

%% Historical decomposition

% Shock decomposition
estim_results = ComputeShockDecomposition(estim_results);
% Historical decomposition
estim_results = ComputeHistoricalDecomposition(estim_results,info,data);

%% Save decomposition to an xls-file

% xlswrite(info.file,...
% [info.dates,sum(estim_results.hist_dec(:,:,end:end-1),3),estim_results.hist_dec(:,:,1),...
% estim_results.hist_dec(:,:,2),estim_results.hist_dec(:,:,3)],'decomp.var','A3');

% Base Larga 1990Q1-2019Q1

%        A=zeros(117,6); %117 observations, 4 variables + constant = 6 cols.
%          for i=1:6
%              A(:,i)=estim_results.hist_dec(:,2,i); %save hist decomp of fourth variable in SVAR.
%         end; 
%         xlswrite('2_HD_BASE.xlsx',A,'GINI','B3');


% Save IRFs to Excel file

      load('low_band_90.mat')
      load('up_band_90.mat')
      load('low_band_68.mat')
      load('up_band_68.mat')

AA=[transpose(irf_point) transpose(low_90) transpose(up_90) transpose(low_68) transpose(up_68)];
BB=[find(lag_test_results.aic==min(lag_test_results.aic))-1; find(lag_test_results.sic==min(lag_test_results.sic))-1; find(lag_test_results.hqc==min(lag_test_results.hqc))-1];

switch model
    
    case 1
        xlswrite('IRFs_VARs.xlsx',AA,'output','C5');
        xlswrite('IRFs_VARs.xlsx',BB,'output','C28');  
    case 2
        xlswrite('IRFs_VARs.xlsx',AA,'output','AD5');
        xlswrite('IRFs_VARs.xlsx',BB,'output','AD28');   
    case 3
        xlswrite('IRFs_VARs.xlsx',AA,'output','BE5');
        xlswrite('IRFs_VARs.xlsx',BB,'output','BE28');   
    case 4
        xlswrite('IRFs_VARs.xlsx',AA,'output','CF5');
        xlswrite('IRFs_VARs.xlsx',BB,'output','CF28');           
end;

