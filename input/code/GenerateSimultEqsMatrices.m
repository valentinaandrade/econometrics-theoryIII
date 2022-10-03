function [B C] = GenerateSimultEqsMatrices(vecb,info)

% GenerateSimultEqsMatrices
%
% Usage:
%   [B C] = GenerateSimultEqsMatrices(vecb,info);
%
% Purpose:
%   Generates matrices of the simulataneous equations form of a VAR(p)
%   with exclusion restrictions from the vector of coefficients.
%
% Input:
%   vecb    n*m+n*p x 1 vector, coefficients in vector form
%   info    structure including
%   .nexo       integer, number of exogenous terms in the VAR
%   .lag_order  integer, VAR lag order p
%   .nvars      integer, number of variables n
%   .rhs_names  n x n matrix of strings, rhs variables per equation
% 
% Output:
%   B       n*p x n matrix, autoregr. coef. in simult. eq. form
%   C       m x n matrix, coefficients on exogenous variables
%
% Author:
%   Markus Kirchner, August 2012

% Get input
nexo      = info.nexo;
nlags     = info.lag_order;
nvars     = info.nvars;
rhs_names = info.rhs_names;

% Allocate memory
C = zeros(nexo,nvars);
B = zeros(nvars*nlags,nvars);

% Initialize counter
ib = 1;

% Go through equations
for ieq=1:nvars
    
    % Coefficients on exogenous variables
    C(:,ieq) = vecb(ib:ib+nexo-1);
    
    % Update counters
    ib = ib + nexo;
    iB = 1;
    
    % Autoregressive coefficients
    for ilag=1:nlags
        for jeq=1:nvars
            if ~isempty(char(rhs_names{ieq,jeq}))   % check restrictions
                B(iB,ieq) = vecb(ib);
                ib = ib + 1;    % next coefficient
            end
            iB = iB + 1;    % next row of B
        end
    end
    
end        