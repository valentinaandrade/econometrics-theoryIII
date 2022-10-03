function F = GenerateCompanionMatrix(A)

% GenerateCompanionForm
%
% Usage:
%   F = GenerateCompanionMatrix(A);
%
% Purpose:
%   Generates VAR(1) companion matrix of a VAR(p) from autoregressive lag 
%   matrices.
%
% Input:
%   A       n x n x p array, autoregressive lag matrices
% 
% Output:
%   F       p*n x p*n matrix, companion matrix
%
% Author:
%   Markus Kirchner, August 2012

% Get info
nlags = size(A,3);
nvars = size(A,1);

% Allocate memory
F = zeros(nvars,nvars*nlags);

% Add lag matrices
j = 1;
for i=1:nlags
    F(:,j:j+nvars-1) = A(:,:,i);
    j = j+nvars;
end

% Add identity and zeroes
F = [F; eye(nvars*(nlags-1)), zeros(nvars*(nlags-1),nvars)];