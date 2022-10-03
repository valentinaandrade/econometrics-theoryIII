function A = GenerateLagMatrices(B)

% GenerateLagMatrices
%
% Usage:
%   A = GenerateLagMatrices(B);
%
% Purpose:
%   Generates matrices of a VAR lag polynominal from autoregressive 
%   coefficients in simultaneous equation form.
%
% Input:
%   B       n*p x n matrix, autoregr. coef. in simult. eq. form
% 
% Output:
%   A       n x n x p array, autoregressive lag matrices
%
% Author:
%   Markus Kirchner, August 2012

% Get info
nvars = size(B,2);
nlags = size(B,1) / nvars;

% Allocate memory
A = zeros(nvars,nvars,nlags);

% Transpose B for convenience
Bprime = B';

% Generate lag matrices
j = 1;
for i=1:nvars:1+nvars*(nlags-1)
    if j>nlags
        break
    end
    A(:,:,j) = Bprime(:,i:i+nvars-1);
    j = j+1;
end