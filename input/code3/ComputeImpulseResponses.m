    function irf = ComputeImpulseResponses(input,info)

% ComputeImpulseResponses
%
% Usage:
%   irf = ComputeImpulseResponses(input,info);
%
% Purpose:
%   Computes structural impulse responses for a reduced-form VAR that 
%   allows for general exclusion restrictions (SUR form). Identification 
%   is made by a Choleski decomposition of the estimated reduced-form 
%   residual variance-covariance matrix. Returns an array of dimension 
%   (i,j,k) with the response of variable i at horizon j due to an initial 
%   shock to variable k.
%
% Input:
%   input   structure including
%   .B          max(k_i)*p x n matrix, autoregr. coef. in simult. eqs. form
%   .Omega      n x n matrix, residual variance-covariance matrix
%   info    structure including
%   .lag_order  integer, VAR lag order p
%   .horizon    integer, horizon h for impulse response functions
%   .do_norm    boolean, 1 for normalized impulse responses, else 0
%   .norm_fac   double, normalization factor in percent
%   .pos_shock  integer, position of shock of interest
%
% Output:
%   irf     n x h+1 array, structural impulse responses to chosen shock
%
% Author:
%   Markus Kirchner, May 2012

% Provide inputs
B         = input.B;
Omega     = input.Omega;
n         = size(Omega,1);
p         = info.lag_order;
h         = info.horizon;
do_norm   = info.do_norm;
norm_fac  = info.norm_fac;
pos_shock = info.pos_shock;

% Allocate memory
red_irf   = zeros(n,n,h+1);
chol_irf  = zeros(n,n,h+1);
order_irf = zeros(n,h+1,n);
lag_mat   = zeros(n,n,p); 

% Generate lag matrices
Bprime = B';
j      = 1;
for i=1:n:1+n*(p-1)
    if j>p
        break
    end
    lag_mat(:,:,j) = Bprime(:,i:i+n-1);
    j = j+1;
end

% Reduced-form impact responses
red_irf(:,:,1) = eye(n);

% Compute remaining reduced-form responses
for i=2:h+1
    temp = zeros(n);
    for j=1:p
        if j>=i
            break
        end
        temp = temp + red_irf(:,:,i-j) * lag_mat(:,:,j);
    end
    red_irf(:,:,i) = temp;
end

% Orthogonalization
cho_fac = chol(Omega)';
for i=1:h+1
   chol_irf(:,:,i) = red_irf(:,:,i) * cho_fac;  
end

% Re-order array dimensions
for i=1:h+1
    for j=1:n  
        order_irf(:,i,j) = chol_irf(:,j,i);  
    end
end

% Prepare output
irf = squeeze(order_irf(:,:,pos_shock));

% Normalization
if do_norm
    irf = norm_fac*irf/irf(pos_shock,1);
end