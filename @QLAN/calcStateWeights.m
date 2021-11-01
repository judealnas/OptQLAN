function lambda_L = calcStateWeights(Q, A_L, C_L)
%CALCSTATEWEIGHTS calculates the relative weights between the maximally
%mixed state and true desired state per link
%
%   lambda_L = calcStateWeights(Q, A_L,C_L) calculates the weight
%   vector as
%       lambda_L = 1 ./ (1 + A_L ./ C_L)
%   where A_L and C_L are the link accidental and true
%   coincidence rates respectively
%
%   lambda_L = calcStateWeights(Q) calcualtes the relative state weights
%   using the A_L and C_L stored in the QLAN instance Q
%
% Note that inputs and outputs of this method are NOT stored within the
% instance Q
narginchk(1,3);
if nargin == 1 && isa(Q,'QLAN')
    A_L = Q.A_L;
    C_L = Q.C_L;
end
lambda_L = 1./(1 + A_L ./ C_L);
lambda_L(isnan(lambda_L)) = 0;
end