function sigma_L = calcNoisyLinkStates(Q, noise_L, rho_L, lambda_L)
%CALCNOISYLINKSTATES calculates the link states after inclusion of noise.
%
%   sigma_L = calcNoisyLinkStates(Q, NOISE, RHO, LAMBDA) returns
%   the sum of the corresponding RHO's and NOISE's weighted by
%   the corresponding LAMBDA and 1-LAMBDA respectively
%       NOISE - DxDxL 3D matrix of noise states
%       RHO - DxDxL 3D matrix of link states
%       LAMBDA - Lx1 vector of weights
%       sigma_L - DxDxL 3D matrix of noisy link states
%           sigma_L(:,:,l) =
%               lambda_L(l)*rho_L(:,:,l) +
%               (1-lambda_L(l))*noise_L(:,:,l)
%
%   sigma_L = calcNoisyLinkStates(Q, NOISE) returns the link state using
%   the RHO and LAMBDA values stored within QLAN instance Q
narginchk(2,4);
if nargin == 2 && isa(Q,'QLAN')
    rho_L = Q.rho_L;
    lambda_L = Q.lambda_L;
end
lambda_L = reshape(lambda_L, 1,1, []);
sigma_L = pagemtimes(lambda_L, rho_L) + pagemtimes(1-lambda_L, noise_L);
end