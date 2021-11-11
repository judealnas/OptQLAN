function R_N = calcUserSingleRate(Q, user_fluxes, eff_N, d_N)
%CALCUSERSINGLERATE calcualtes the single count rate at a user's detector
%
%   R_N = calcUserSingleRate(Q, USER_FLUX, EFF_N, D_N) returns
%   the single count rate given the user flux vector USER_FLUX,
%   user efficiency vector EFF_N, and user dark count rate
%   vector D_N. All arguments must have same number of
%   elements.
%
%   R_N - calcUserSingleRate(Q) calculates the user single rates using the
%   USER_FLUX, EFF_N, D_N values stored within the instance Q.
%
% Note that inputs and outputs of this method are NOT stored within the
% instance Q.
narginchk(1,4);
if nargin == 1 && isa(Q,'QLAN')
    user_fluxes = Q.user_fluxes;
    eff_N = Q.eff_N;
    d_N = Q.d_N;
end
R_N = eff_N.*user_fluxes./2 + d_N;
end