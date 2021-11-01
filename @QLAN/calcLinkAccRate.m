function A_L = calcLinkAccRate(Q, inc_mat, R_N, coinc_wind_s)
%CALCLINKACCRATE calculates the accidental
%coincidence rate of network links
%
%   A_L = calcLinkAccRate(Q, INCID, R_N, TAU) calculates the link
%   accidental rates A_L given:
%       INCID - NxL undirected incidence matrix
%       R_N - Nx1 vector of user single count rates
%       TAU - scalar double representing coincidence window in
%       seconds.
%       A_L - Lx1 vector of accidental coincidence rates for
%       each link
%
%   A_L = calcLinkAccRate(Q) calculates the link accidental rates A_L using
%   the values of INCID, R_N, and TAU stored within the QLAN object Q
%
%   Note that provided and calcualted values are NOT stored within the QLAN
%   object

narginchk(1,4);
if nargin == 1 && isa(Q,'QLAN')
    inc_mat = Q.inc_mat;
    R_N = Q.R_N;
    coinc_wind_s = Q.coinc_wind_s;
end
L = size(inc_mat,2);
A_L = zeros(L, 1);
full_inc_mat = logical(full(inc_mat));
for i = 1:L
    A_L(i) = ...
        8*coinc_wind_s.*prod(R_N(full_inc_mat(:,i)));
end
end