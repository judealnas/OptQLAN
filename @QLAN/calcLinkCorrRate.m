function C_L = calcLinkCorrRate(Q, inc_mat, link_fluxes, eff_N)
%CALCLINKCORRRATE calculates the true,
%correlated coincidence rate of each link in a network.
%
%   C_L = calcLinkCorrRate(Q, INCID, LINK_FLUX, EFF_N) returns
%   the link true coincidence rate vector C_L.
%       INCID - NxL undirected incidence matrix
%       LINK_FLUX - (L+1)x1 vector of L non-null-link fluxes
%                   and one null link flux at index 1
%       EFF_N - Nx1 vector of user detector efficiencies
%       C_L - Lx1 vector of true coincidence rates for each
%               link
%
%   C_L = calcLinkCorrRate(Q) calculates the true coincidence rates using
%   the INCID, LINK_FLUX, and EFF_N values stored in the QLAN instance Q
%
% Note that inputs and outputs of this methods are NOT stored within the
% instance Q
narginchk(1,4);
if nargin == 1 && isa(Q,'QLAN')
    inc_mat = Q.inc_mat;
    eff_N = Q.eff_N;
    link_fluxes = Q.link_fluxes;
end

L = size(inc_mat,2);
C_L = zeros(L, 1);
full_inc_mat = logical(full(inc_mat));
for i = 1:L
    C_L(i) = ...
        link_fluxes(i+1)*prod(eff_N(full_inc_mat(:,i)));
end
end