function user_fluxes = calcUserFluxes(Q, incidence, link_fluxes)
%CALCUSERFLUXES calculates flux received by each user in network
%
%   user_fluxes = calcUserFluxes(Q, INCID, LINK_FLUX) calculates
%   the flux received by network users via matrix
%   multiplication of the incidene matrix and link flux vector.
%       INCID - NxL incidence matrix
%       LINK_FLUX - (L+1)x1 vector of L non-null-link fluxes
%       and 1 null link flux at index 1
%
%   user_fluxes = calcUserFluxes(Q) calculates the user fluxes using the
%   INCID and LINK_FLUX values stored within Q
%
% Inputs and outputs of this method are NOT stored within the instance Q

narginchk(1,3);
if nargin == 1 && isa(Q,'QLAN')
    incidence = Q.inc_mat;
    link_fluxes = Q.link_fluxes;
end
user_fluxes = incidence * link_fluxes(2:end);
end