function link_fluxes = calcLinkFluxes(Q, alloc_mat, chan_fluxes)
%CALCLINKFLUXESreturns the flux per
%network link
%
%   link_fluxes = calcLinkFluxes(Q, ALLOC_MAT, CHAN_FLUX) calculates
%   the flux through the links as the matrix multiplication of the
%   allocation  matrix and the channel pair flux vector:
%       link_fluxes = ALLOC_MAT * CHAN_FLUX
%
%   link_fluxes = calcLinkFluxes(Q) calculates the link fluxes using the
%   ALLOC_MAT and CHAN_FLUX stored wihtin Q
%
% Note that inputs and outputs of this method are NOT stored within Q

if nargin == 1 && isa(Q,'QLAN')
    alloc_mat = Q.alloc_mat;
    chan_fluxes = Q.chan_fluxes;
end
link_fluxes = alloc_mat * chan_fluxes;
end