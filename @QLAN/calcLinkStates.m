function rho_L = calcLinkStates(Q, alloc_mat, link_fluxes,...
    chan_fluxes, U_K, rho_K)
%CALCLINKSTATES_ is a static method that calculates the state
%shared on a network link prior to the introduction of noise
%
%   rho_L = calcLinkStates(Q, ALLOC_MAT, LINK_FLUX, CHAN_FLUX,
%                         U_K, RHO_K)
%       ALLOC_MAT - LxK allocation matrix
%       LINK_FLUX - (L+1)x1 vector of L non-null-link fluxes
%                   and one null link flux at index 1
%       CHAN_FLUX - Kx1 vector of channel pair fluxes
%       U_K       - DxDxK 3D matrix of channel unitaries
%       RHO_K     - DxDxK 3D matrix of initial channel states
%   The allocation matrix is used to deterine which channels
%   are allocated to a link. The unitaries U_k of those
%   channels are applied to their respective intial channel
%   states rho_k. These states are summed, weighting each
%   channel state by the channel's flux contribution to the
%   link flux.
%
%   rho_L = calcLinkStates(Q) calculates the link states using the
%   ALLOC_MAT, LINK_FLUX, CHAN_FLUX, U_K, and RHO_K values stored within
%   the QLAN instance Q.
%
% Note that the inputs and outputs of this method are NOT stored within the
% QLAN instance Q.
narginchk(1,6);
if nargin == 1 && isa(Q,'QLAN')
    alloc_mat = Q.alloc_mat;
    link_fluxes = Q.link_fluxes;
    chan_fluxes = Q.chan_fluxes;
    U_K = Q.U_K;
    rho_K = Q.rho_K;
end
L = size(alloc_mat,1) - 1; % get number of non null links
rho_L = zeros([size(rho_K,[1 2]) L]);
for i = 2:(L+1) % iterate over non-null links
    if link_fluxes(i) ~= 0 % if link has no flux, leave rho_L = 0 matrix
        alloc_chans = logical(alloc_mat(i,:));
        U = U_K(:,:,alloc_chans); % unitaries of channels allocated to link i-1
        rho = rho_K(:,:,alloc_chans); % initial states of channels allocated to link i-1
        chan_flux = chan_fluxes(alloc_chans); % get channel fluxes allocated to link i-1
        chan_flux = reshape(chan_flux,1,1,numel(chan_flux)); % reshape to 3D vector for quicker math
        
        % Calculate sum of all chan_flux*U*rho*U' density matrices
        % associated with link and normalize by total link flux
        rho_L(:,:,i-1) = sum(pagemtimes(chan_flux,pagemtimes(pagemtimes(U,rho),'none',U,'ctranspose')),3);
        rho_L(:,:,i-1) = rho_L(:,:,i-1)./link_fluxes(i);
    end
end
end