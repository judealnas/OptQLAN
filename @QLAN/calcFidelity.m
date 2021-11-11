function F_L = calcFidelity(Q, psi_L, rho_L)
%CALCFIDELITY calculates the fidelities between pages of provided states
%
%   F_L = calcFidelity(Q, PSI, RHO) returns a vector of
%   fidelities between the corresponding matrices of the 3D
%   matrices RHO and PSI. RHO and PSI must have the same size
%   third dimension (i.e. same number of pages)
%       RHO - 3D matrix of L pages of matrices
%       PSI - 3D matrix of L pages of matrices
%       F_L - Lx1 vector of fidelities between corresponding
%             RHO's and PSI's
%
%   F_L = calcFidelity(Q, PSI) calculates the link fidelities with respect
%   to PSI using the noisy link states stored in the QLAN instance Q
narginchk(2,3);
if nargin == 2 && isa(Q,'QLAN')
    rho_L = Q.sigma_L;
end

if size(rho_L,3) ~= size(psi_L,3)
    error("3rd dimension of args must be equal")
end

L = size(rho_L,3);

F_L = zeros(L, 1);
for i = 1:L
    F_L(i) = fidelity(rho_L(:,:,i), psi_L(:,:,i));
end
end