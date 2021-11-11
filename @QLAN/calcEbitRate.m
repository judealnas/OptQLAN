function ebr_L = calcEbitRate(Q, A_L, C_L, sigma_L)
%CALCEBITRATE Calculates the entangled bit rate over network links
%
%   ebr_L = calcEbitRate(Q, A_L, C_L, SIGMA_L) calculates the
%   link entangled bit rate using the link coincidence rates
%   A_L, C_L and the log negativity of the noisy link states
%   SIGMA_L
%       A_L - Lx1 vector of link accidental coincidence rates
%       C_L - Lx1 vector of link truue coincidence rates
%       SIGMA_L = DxDxL 3D matrix of noisy link states used for
%       log negativity calculation
%
%   ebr_L = calcEbitRate(Q) calcualtes the link entangled bit rates using
%   the A_L, C_L, and SIGMA_L stored in QLAN instance Q
narginchk(1,4);
if nargin == 1 && isa(Q,'QLAN')
    A_L = Q.A_L;
    C_L = Q.C_L;
    sigma_L = Q.sigma_L;
end
L = size(sigma_L,3);
ebr_L = zeros(L,1);
for i = 1:L
    ebr_L(i) = logneg(sigma_L(:,:,i));
end

ebr_L = (A_L + C_L).*ebr_L;
end