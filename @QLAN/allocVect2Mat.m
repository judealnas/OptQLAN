function M = allocVect2Mat(V,L)
%ALLOCVECT2MAT converts an allocation vector to an allocation
%matrix
%
%   M = allocVect2Mat(V,L) converts allocation vector V to
%   allocation matrix M given L non-null links. Each row of M
%   corresponds to a link. Each column of M corresponds to a
%   channel pair. If M(l,k) == 1, then the kth channel is
%   allocated to the lth link, where the 0th link is the null
%   link.

if ~(isvector(V) && all(mod(V,1) == 0, 'all'))
    error("V must be a vector of integers >= 0")
end
M = logical((0:L)' == V(:)');
end

