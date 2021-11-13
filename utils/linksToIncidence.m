function G = linksToIncidence(edges)
%LINKSTOINCIDENCE returns an incidence matrix given a matrix of edge nodes
%
%   G = linksToIncidence(EDGES) returns incidence matrix G where EDGES is
%   an Ex2 matrix where E is the number of edges and each row contains the
%   vertex indices of the edge end points. The size of G is
%   max(edges)x(size(edges,1))
E = size(edges,1);
V = max(edges, [], 'all');
G = zeros(V, E);

for e = 1:E
    G(edges(e,:), e) = 1;
end
end

