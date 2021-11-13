function nonadj_edge_set = nonAdjEdgeSet(edges)
%NONADJedgeSET returns a cell array of vectors. Each cell contains row
%indices of non-adjacent edges, i.e. edges that do not share nodes.
%
%   SET = nonAdjedgeSet(edgeS) returns a cell array of vectors SET. Each
%   vector in SET contains row indices of the Nx2 matrix edgeS. The rows of
%   edgeS contains the node indices constituting an edge/edge.s
N = max(edges,[],'all'); % Get max number of nodes
nonadj_edge_set = {};
b_edges_used = false(size(edges,1),1); % True at edge index if edge has been used
while ~all(b_edges_used)
    % Get first available edge
    edges_avail = find(~b_edges_used);
    curr_edge = edges_avail(1); % Choose first of available edges
    
    % mark edge as used
    b_edges_used(curr_edge) = true;
    
    % Create list of currently used nodes
    % and add nodes of current edge
    b_node_used = false(N,1);
    curr_nodes = edges(curr_edge,:);
    b_node_used(curr_nodes) = true;
    
    % Initialize set vector to current edge 
    nonadj_edges = curr_edge;
    
    % Loop over remaining available edges
    edges_avail = edges_avail(2:end); % Exclude current edges
    for i = 1:numel(edges_avail)
        % Get nodes of current edge
        edge_nodes = edges(edges_avail(i),:);
        
        % If any nodes used, skip to next available link 
        if any(b_node_used(edge_nodes))
            continue
        else % Else mark nodes as used and add current link to set vector
            b_node_used(edge_nodes) = true;
            b_edges_used(edges_avail(i)) = true;
            nonadj_edges = [nonadj_edges ; edges_avail(i)];
        end
    end
    
    % Add vector to set
    nonadj_edge_set = [nonadj_edge_set ; nonadj_edges];
end
end