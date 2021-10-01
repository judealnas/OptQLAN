function nonadj_link_set = nonAdjLinkSet(links)
used_edges = [];
X = 1:size(links,1);
nonadj_link_set = {};
while ~isequal(sort(used_edges), X)
    links_avail = setdiff(X, used_edges);
    curr_link = links_avail(1);
    links_avail = links_avail(2:end);
    
    curr_users = links(curr_link,:);
    used_edges = union(used_edges, curr_link, 'stable');
    nonadj_links = curr_link;
    for i = 1:numel(links_avail)
        link_users = links(links_avail(i),:);
        if any(ismember(link_users, curr_users))
            continue
        else
            curr_users = union(curr_users, link_users, 'stable');
            used_edges = union(used_edges, [curr_link links_avail(i)], 'stable');
            nonadj_links = [nonadj_links ; links_avail(i)];
        end
    end
    
    nonadj_link_set = [nonadj_link_set ; nonadj_links];
end
end