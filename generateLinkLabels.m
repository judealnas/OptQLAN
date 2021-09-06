function link_labels = generateLinkLabels(user_names, links)
%GENERATELINKLABELS Generate list of suggested link labels
%   labels = generateLinkLabels(USERS, LINKS) returns a suggested list l of link labels
%   from the list of user names USERS and network links LINKS. 
%       USERS = a string or cell array of user names
%       LINKS = a numel(N)x2 numeric array; each row contains the indices of
%       the two users constituting the link
%       labels = a size(L,1)x1 cell array of character vectors; suggested link
%       labels

user_names = cellstr(user_names); % ensure names in cell array form
N = numel(user_names); % get # of users
L = size(links,1); % get # of links

% Generate array of unique characters
unique_chars = cell(N,1);
for i = 1:N
    name = user_names{i};
    rem_names = user_names(setdiff(1:N,i));
    rem_chars = unique(cat(2,rem_names{:}));
    unique_chars{i} = setdiff(name, rem_chars,'stable');
    unique_chars{i} = unique_chars{i}(1);
end

% Using provided links, generate combination of unique chars
link_labels = cell(L,1);
for i = 1:L
    link_labels{i} = [unique_chars{links(i,1)} unique_chars{links(i,2)}];
end
