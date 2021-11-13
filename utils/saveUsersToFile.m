function saveUsersToFile(filename, eff_N, d_N, names)
%SAVEUSERSTOFILE Save user information to text file in human-readable JSON
%format
%
%  saveUsersToFile(FILE, EFF, DARKRATE, NAMES) writes to FILE in JSON
%  format the users described by the EFF, DARKRATE, and optional NAMES
%  vectors.
%       EFF - Nx1 vector of detector efficiencies between 0 and 1
%       DARKRATE - Nx1 vector of detector dark count rates
%       NAMES - vector of char arrays or string vector of user names

narginchk(3,4);

% get number of users
N = min([numel(eff_N), numel(d_N)]);
if numel(eff_N) ~= numel(d_N)
    warning('Mismatched number of efficiencies and dark count rates. Ignoring some values')
end

% If names not provided, generate names
% Else check if enough names passed and pad if not
if nargin == 3
    names = arrayfun(@(x) sprintf('n%d', x), 1:N);
else
    names = cellstr(names);
    num_names = numel(names);
    
    % Add extra names if not enough passed
    if num_names < N
        warning('Insufficient number of names passed. Generating additional names...')
        names = [names(:);compose('n%d',((num_names+1):N).')];
    end
end
    

% construct user struct array
user_arr = [];
for i = 1:N
    user_arr = [user_arr
                struct('name', names(i),...
                       'eff', eff_N(i),...
                       'dark_count_rate', d_N(i))];
end

% Encode struct array to JSON and write to file
fp = fopen(filename, 'w');
json = jsonencode(user_arr, 'PrettyPrint', true);
fprintf(fp, "%s", json);
fclose(fp);

end


