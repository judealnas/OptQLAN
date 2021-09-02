function saveUsersToFile(filename, eff_N, d_N, names)
%SAVEUSERSTOFILE Save user information to text file in human-readable JSON
%format
%
%  filename = File to write; pre-existing contents will be overwritten
%  eff_N = array of user detector efficiencies
%  d_N = array of user dark count rates
%  names = [OPTIONAL] string array of user names; if not provided, default
%          names of style "nx" with integer x is used

% get number of users
N = min([numel(eff_N), numel(d_N)]);

% If names not provided, generate names
if nargin == 3
    names = arrayfun(@(x) sprintf("n%d", x), 1:N);
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


