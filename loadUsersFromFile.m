function [user_names, eff_N, d_N] = loadUsersFromFile(filename)
%[user_names, eff_N, d_N] = loadUsersFromFile(filename)
%   returns corresponding arrays of user parameters
%   contained in filename
%
%   user_names = string array of user names
%   eff_N = array of user detector efficiencies
%   d_N = array of user detector dark count rates.

% Decode struct array from file
user_arr = jsondecode(fileread(filename));

% Initializing arrays
N = numel(user_arr);
eff_N = zeros(N,1);
d_N = zeros(N,1);
user_names = strings(N,1);

% Populating arrays
for i = 1:N
    eff_N(i) = user_arr(i).eff;
    d_N(i) = user_arr(i).dark_count_rate;
    user_names(i) = user_arr(i).name;
end

