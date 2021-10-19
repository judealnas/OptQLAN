clear
clc
%% Setup
chan_fluxes = [3;1;1;1;1;1;1;2]*1e6;
psi_m = [0;1;-1;0]/sqrt(2);
rho_psi_m = psi_m*psi_m';

D = size(rho_psi_m,1);
K = numel(chan_fluxes);
U_K = repmat(eye(D), [1 1 K]);
rho_K = repmat(rho_psi_m, [1 1 K]);

[names, eff, d] = loadUsersFromFile('optica_users.json');

q = QLAN([], eff, d, chan_fluxes, U_K, rho_K, 5e-9, 'UserNames', names);

alloc = [1;3;3;3;3;3;3;2];
%% Profiling
num_trials = 100;
totaltime = 0;
totalbatchtime = 0;
for i = 1:num_trials
profile on
q.evalAlloc(alloc,q.rho_K(:,:,1));
p_full = profile('info');

profile on
q.evalIsoAlloc(alloc, q.rho_K(:,:,1));
p = profile('info');

profile on
q.evalBatchIsoAlloc(alloc,q.rho_K(:,:,1));
p_batch = profile('info');

totaltime = totaltime + sum([p.FunctionTable(end).TotalTime]);
totalbatchtime = totalbatchtime + sum([p_batch.FunctionTable(end).TotalTime]);
end
avgtime = totaltime/num_trials
avgbatchtime = totalbatchtime/num_trials

profile_names = {'Full', 'Non-batch', 'Batch'};
profile_arr = [p_full p p_batch];
%% Build table comparing function execution times between profiles
% 1xP cell array of profile function names
fnames = arrayfun(@(s)cellstr(char(s.FunctionTable.FunctionName)), profile_arr, 'UniformOutput', false);
fnames_unique = unique(vertcat(fnames{:})); % list of all unique function names
exec_times = NaN(numel(fnames_unique),numel(profile_arr));

% iterate over all unique function names
for i = 1:numel(fnames_unique)
    % Current function to get total exec time
    fname = fnames_unique{i};
    
    % Iterate over available profiles
    for j = 1:numel(profile_arr)
        % Get functions in current profile
        prof_fnames = fnames{j};
        
        % find in current profile index of current function name 
        % if found, add execution time to table
        fidx = find(contains(prof_fnames, fname));
        if ~isempty(fidx)
            exec_times(i,j) = profile_arr(j).FunctionTable(fidx).TotalTime;
        end
    end
end

T = [table(fnames_unique,'VariableNames',{'FuncName'}),...
    array2table(exec_times, 'VariableNames', profile_names)];

tots = nanmax(exec_times);
T = [T;{'Total'}, mat2cell(tots, ones(size(tots,1),1) , ones(size(tots,2),1))]