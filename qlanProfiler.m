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

[names, eff, d] = loadUsersFromFile('arxiv_users.json');

q = QLAN([], eff, d, chan_fluxes, U_K, rho_K, 10e-9, 'UserNames', names);

alloc = [1;3;3;3;3;3;3;2;5e8];
%% Profiling
num_trials = 1000;
eval_modes = {QLAN.EVALMODE_ISO, QLAN.EVALMODE_SIMULT, QLAN.EVALMODE_BATCH};

for i = 1:numel(eval_modes)
    mode = eval_modes{i};
    profile on
    for j = 1:num_trials
        gaFitFunc(alloc, q, mode);
    end
    
    % Build vector of profiles containing total execution times
    ps(i) =  profile('info');
end
profile_names = eval_modes;
%% Build table comparing function execution times between profiles
% Get unique function names
all_profs = cell2mat(reshape({ps.FunctionTable},[],1));
fnames_unique = unique({all_profs.FunctionName});

% Matrix of average function execution times. Row index selects function
% and column index selects evaluation method.
exec_times = NaN(numel(fnames_unique),numel(ps)); 

% iterate over all unique function names
for i = 1:numel(fnames_unique)
    % Current function to get total exec time
    fname = fnames_unique{i};
    
    % Iterate over available profiles
    for j = 1:numel(ps)
        % Get all trials of evaluation method
        prof = ps(j);
        
        % Get functions in current profile
        prof_fnames = {prof.FunctionTable.FunctionName};
        
        % find in current profile index of current function name 
        % if found, add execution time to table
        fidx = find(strcmp(prof_fnames, fname));
        if ~isempty(fidx)
            exec_times(i,j) = max([prof.FunctionTable(fidx).TotalTime])/num_trials;
        end
    end
end

T = [table(fnames_unique(:),'VariableNames',{'FuncName'}),...
    array2table(exec_times, 'VariableNames', profile_names)];

tots = nanmax(exec_times);
T = [T;{'Total'}, mat2cell(tots, ones(size(tots,1),1) , ones(size(tots,2),1))]

%% Local Functions
function fitness = gaFitFunc(alloc, q, mode)

% Re-configure QLAN object channel fluxes
tot_flux = alloc(end);
alloc = alloc(1:end-1);
K = q.K;
new_chans = tot_flux.*ones(K,1)./K;
q.configChans(new_chans, q.U_K, q.rho_K);

% Evaluate
q.evalAlloc(alloc, q.rho_K(:,:,1), mode);

% Calcuate fitness
ebr_norm = q.ebr_L./q.ebr_Lmax;

% fitness = sum(abs(1-ebr_norm));
fitness = -sum(ebr_norm.*(q.F_L>=0.75));

end