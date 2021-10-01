%% QLAN Optimization
try
    close('Genetic Algorithm')
catch 
end
clear
%% Network Parameters
P = sqrSincPowerSpect(192.3125e12, 310e9, 17*1e-8*0.7274);
chan_fluxes = [3;1;1;1;1;1;1;2]*1e6;
psi_m = [0;1;-1;0]/sqrt(2);
rho_psi_m = psi_m*psi_m';

D = size(rho_psi_m,1);
K = numel(chan_fluxes);
U_K = repmat(eye(D), [1 1 K]);
rho_K = repmat(rho_psi_m, [1 1 K]);

q = QLAN([], 'arxiv_users.json');
q.coinc_wind_s = 5e-9;
q.configChans(chan_fluxes, U_K, rho_K);
L = q.L;
ub = L*ones(K,1);
intcon = 1:K;
w = 0.5;
%% Optimize

opts = optimoptions('ga',...
    "CreationFcn", @boundIntegerCreationFcn,...
    "CrossoverFcn", @boundIntegerXOverFcn,...
    "MutationFcn", @boundIntegerMutationFcn,...
    "PlotFcn",{'gaplotbestf','gaplotbestindiv','gaplotscores'},...
    "Display",'iter',...
    "FunctionTolerance", 1e-12);

[solns, fval, exit, output, pop, scores] = ...
    ga(@(x)fitnessFunc(x,q,w),K,[],[],[],[],zeros(1,K),L*ones(1,K), [], opts);

q.evalAlloc(solns, q.rho_K(:,:,1))
q.F_L
q.ebr_L
%%
function pop = boundIntegerCreationFcn(GenomeLength, FitnessFunc, opts)
    total_pop_sz = sum(opts.PopulationSize);
    init_pop_sz = size(opts.InitialPopulation,1);
    num_new_genes = total_pop_sz - init_pop_sz;
    pop_bounds = opts.PopInitRange;
    
    pop = zeros(total_pop_sz, GenomeLength);
    % Add initial population first
    if init_pop_sz > 0
        pop(1:init_pop_sz,:) = opts.InitialPopulation;
    end    
    
    % Generate remaining genes
    new_genes = zeros(num_new_genes, GenomeLength);
    for i = 1:num_new_genes
        gene = zeros(1,GenomeLength);
        for j = 1:GenomeLength
            gene(j) = randi(pop_bounds(:,j));
        end
        new_genes(i,:) = gene;
    end
    
    pop((init_pop_sz+1):end,:) = new_genes;
end

function kids = boundIntegerXOverFcn(parents, options, nvars, FitnessFcn, unused, currPop)
%%XOVERFUNC: 
%   parents - vector of indices pointing to parents in currPop
% Two point crossover; From selected parents, create children from
% non-overlapping pairs of parents
    num_kids = length(parents)/2;
    kids = zeros(num_kids,nvars);
    parent_idx = 1;
    for i = 1:num_kids
        parent1 = currPop(parents(parent_idx),:);
        parent_idx = parent_idx + 1;
        parent2 = currPop(parents(parent_idx),:);
        parent_idx = parent_idx + 1;
        
        % Choosing two indices at which to splice parent genes
        split_idx = randi(nvars-1, [1 2]);
        while (split_idx(1) == split_idx(2))
            split_idx = randi(nvars-1, [1 2]);
        end
        
        % Handle case if split indices wrap
        if split_idx(1) < split_idx(2)
            left = split_idx(1);
            right = split_idx(2);
        else
            left = split_idx(2);
            right = split_idx(1);
            swap = parent1;
            parent1 = parent2;
            parent2 = swap;
        end
        
        kids(i,:) = [parent1(1:left), parent2((left+1):right), parent1((right+1):end)];
    end
end

function mutationChildren = boundIntegerMutationFcn(parents,...
    options, GenomeLength, fitnessfunc, state, currScore, currPop, mutationRate)

    if nargin < 8 || isempty(mutationRate)
        mutationRate = 0.01; % default mutation rate
    end
    
    lin_constrs = options.LinearConstr;
    
    mutationChildren = zeros(length(parents),GenomeLength);
    for i=1:length(parents)
        child = currPop(parents(i),:);
        % Each element of the genome has mutationRate chance of being
        % mutated.
        mutationPoints = find(rand(1,length(child)) < mutationRate);
        
        % Generate new genome within boundaries
        new_gene = zeros(1,length(mutationPoints));
        for j = 1:length(mutationPoints)
            new_gene(j) = randi([lin_constrs.lb(j) lin_constrs.ub(j)]);
        end
        child(mutationPoints) = new_gene;
        mutationChildren(i,:) = child;
    end
end

function fitness = fitnessFunc(alloc, q, w)
    alloc = round(alloc);
    alloc(alloc < 0) = 0;
    alloc(alloc > q.L) = q.L;
    
    q.evalIsoAlloc(alloc, q.rho_K(:,:,1));
    
    min_fid = mean(q.F_L);
    min_ebr = mean(q.ebr_L);
    
    % fitness is average fidelity and average ebit rate
    % Negative as algorithm must maximize values
    fitness = w/min_fid + (1-w)/min_ebr + q.link_fluxes(1)/sum(q.chan_fluxes);
end