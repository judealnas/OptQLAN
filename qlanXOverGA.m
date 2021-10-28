%% qlanXOverGA.m
% Script to benchmark genetic algorithm optimization of QLAN network with
% varying crossover fraction.
%
% Crossover fraction is the fraction of non-elite genes are used in
% crossover gene generation. (1 - (crossover fraction)) fraction of genes
% are then used in mutatation generation.
%
% To best assess influence of crossover fraction, a QLAN with a very
% specific solution must be used. Here, a QLAN with large K and K = L is
% optimized to maximize entangled bit rate. The ideal solution is one
% channel assigned to each link, i.e. an allocation vector of 1:L. K is
% large to increase difficulty in finding this solution.

%% QLAN object instantiation
N = 7;
links = nchoosek(1:7,2);
K = size(links,1);

psi_m = [0;1;-1;0]/sqrt(2);
rho_psi_m = psi_m*psi_m';
tau = 5e-9;

D = size(rho_psi_m,1);
U_K = repmat(eye(D), [1 1 K]);
rho_K = repmat(rho_psi_m, [1 1 K]);

eta_N = 1e-2*ones(N,1);
d_N = 50.*ones(N,1);

% Instantiate QLAN object with dummy channel flux values
q = QLAN([], eta_N, d_N, ones(K,1), U_K, rho_K, 5e-9);

% Ensure enough flux to reach optimal solution
chan_fluxes = q.L*max(q.mu_L_ebrmax,[],'all').*finiteUnitGeoSum(K,0.925);

q.configChans(chan_fluxes, U_K, rho_K);
%% Sweeping Xover fraction
xover_fracs = linspace(0,1,10);
solns = zeros(numel(xover_fracs),K);
fvals = zeros(numel(xover_fracs),1);

% GA parameters
func = @(x) fitnessFunc(x,q); % optimization objective function
num_vars = K; % number of optimization variables
intcon = 1:K; % indices of integer optim variables
lb = [zeros(1,K)]; % lower bountds
ub = [q.L*ones(1,K)]; % Upper bounds
nonlincon = @(x)nonlinConstr(x,q); % Non-linear constraint function
opts = optimoptions('ga',...
    "PlotFcn",{'gaplotbestf','gaplotbestindiv','gaplotscores','gaplotdistance'},...
    "Display",'iter',...
    "FunctionTolerance", 1e-12,...
    "MaxStallGenerations", 100);

for i = 1:numel(xover_fracs)
    opts.CrossoverFraction = xover_fracs(i);
    [solns(i,:), fvals(i)] = ga(func, num_vars,[],[],[],[],lb,ub,nonlincon,intcon, opts);
end

%% GA functions
function fitness = fitnessFunc(alloc, q)
    q.evalIsoAlloc(alloc, q.rho_K(:,:,1));
    ebr_norm = q.ebr_L./q.ebr_Lmax;
    fitness = sum(abs((ebr_norm - 1)));
end

function [c, ceq] = nonlinConstr(alloc,q)
    % Enforces the constraint that each link flux be less than or equal to
    % the link flux of maximum EBR. This is because both EBR and fidelity
    % decrease monotonically for fluxes larger than the flux of maximum EBR
    alloc_mat = QLAN.allocVect2Mat(alloc,q.L);
    link_fluxes = QLAN.calcLinkFluxes_(alloc_mat, q.chan_fluxes);
    c = link_fluxes(2:end) - q.mu_L_ebrmax;
    ceq = [];
end

