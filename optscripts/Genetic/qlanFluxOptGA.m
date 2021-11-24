%% QLAN Optimization
% Optimize QLAN with fixed channel distribution but variable total flux
try
    close('Genetic Algorithm')
    close all
catch
end
clear
files = genpath('../..')
files = 
addpath(genpath('../..'))
%% Network Parameters
psi_m = [0;1;-1;0]/sqrt(2);
rho_psi_m = psi_m*psi_m';
tau = 5e-9;

D = size(rho_psi_m,1);
K = 10;
chan_fluxes = rand(K,1); % random initialization of chan fluxes for construction of QLAN object
U_K = repmat(eye(D), [1 1 K]);
rho_K = repmat(rho_psi_m, [1 1 K]);

[names, eta_N, d_N] = loadUsersFromFile('arxiv_users.json');

q = QLAN([], eta_N, d_N, chan_fluxes, U_K, rho_K, 5e-9, 'UserNames', names);

q.configChans(chan_fluxes, U_K, rho_K);
L = q.L;

w = 0.5;

%% Optimize with varying total channel flux
% GA parameters
func = @(x) fitnessFunc(x,q); % optimization objective function
num_vars = K+1; % number of optimization variables
intcon = 1:K; % indices of integer optim variables
lb = [zeros(1,K) 1e6]; % lower bountds
ub = [L*ones(1,K) Inf]; % Upper bounds
nonlincon = @(x)nonlinConstr(x,q); % Non-linear constraint function
opts = optimoptions('ga',...
    "PlotFcn",{'gaplotbestf','gaplotbestindiv','gaplotscores','gaplotdistance', 'gaplotselection'},...
    "Display",'iter',...
    "FunctionTolerance", 1e-12,...
    "MaxStallGenerations", 100,...
    "OutputFcn", @gaoutputfunc,...
    "CrossoverFraction", 0.8);

num_opts = 10;
solns = zeros(num_opts,num_vars);
fvals = zeros(num_opts,1);

for i = 1:num_opts
    % Optimize to find optimal total channel flux and channel flux distribution
    [solns(i,:), fvals(i)] = ...
        ga(func,num_vars,[],[],[],[],lb,ub,nonlincon, intcon, opts);
end
%% Histogram objective functions and solutions
f_hist = figure();
[g, grps] = groupcounts(fvals);
stem(grps,g)
xticks(grps)
title("Objective Value Histogram")

% For each solution, count occurences of each link index. Then histogram
% occurances of each unique link index


%% Plotting results
f_plots = figure();
ax_plots = axes(f_plots);
domain = [0 3e8];
F = QLANModel.fidUniformState();
RE = QLANModel.ebrUniformState();

% Plotting theory curves for links
for i = 1:q.L
    ax_sub = subplot(q.L,1,i);
    ax_subs(i) = ax_sub;
    
    sub_vars = [QLANModel.sym_etaA
        QLANModel.sym_etaB
        QLANModel.sym_dA
        QLANModel.sym_dB
        QLANModel.sym_tau
        QLANModel.sym_D];
    
    sub_vals = [q.eff_N(q.links(i,:))
                q.d_N(q.links(i,:))
                q.coinc_wind_s
                q.D];
    
    hold(ax_sub, 'on')
    yyaxis(ax_sub ,'left')
    fplot(subs(F,sub_vars,sub_vals), domain)
    ylabel('Fidelity')
    ylim([0 1])
    yyaxis(ax_sub ,'right')
    fplot(subs(RE,sub_vars,sub_vals),domain)
    ylabel('EBR')
    ylim([0 q.ebr_Lmax(i)])
    xlabel(sprintf("\\mu_{%d}",i))
    hold(ax_sub, 'off')
end


% Plotting opt points
mkrs = 'o+*xsd^v<>ph';
num_mkrs = numel(mkrs);
for i = 1:size(solns,1)
    soln = solns(i,:);
    fitnessFunc(soln, q); % evaluate fitnessFunc to re-eval QLAN object
    mkr = mkrs(mod(i,num_mkrs));
    
    for j = 1:numel(ax_subs)
        ax = ax_subs(j);
        mu = q.link_fluxes(j+1);
        f = q.F_L(j);
        ebr = q.ebr_L(j);
        
        hold(ax,'on')
        yyaxis(ax,'left')
        scatter(ax, mu, f, mkr)
        yyaxis(ax,'right')
        scatter(ax, mu, ebr, mkr)
        hold(ax,'off')
    end
end
%% GA functions
function fitness = fitnessFunc(alloc, q)
tot_flux = alloc(end);
K = q.K;
chan_fluxes = tot_flux.*ones(K,1)./K;
q.configChans(chan_fluxes,q.U_K,q.rho_K);
q.evalAlloc(alloc(1:end-1), q.rho_K(:,:,1),QLAN.EVALMODE_ISO);

ebr_norm = q.ebr_L./q.ebr_Lmax;
fitness = sum(abs(ebr_norm - 1));
end

function [c, ceq] = nonlinConstr(alloc,q)
% Enforces the constraint that each link flux be less than or equal to
% the link flux of maximum EBR. This is because both EBR and fidelity
% decrease monotonically for fluxes larger than the flux of maximum EBR
alloc = alloc(1:end-1); % exclude total flux
alloc_mat = QLAN.allocVect2Mat(alloc,q.L);
link_fluxes = QLAN.calcLinkFluxes_(alloc_mat, q.chan_fluxes);
c = link_fluxes(2:end) - q.mu_L_ebrmax;
ceq = [];
end

function [state, options, optchanged] = gaoutputfunc(options,state,flag)
persistent pop_history
if isempty(pop_history)
    pop_history = {};
    pop_history(1) = {unique(state.Population,'stable', 'rows')};
else
    pop_history(end+1) = {unique(state.Population,'stable', 'rows')};
end
assignin('base','pop_history',pop_history)
optchanged = false;
end