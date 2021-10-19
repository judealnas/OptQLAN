%% QLAN Optimization
try
    close('Genetic Algorithm')
    close all
catch 
end
clear
%% Network Parameters
K = 10;
psi_m = [0;1;-1;0]/sqrt(2);
rho_psi_m = psi_m*psi_m';
tau = 5e-9;

D = size(rho_psi_m,1);
U_K = repmat(eye(D), [1 1 K]);
rho_K = repmat(rho_psi_m, [1 1 K]);

[names, eta_N, d_N] = loadUsersFromFile('arxiv_users.json');

% Instantiate QLAN object with dummy channel flux values
q = QLAN([], eta_N, d_N, ones(K,1), U_K, rho_K, 5e-9, 'UserNames', names);

% Ensure enough flux to reach optimal solution
chan_fluxes = q.L*max(q.mu_L_ebrmax,[],'all').*ones(K,1)./K;

% Halve the available flux to emulate resource starvation
chan_fluxes = chan_fluxes./2;

% Double available flux to emulate resource surplus
% chan_fluxes = chan_fluxes.*2;

q.configChans(chan_fluxes, U_K, rho_K);
L = q.L;
ub = L*ones(K,1);
intcon = 1:K;
w = 0.5;
%% Optimize Several Times
num_opts = 100;
solns = cell(num_opts,1);
fval = zeros(num_opts,1);

opts = optimoptions('ga',...
    "PlotFcn",{},...{'gaplotbestf','gaplotbestindiv','gaplotscores','gaplotdistance'},...
    "Display",'iter',...
    "FunctionTolerance", 1e-12,...
    "MaxStallGenerations", 100,...
    "CrossoverFraction", 0.8);

for i = 1:num_opts
    [solns{i}, fval(i), exit, output, pop, scores] = ...
        ga(@(x)fitnessFunc(x,q,w),K,[],[],[],[],zeros(1,K),L*ones(1,K), @(x)nonlinConstr(x,q), intcon, opts);
end

mean_fval = mean(fval);
% histogram(fval);
%% Histogram objective functions and solutions
figure()
[counts, groups] = groupcounts(fval);
stem(groups, counts,'.')
title("Objective Value Histogram")

% For each solution, count occurences of each link index. Then histogram
% occurances of each unique link index
soln_hist_data = zeros(numel(solns), q.L+1);
for i = 0:q.L
    for j = 1:numel(solns)
        soln = solns{j};
        soln_hist_data(j,i+1) = sum(soln(:) == i);
    end
end

% Histogram link indices
figure()
for i = 1:(q.L+1)
    subplot(q.L+1,1,i)
    [counts, groups] = groupcounts(soln_hist_data(:,i));
    stem(groups, counts,'.','MarkerSize',10)
    if i == 1
        title("Null Link Channel Count Histogram");
    else
        title(sprintf("Link %s Channel Count Histogram", q.link_labels{i-1}));    
    end
end

%% Comparing optimization solution to model maximums
symRE = QLANModel.ebrUniformState();
symF = QLANModel.fidUniformState();
[~,best_soln_idx] = min(fval,[],'all','linear'); 
best_soln = solns{best_soln_idx};

fitnessFunc(best_soln, q, w)
figure()
for i = 1:q.L
    subplot(q.L,1,i)
    hold on
    sub_vals = [reshape(q.eff_N(q.links(i,:)),1,[]),...
        reshape(q.d_N(q.links(1,:)),1,[]),...
        q.coinc_wind_s, q.D];
    symRE_link = subs(symRE, QLANModel.getParamVars(), sub_vals);
    symF_link = subs(symF, QLANModel.getParamVars(), sub_vals);
    
    % Plotting EBR
    yyaxis right
    fplot(symRE_link, [0 3e8]); % Plot model curve
    ylabel('EBR')
    scatter(q.link_fluxes(i+1), q.ebr_L(i)) % Plot optim solution
    temp = ylim;
    ylim([0 q.ebr_Lmax(i)]) % Set lower ylim to 0
    
    yyaxis left
    ylabel('Fidelity')
    fplot(symF_link, [0 3e8]); % Plot model fidelity curve
    scatter(q.link_fluxes(i+1), q.F_L(i)) % Plot solution point
    ylim([0 1])
    
    xlabel(sprintf("\\mu_%d",i))
    title(sprintf("Link %s",q.link_labels{i}))
    legend({'Fidelity (Model)', 'EBR (Model)', 'Fidelity (Opt)', 'EBR (Opt)'}, 'Location', 'bestoutside')
    hold off
end

%% GA functions
function fitness = fitnessFunc(alloc, q, w)
    q.evalIsoAlloc(alloc, q.rho_K(:,:,1));
    
    fid = q.F_L;
    ebr_norm = q.ebr_L./q.ebr_Lmax;
    f = q.ebr_L;
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