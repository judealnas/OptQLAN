%% QLAN Optimization
try
    close('Genetic Algorithm')
    close all
catch
end
clear
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
ub = L*ones(K,1);
intcon = 1:K;
w = 0.5;

%% Optimize with Varying K
allK = 1:25;
allQLANS = [];
num_opts = numel(allK);
solns = cell(numel(allK),1);
fval = zeros(numel(allK),1);

opts = optimoptions('ga',...
    "PlotFcn",{},...{'gaplotbestf','gaplotbestindiv','gaplotscores','gaplotdistance', 'gaplotselection'},...
    "Display",'iter',...
    "FunctionTolerance", 1e-12,...
    "MaxStallGenerations", 100,...
    "OutputFcn", @gaoutputfunc,...
    "CrossoverFraction", 0.8);

% Generate new QLAN configuration and save to array to later re-evaluate
% optimization solutions
for i = 1:num_opts
    k = allK(i);
    fprintf("Using %d channels", k);
    
    % Ensure enough flux to reach optimal solution
    chan_fluxes = q.L*max(q.mu_L_ebrmax,[],'all').*ones(k,1)./k;

    % Double available flux to emulate resource surplus
    chan_fluxes = 2.*chan_fluxes;
    
    U_K = repmat(eye(D), [1 1 k]);
    rho_K = repmat(rho_psi_m, [1 1 k]);
    
    q.configChans(chan_fluxes, U_K, rho_K);
    
    allQLANS = [allQLANS copy(q)];
    
    [solns{i}, fval(i), exit, output, pop, scores] = ...
        ga(@(x)fitnessFunc(x,q),k,[],[],[],[],zeros(1,k),L*ones(1,k),...
           @(x)nonlinConstr(x,q), 1:k, opts);
end
%% Histogram objective functions and solutions
figure()
histogram(fval)
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

% Plot # of allocated channels per link
figure()
hold on
for i = 1:(q.L+1)
%     subplot(q.L+1,1,i)
    plot(allK,soln_hist_data(:,i),'DisplayName',sprintf("Link %d Allocated Channels", i-1))
end
title("No. of Allocated Channels")
xlabel("K")
ylabel("# of Channels")
legend
hold off

%% Comparing optimization solution to model maximums and generating GIF
symRE = QLANModel.ebrUniformState();
symF = QLANModel.fidUniformState();

h = figure();
h.Position(3:4) = [1920/1.5 1080]./2;

for j = 1:numel(solns)
    q = allQLANS(j);
    fitnessFunc(solns{j},q);
    clf(h)
    for i = 1:q.L
        subplot(q.L,1,i)
        sub_vals = [reshape(q.eff_N(q.links(i,:)),1,[]),...
            reshape(q.d_N(q.links(1,:)),1,[]),...
            q.coinc_wind_s, q.D];
        symRE_link = subs(symRE, QLANModel.getParamVars(), sub_vals);
        symF_link = subs(symF, QLANModel.getParamVars(), sub_vals);
        
        % Plotting EBR
        yyaxis right
        hold off
        fplot(symRE_link, [0 3e8]); % Plot model curve
        hold on
        ylabel('EBR')
        scatter(q.link_fluxes(i+1), q.ebr_L(i)) % Plot optim solution
        temp = ylim;
        ylim([0 q.ebr_Lmax(i)]) % Set lower ylim to 0
        
        yyaxis left
        ylabel('Fidelity')
        fplot(symF_link, [0 3e8]); % Plot model fidelity curve
        scatter(q.link_fluxes(i+1), q.F_L(i)) % Plot solution point
        ylim([0 1])
        
        title(sprintf("Link %s",q.link_labels{i}))
        xlabel(sprintf("\\mu_%d",i))
    
        % Add legend beside just the middle graph
        if true %i == ceil(q.L/2)
            legend({'Fidelity (Model)', 'Normalized EBR (model)', 'Fidelity (Opt)', 'Normalized EBR (Opt)'}, 'Location', 'bestoutside')
        end
    end
    
    sgtitle(sprintf("$K=%d$",allK(j)), 'interpreter', 'latex')
    drawnow
    
    % Capturing plot as image
    frame = getframe(h);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im,256);
    
    % Saving frame to GIF
    if j == 1
        imwrite(imind, cm, 'koptim.gif', 'gif', 'Loopcount',inf, 'DelayTime',1);
    else
        imwrite(imind, cm, 'koptim.gif', 'gif', 'WriteMode', 'append', 'DelayTime',0.5);
    end
    
end
%% GA functions
function fitness = fitnessFunc(alloc, q)
q.evalAlloc(alloc, q.rho_K(:,:,1),'EvalMode',QLAN.EVALMODE_ISO);

ebr_norm = q.ebr_L./q.ebr_Lmax;
fitness = sum(abs(ebr_norm - 1));
end

function [c, ceq] = nonlinConstr(alloc,q)
% Enforces the constraint that each link flux be less than or equal to
% the link flux of maximum EBR. This is because both EBR and fidelity
% decrease monotonically for fluxes larger than the flux of maximum EBR
alloc_mat = QLAN.allocVect2Mat(alloc,q.L);
link_fluxes = q.calcLinkFluxes(alloc_mat, q.chan_fluxes);
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