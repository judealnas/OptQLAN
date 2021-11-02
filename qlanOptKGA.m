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

% Add a user
names = [names(:);'Dave'];

eta_N = [eta_N(:); 0.0001];
d_N = [d_N(:); 2000];

q = QLAN([], eta_N, d_N, chan_fluxes, U_K, rho_K, 5e-9, 'UserNames', names);

ys = q.coinc_wind_s.*q.d_N./q.eff_N;

q.configChans(chan_fluxes, U_K, rho_K);
%% Optimize with Varying K
allK = 1:25;
allQLANS = [];
num_opts = numel(allK);
soln_allocs = cell(numel(allK),1);
soln_totfluxes = zeros(numel(allK),1);
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
    
    % Reconfiguring QLAN object for current number of channels
    chan_fluxes = ones(k,1)./k; % Real channel fluxes assigned in fitness func
    U_K = repmat(eye(D), [1 1 k]);
    rho_K = repmat(rho_psi_m, [1 1 k]);
    q.configChans(chan_fluxes, U_K, rho_K);
    
    % GA settings
    num_vars = k+1; % optimize k allocations and total flux
    intcon = 1:k; % First k variables are integer allocations
    lb = [zeros(1,k) min(q.mu_L_ebrmax,[],'all')];
    ub = [q.L.*ones(1,k) 3*q.L*max(q.mu_L_ebrmax,[],'all')];
    func = @(x)gaFitFunc(x,q);
    nonlinconstr = [];%@(x)nonlinConstr(x,q);
    
    % Actual optimization
    [soln, fval(i)] = ...
        ga(func,num_vars,[],[],[],[],lb,ub,nonlinconstr, intcon, opts);
    soln_allocs{i} = soln(1:end-1);
    soln_totfluxes(i) = soln(end);
    
    % After optimization, put QLAN object in state with best found
    % allocation and save copy to an array
    gaFitFunc(soln,q);
    allQLANS = [allQLANS copy(q)];
end
%% Histogram of objective functions
figure()
histogram(fval)
title("Objective Value Histogram")

% For each solution, count occurences of each link index. Then histogram
% occurances of each unique link index
soln_hist_data = zeros(numel(soln_allocs), q.L+1);
for i = 0:q.L
    for j = 1:numel(soln_allocs)
        soln_alloc = soln_allocs{j};
        soln_hist_data(j,i+1) = sum(soln_alloc(:) == i);
    end
end

%% Plot of final objective values against K
figure()
semilogy(allK,fval)
xlabel("Number of Channels K")
ylabel("Objective Value")
title("Objective: Maximum Normalized EBR")
%% Plot of total flux against K
fig_totflux = figure();
ax_totflux = axes(fig_totflux);
plot(ax_totflux,allK,soln_totfluxes)
xlabel("# of Channels K")
title("Total Channel Flux")
%% Plot of # of allocated channels per link against K
fig_alloc = figure();
ax_alloc = axes(fig_alloc);
area(ax_alloc,allK,soln_hist_data)
grid(ax_alloc,'minor')
grid(ax_alloc,'on')
ax_alloc.Layer = 'top';
title("# of Allocated Channels")
xlabel("# of Channels K")
legend(["Null Link";compose("Link %s",cell2mat(q.link_labels))])

%% Comparing optimization solution to model maximums and generating GIF
fig = figure();
fig.Position = [800 200 1920/1.1 1080]./2;
ax_fig = axes(fig);

% Create subplots with theoretical link curves
L = q.L;
w = ceil(sqrt(L));
h = ceil(L/w);
ys = q.coinc_wind_s.*q.d_N(q.links)./q.eff_N(q.links);
domain = [0 2/q.coinc_wind_s];
subaxes = gobjects(L,1);
for i = 1:L
    curr_ax = subplot(h,w,i);
    subaxes(i) = curr_ax;
    
    model = QLANLinkModel(ys(i,:),q.D);
    fid = model.fidelityFunc();
    fid = @(mu)fid(q.coinc_wind_s.*mu);
    ebr = model.scaledEBRFunc();
    ebr = @(mu)prod(q.eff_N(q.links(i,:)))./q.coinc_wind_s.*ebr(q.coinc_wind_s.*mu);
    
    yyaxis(curr_ax,'left')
    fplot(curr_ax, fid, domain)
    ylabel("Fidelity")
%     yline(0.75,'--')
    ylim(curr_ax,[0 1]);
    
    yyaxis(curr_ax,'right')
    fplot(curr_ax,ebr,domain)
    ylabel("EBR")
    ylim(curr_ax,[0 q.ebr_Lmax(i)])
    
    title(sprintf("Link %s",q.link_labels{i}))
    xlabel(sprintf("\\mu_%d",i))
end
%%
% Add scatter points to theoretical curves
% scatter_axes = gobjects(L,2);
for j = 1:num_opts
    q = allQLANS(j);
    for i = 1:q.L
        curr_ax = subaxes(i);
        try
            delete(scatter_axes(i,:));
        catch
        end
        
        % Plotting EBR
        hold(curr_ax,'on')
        yyaxis(curr_ax, 'right')
        scatter_axes(i,1) = scatter(curr_ax,q.link_fluxes(i+1), q.ebr_L(i)); % Plot optim solution
        yyaxis(curr_ax,'left')
        scatter_axes(i,2) = scatter(curr_ax,q.link_fluxes(i+1), q.F_L(i)); % Plot solution point
        hold(curr_ax,'off')
    end
    
    sgtitle(sprintf("$K=%d$",allK(j)), 'interpreter', 'latex')
    drawnow
    
    % Capturing plot as gif frame
%     if j == 1
%         writeToGIF(fig, 'kopt.gif', true)
%     else
%         writeToGIF(fig,'kopt.gif', false, 'DelayTime', 1)
%     end
end
%% GA functions
function fitness = gaFitFunc(alloc, q)

% Re-configure QLAN object channel fluxes
tot_flux = alloc(end);
alloc = alloc(1:end-1);
K = q.K;
new_chans = tot_flux.*ones(K,1)./K;
q.configChans(new_chans, q.U_K, q.rho_K);

% Evaluate
q.evalAlloc(alloc, q.rho_K(:,:,1), QLAN.EVALMODE_ISO);

% Calcuate fitness
ebr_norm = q.ebr_L./q.ebr_Lmax;

fitness = sum(abs(1-ebr_norm));

end

function [c, ceq] = nonlinConstr(alloc,q)
% Enforces the constraint that each link flux be less than or equal to
% the link flux of maximum EBR. This is because both EBR and fidelity
% decrease monotonically for fluxes larger than the flux of maximum EBR
alloc_mat = QLAN.allocVect2Mat(alloc(1:end-1),q.L);
link_fluxes = QLAN.calcLinkFluxes(alloc_mat, q.chan_fluxes);
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