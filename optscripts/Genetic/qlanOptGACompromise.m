%% qlanOptGACompromise.m
% This script uses the genetic algorithm to optimize a QLAN network
% consisting of 6 varying quality (different noise parameters)
% non-overlapping links. In addition to the integer allocation vector, the
% total channel flux is included as a continuous-valued optimization
% parameter. The channel flux distribution is uniform throughout the
% optimization.
%
% Optimizations are done for 4 different values of K (6, 12, 24, 48).
%
% Objective is normalized EBR = fidelity
%
% The objective value, # of channels allocated to each link, and the total
% channel flux are plotted against the # of channels K
try
    close('Genetic Algorithm')
    close all
catch
end
clear
%% Network Parameters
psi_m = [0;1;-1;0]/sqrt(2);
rho_psi_m = psi_m*psi_m';
tau = 10e-9;
fid_thresh  = 0.9;

D = size(rho_psi_m,1);
K = 10;
chan_fluxes = rand(K,1); % random initialization of chan fluxes for construction of QLAN object
U_K = repmat(eye(D), [1 1 K]);
rho_K = repmat(rho_psi_m, [1 1 K]);

% Calculating noise parameters for positive EBR roots
zh = linspace(2,1,6); % upper root of dimless EBR
y1_max = 1-zh./2;
y1s = 0.5.*y1_max;
y2s = (2*zh.*(1-y1s)-zh.^2)./(4*y1s+2*zh);

% pre-defined user parameters
user_ys = [y1s(:) y2s(:)];

names = cellstr(('A':char('A'+numel(user_ys)-1)).');
eta_N = 1e-4 + (1e-2-1e-4).*rand(size(user_ys));
eta_N = reshape(eta_N.',[],1);
d_N = reshape(user_ys.',[],1)./tau.*eta_N;
       
% Constructing one-to-one links
links = reshape(1:numel(user_ys),2,[]).';

q = QLAN(links, eta_N, d_N, chan_fluxes, U_K, rho_K, tau, 'UserNames', names);

q.configChans(chan_fluxes, U_K, rho_K);
%% Optimize with Varying K
allK = floor(q.L.*[1 2 4 8]);
allQLANS = [];
num_opts = numel(allK);
soln_allocs = cell(numel(allK),1);
soln_totfluxes = zeros(numel(allK),1);
fvals = zeros(numel(allK),1);

opts = optimoptions('ga'...
    ,"PlotFcn",{}...{'gaplotbestf'}...
    ...,"Display",'iter'...
    ,"FunctionTolerance", 1e-12...
    ,"MaxStallGenerations", 100 ...
    ...,"OutputFcn", @gaoutputfunc...
    ,"CrossoverFraction", 0.8);

% Generate new QLAN configuration and save to array to later re-evaluate
% optimization solutions
for i = 1:num_opts
    k = allK(i);
    fprintf("Using %d channels\n", k);
    
    % Reconfiguring QLAN object for current number of channels
    chan_fluxes = ones(k,1)./k; % Real channel fluxes assigned in fitness func
    U_K = repmat(eye(D), [1 1 k]);
    rho_K = repmat(rho_psi_m, [1 1 k]);
    q.configChans(chan_fluxes, U_K, rho_K);
    
    % GA settings
    num_vars = k+1; % optimize k allocations and total flux
    intcon = 1:k; % First k variables are integer allocations
    lb = [zeros(1,k) 0];
    ub = [q.L.*ones(1,k) 3*q.L*max(q.mu_L_ebrmax,[],'all')];
    func = @(x)gaFitFunc(x,q, fid_thresh);
    nonlinconstr = [];%@(x)nonlinConstr(x,q);
    
    % Run 5 optimizations, keeping best one
    last_fval = Inf;
    for n = 1:5
        [soln, fval] = ...
            ga(func,num_vars,[],[],[],[],lb,ub,nonlinconstr, intcon, opts);
        
        if fval < last_fval
            last_fval = fval;
            fvals(i) = fval;
            soln_allocs{i} = soln(1:end-1);
            soln_totfluxes(i) = soln(end);
        end
    end

    
    % After optimization, put QLAN object in state with best found
    % allocation and save copy to an array
    gaFitFunc(soln,q,fid_thresh);

    allQLANS = [allQLANS copy(q)];
end
%% Histogram of objective functions
figure(1)
histogram(fvals)
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
fig_fvals = figure(2);
ax_fvals = fig_fvals.CurrentAxes;
if isempty(ax_fvals)
    ax_fvals = axes(fig_fvals);
end
plot(ax_fvals, allK,fvals)
ax_fvals.XMinorGrid = true;
grid(ax_fvals,'on')
xlabel("Number of Channels K")  
ylabel("Objective Value")
title("Objective: Maximum Normalized EBR w/ F>="+num2str(fid_thresh))
%% Plot of total flux against K
fig_totflux = figure(3);
ax_totflux = fig_totflux.CurrentAxes;
if isempty(ax_totflux)
    ax_totflux = axes(fig_totflux);
end
plot(ax_totflux,allK,soln_totfluxes)
grid(ax_totflux,'on');
ax_totflux.XMinorGrid = true;
xlabel("# of Channels K")
% ylim([lb(end),ub(end)])
title("Total Channel Flux")
%% Plot of # of allocated channels per link against K
fig_alloc = figure(4);
ax_alloc = fig_alloc.CurrentAxes;
if isempty(ax_alloc)
    ax_alloc = axes(fig_alloc);
end
area(ax_alloc,allK,soln_hist_data)
grid(ax_alloc,'minor')
grid(ax_alloc,'on')
ax_alloc.Layer = 'top';
title("# of Allocated Channels")
xlabel("# of Channels K")
legend(["Null Link";compose("Link %s",cell2mat(q.link_labels))]...
    ,'Location','northwest')

%% Comparing optimization solution to model maximums
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
fig = figure(5);
ax_fig = fig.CurrentAxes;
if isempty(ax_fig)
    ax_fig = axes(fig);
end
colors = colororder; % get default MATLAB color order
fig.Position = [800 200 1920/1.1 1080]./2;


% Create subplots with theoretical link curves and optimization solutions
L = q.L;
w = ceil(sqrt(L));
h = ceil(L/w);
ys = q.coinc_wind_s.*q.d_N(q.links)./q.eff_N(q.links);

% Collated results of K sweeps; Column index = K sweep
norm_ebrs = [allQLANS.ebr_L]./[allQLANS.ebr_Lmax];
fids = [allQLANS.F_L];
fluxes = [allQLANS.link_fluxes];

% X-axis range
domain = [0 max(fluxes,[],'all')*1.15];
mrkrs = {'o','^','s','d','p','h'};
for i = 1:L
    plot_ind = i + floor(i/(w+1));
    curr_ax = subplot(h,w+1,plot_ind);
    
    % Get model function handles
    model = QLANLinkModel(ys(i,:),q.D);
    fid = model.fidelityFunc();
    fid = @(mu)fid(q.coinc_wind_s.*mu);
    dimless_ebr = model.dimlessEBRFunc();
    ebr = @(mu)prod(q.eff_N(q.links(i,:)))./q.coinc_wind_s.*dimless_ebr(q.coinc_wind_s.*mu);
    norm_ebr = @(mu)ebr(mu)./q.ebr_Lmax(i);
    
     % Plot fidelity curve
    fplot(curr_ax, fid, domain,'DisplayName','Fidelity (Model)');
    
    % Threshold line plot; exclude from legend
    yl = yline(fid_thresh,'--'); 
    yl.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    hold(curr_ax,'on')
    
    % plot normalized EBR
    fplot(curr_ax,norm_ebr,domain,'DisplayName','$R/R_{max}$ (Model)')
    
    % Construct scatter plot data; for each link, row idx = K val idx, col
    % idx = fidelity or EBR
    x = fluxes(i+1,:);
    y = [fids(i,:).' norm_ebrs(i,:).'];
    
    % Group data by number of channels K, represented by point shape
    % Do not display in legend (manually construct later)
    gls = gscatter(x,y,allK,repmat('k',numel(allK),1),[mrkrs{1:numel(allK)}],6,'off'); %gls = GxN matrix of line objects; G=groups, N=points in group
    gls_ann = [gls.Annotation];
    gls_legend = [gls_ann.LegendInformation];
    [gls_legend.IconDisplayStyle] = deal('off');
    
    % iterate over groups
    for k = 1:size(gls,1)
        gr = gls(k,:);
        
        % Add empty to data to add custom symbols/colors to legend.
        % Add entry for each scatter plot shape
        plot(nan,nan,[mrkrs{k}, 'k'],'DisplayName',sprintf("$K=%d$ (GA)",allK(k)))
        
        % iterate over data members of group, changing marker shape, color,
        % and size
        for j = 1:numel(gr)
            g = gr(j);
            g.Marker = mrkrs{k};
            g.MarkerSize = 7;
            g.Color = colors(j,:);
        end
    end
    
    
    title(sprintf("Link %s",q.link_labels{i}))
    ylim(curr_ax,[0 1]);
    xlabel(sprintf("\\mu_{%s}",q.link_labels{i}))
    ylabel('')
    
    
    hold(curr_ax,'off')
end
legend('Interpreter','latex')
%% Plotting optimization results using bar graph
link_cats = categorical(q.link_labels);
%% GA Functions
function fitness = gaFitFunc(alloc, q, fid_thresh)
% Re-configure QLAN object channel fluxes
tot_flux = alloc(end);
alloc = alloc(1:end-1);
K = q.K;
new_chans = tot_flux.*ones(K,1)./K;
q.configChans(new_chans, q.U_K, q.rho_K);

% Evaluate
q.evalAlloc(alloc, q.rho_K(:,:,1), QLAN.EVALMODE_SIMULT);

% Calcuate fitness
ebr_norm = q.ebr_L./q.ebr_Lmax;
if any(q.link_fluxes(2:end) > q.mu_L_ebrmax)
    fitness = Inf;
else
    fitness = max(abs(ebr_norm - q.F_L)./ebr_norm);
end

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