clear
close all
%%
[names, eff_N, d_N] = loadUsersFromFile('arxiv_users.json');
alloc = [1;2];
K = numel(alloc);
D = 4;
tau = 5e-9;

psi_m = [0;1;-1;0]/sqrt(2);
rho_K = repmat(psi_m*psi_m',[1 1 K]);
U_K = repmat(eye(D),[1 1 K]);
q = QLAN([1 2; 2 3], eff_N, d_N, 1:K, U_K, rho_K, tau,...
    'UserNames', names);

[flux1, flux2] = meshgrid(linspace(0,4e8,200));

fids = zeros([size(flux1) q.L 2]);
ebrs = zeros([size(flux1) q.L 2]);

%% Evaluate

for i = 1:size(flux1,1)
    for j = 1:size(flux1,2)
        %         fprintf("Iteration %d of %d\n", j+(i-1)*numel(flux), numel(flux)^2)
        
        q.configChans([flux1(i,j), flux2(i,j)] ,U_K,rho_K);
        
        q.evalAlloc(alloc, rho_K(:,:,1));
        fids(i,j,:, 1) = q.F_L;
        ebrs(i,j,:, 1) = q.ebr_L;
        
        q.evalIsoAlloc(alloc, rho_K(:,:,1));
        fids(i,j,:,2) = q.F_L;
        ebrs(i,j,:,2) = q.ebr_L;
    end
end

% ebrmax = max(ebrs,[],[1,2]);
ebrmax = reshape(q.ebr_Lmax, 1,1,[]);
ebrs_norm = ebrs./ebrmax;
%% Plotting Fidelities, EBR
% Generates a figure for each link for each evaluation mode containing the
% fidelity an ebr as a function of link fluxes

plot_titles = {'Fidelity', 'EBR'};
for i = 1:size(fids,4)
    link_quantites = cat(4, fids(:,:,:,i), ebrs(:,:,:,i));
    for j = 1:q.L
        figure()
        for k = 1:size(link_quantites,4)
            subplot(1,2,k)
            surf(flux1,flux2,link_quantites(:,:,j,k),'LineStyle','none', 'FaceLighting', 'gouraud');
            xlabel("\mu_1")
            ylabel("\mu_2")
            title(plot_titles{k})
        end
        sgtitle(sprintf("Link %d Metrics", j))
    end
end
%% Plotting Candidate Objective Functions
% Generate a 3D plot for every function handle in funcs

w = 0.5; % fidelity-EBR weight

ebrmax_min = min(ebrmax);
% Cell array of objective function handles
funcs = {@(F,R_E) -(w*F + (1-w)*log10(R_E));
    @(F,R_E) w*(1-F) + (1-w)*(1-R_E);
    @(F,R_E) 1 - F + 1./R_E;
    @(F,R_E) ((F-R_E).^2)./(F.*R_E);
    @(F,R_E) (F.*R_E).*(abs((F-R_E))-1);
    @(F,R_E) -F;
    @(F,R_E) -R_E};

% Converting function handles into text, omitting @()
funcstrs = cellfun(@(s)sprintf('%s',func2str(s)),funcs,'UniformOutput',false);
funcstrs = cellfun(@(s)sprintf("$%s$",s(9:end)), funcstrs, 'UniformOutput', false);

% Matrices to hold best fidelities, EBRs, and fluxes to achieve said
% metrics for each cost function
best_fids = zeros(numel(funcs), q.L, 2);
best_ebrs = zeros(numel(funcs), q.L, 2);
best_fluxes = zeros(numel(funcs), q.L, 2);

for i = 1:numel(funcs)
    % Costs = maximum of link objective values for each link pair
    curr_func = funcs{i};
    allf = squeeze(max(curr_func(fids, ebrs_norm), [], 3));
    for j = 1:2 % iterate over evaluation modes
        for l = 1:q.L % iterate over links
            % Find location of minimum
            f_link = allf(:,:,j); % Get link objective values
            [~, idx] = min(f_link,[],'all','linear'); % get linear index of min
            [r_min, c_min] = ind2sub(size(f_link), idx); % Convert linear to 2D
            
            best_fids(i,l,j) = fids(r_min,c_min,l,j);
            best_ebrs(i,l,j) = ebrs_norm(r_min,c_min,l,j);
            best_fluxes(i,:,j) = [flux1(r_min,c_min) flux2(r_min,c_min)];
        end
        
        figure();
        surf(flux1,flux2,allf(:,:,j),'LineStyle','none', 'FaceLighting', 'gouraud')
        hold on
        scatter3(flux1(r_min,c_min), flux2(r_min,c_min), allf(r_min,c_min,j),50, 'dr','filled')
        hold off
        xlabel('\mu_1')
        ylabel('\mu_2')
        title(sprintf("$%s$",funcstrs{i}), 'interpreter', 'latex')
    end
end

%% Plotting model curves and achieved poitns
F = QLANModel.fidUniformState();
RE = QLANModel.ebrUniformState();


for evalmode = 1:2
    figure()
    for i = 1:q.L
        % Get symbolic model functions
        subs_vals = [reshape(q.eff_N(q.links(i,:)),1,[]),...
            reshape(q.d_N(q.links(i,:)),1,[]),...
            q.coinc_wind_s, q.D];
        F = QLANModel.subsParamVars(F, subs_vals);
        REnorm = QLANModel.subsParamVars(RE/ebrmax(i), subs_vals);
        
        % Plot theoretical model bests
        subplot(q.L,1,i)
        fplot(F, [0 3e8]);
        hold on
        fplot(REnorm, [0 3e8]);
        
        % Scatter plot best fidelity, EBR of each objective function
        mkr_shape = {'o','+','*','x','s','d','^','v'};
        for j = 1:numel(funcs)
            X = best_fluxes(j,i,evalmode).*ones(2,1);
            Y = [best_fids(j,i,evalmode); best_ebrs(j,i,evalmode)];
            scatter(X,Y,mkr_shape{mod(j,numel(mkr_shape))})
        end
        ylim([0 1])
        xlabel(sprintf("\\mu_%d",i))
        title(sprintf("Link %s", q.link_labels{i}));
        
    end
    if evalmode == 1
        sgtitle("Simultaneous Link Evaluation")
    else
        sgtitle("Isolated Link Evaluation")
    end
    
    legend([{'Fidelity (Model)'}, {'EBR (Model)'}, funcstrs(:)'], 'location', 'northeastoutside','Interpreter', 'latex')
end
