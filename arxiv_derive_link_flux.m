%% Using reported ArXiv fidelities, calculate necessary link fluxes 
K = 8;
D = 4;
U_K = repmat(eye(D),1,1,K);
rho_psi_m = [0;1;-1;0]*[0 1 -1 0]./2;
rho_K = repmat(rho_psi_m,1,1,K);
chan_fluxes = 1:K;

q = QLAN([], 'arxiv_users.json');
q.coinc_wind_s = 5e-9;
q.configChans(chan_fluxes, U_K, rho_K);
q.link_fluxes = [1e7 1e7 1e7]';

alloc1 = [1;3;3;3;3;3;3;2];
alloc2 = [3;3;1;2;0;0;0;0];
allocs = [alloc1 alloc2];
%% Using multi-objective genetic algo

alloc1_goal_fids = [0.75 0.90 0.55];
alloc2_goal_fids = [0.75 0.84 0.69];
alloc_goals = [alloc1_goal_fids(:) alloc2_goal_fids(:)];

link_fluxes = linspace(1e6, 1e9,10000)';
fids = zeros(numel(link_fluxes),q.L,size(allocs,2));
ebrs = zeros(numel(link_fluxes),q.L,size(allocs,2));
for k = 1:size(allocs,2) % iterate over allocations
    for i = 1:q.L %iterate over links
        for j = 1:numel(link_fluxes) %iterate over link fluxes
            iso_alloc = allocs(:,k);
            iso_alloc(iso_alloc ~= i) = 0;
            
            q.alloc_mat = QLAN.allocVect2Mat(iso_alloc, q.L);
            q.link_fluxes = zeros(q.L+1,1);
            q.link_fluxes(i+1) = link_fluxes(j);
            q.chan_fluxes = q.link_fluxes(i+1)/sum(sign(iso_alloc))*sign(iso_alloc);
            q.calcUserFluxes();
            q.calcUserSingleRate();
            q.calcLinkAccRate();
            q.calcLinkCorrRate();
            q.calcStateWeights();
            q.calcLinkStates();
            q.calcNoisyLinkStates();
            
            ref_rho = repmat(q.rho_K(:,:,1),1,1,q.L);
            q.calcFidelity(ref_rho);
            q.calcEbitRates();
            fids(j,i,k) = q.F_L(i);
            ebrs(j,i,k) = q.ebr_L(i);
        end
    end
end
%%
figure()
plot(link_fluxes, fids(:,:,1), link_fluxes, fids(:,:,2))

%%
for k = 1:size(allocs,2)
    fprintf(1,"Predicted Allocation %d Link Fluxes\n", k);
    for i = 1:q.L
        f = interp1(fids(:,i),link_fluxes,alloc_goals(i,k));
        str = sprintf("Link %s:\t%e\n",q.link_labels{i}, f);
        disp(str)
    end
end
%%
function m = fitness(Q,link_fluxes, alloc, goal)
% for i = 1:Q.L
i = 1;
iso_alloc = alloc;
iso_alloc(iso_alloc ~= i) = 0;

Q.alloc_mat = QLAN.allocVect2Mat(iso_alloc, Q.L);
Q.link_fluxes = [0;link_fluxes(:)];
Q.chan_fluxes = link_fluxes(i)/sum(sign(iso_alloc))*sign(iso_alloc);
Q.calcUserFluxes();
Q.calcUserSingleRate();
Q.calcLinkAccRate();
Q.calcLinkCorrRate();
Q.calcStateWeights();
Q.calcLinkStates();
Q.calcNoisyLinkStates();

ref_rho = repmat(Q.rho_K(:,:,1),1,1,Q.L);
Q.calcFidelity(ref_rho);
Q.calcEbitRates();
m = [Q.F_L(i), Q.ebr_L(i)];


% end

% m = mean((goal(i,:) - eval).^2);
end