profile on
%% Arxiv values
% Column index = allocation
Arxiv_alloc = [1 3
               3 3
               3 1
               3 2
               3 0
               3 0
               3 0
               2 0];
           
Arxiv_F = [0.75 0.75
           0.90 0.84
           0.55 0.69];
       
Arxiv_RE = [56 57
            206 320
            30 26];
        
%% QLAN Instantiation
[names, eff, d] = loadUsersFromFile('arxiv_users.json');
K = size(Arxiv_alloc,1);
psi_m = [0;1;-1;0]/sqrt(2);
rho_psi_m = psi_m*psi_m';
D = size(rho_psi_m,1);
rho_K = repmat(rho_psi_m,[1 1 K]);
U_K = repmat(eye(D)./D, [1 1 K]);
chan_fluxes = rand(K,1);
tau = 5e-9;

q = QLAN([], eff, d, chan_fluxes, U_K, rho_K, tau, 'UserNames', names);
%% Optimize
opts = optimoptions('ga',...
    "PlotFcn", {'gaplotbestf','gaplotbestindiv','gaplotscores'},...
    "CrossoverFraction",0.2);

% Linear Constraint definition
% Requiring chan_flux(i) > chan_flux(i+1)
% A*chan_fluxes returns vector of form [chan_flux(i) - chan_flux(i+1)]
A = zeros(K-1,K);
for i = 1:size(A,1)
    A(i,i:i+1) = [1 -1];
end

[soln, fval] = ga(@(x)fitness(x,Arxiv_alloc,q,Arxiv_F, Arxiv_RE),...
        K,...
        -A,zeros(size(A,1),1),... % Negate A as ga() calculates A*x <= b; Want A*x >= 0
        [],[],... % No Equality constraints
        zeros(K,1),max(q.mu_L_ebrmax,[],'all').*ones(K,1),... % bounds
        [],... % no nonlinear constraints
        opts); 
%%
fitness(soln, Arxiv_alloc,q,Arxiv_F,Arxiv_RE);
q.evalIsoAlloc(Arxiv_alloc(:,1),rho_K(:,:,1));
[q.F_L Arxiv_F(:,1)]
[q.ebr_L Arxiv_RE(:,1)]
%% GA Fitness function
function f = fitness(chan_fluxes,alloc,q,goal_F,goal_RE)
    q.configChans(chan_fluxes, q.U_K, q.rho_K);
    
    F = zeros(size(goal_F));
    RE = zeros(size(goal_RE));
    
    for i = 1:size(alloc,2)
        q.evalIsoAlloc(alloc(:,i),q.rho_K(:,:,1));
        F(:,i) = q.F_L;
        RE(:,i) = q.ebr_L;
    end
    
    f = sum(abs(F-goal_F) + abs(RE-goal_RE), 'all');
end