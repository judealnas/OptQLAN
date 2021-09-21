%% Calculating Channel Fluxes from Single Counts Reported in arXiv Experiment
clc
%% Parameters
K = 8; % # of channel pairs
D = 4; % System dimensions
h = 6.626e-34; % Planck's constant

% Each channel has an initial state rho_k
psi_m = [0;1;-1;0]./sqrt(2); % EPR state vector |HV> - |VH>
rho_psi_m = psi_m*psi_m'; % construct density matrix for initial entangled states
rho_K = repmat(rho_psi_m, 1, 1, K); % initial channel states are all |HV>-|VH>

% Each channel has unitary describing change in quantum states after
% transmission. U_K(:,:,k) returns the unitary of the kth channel pair
U_K = repmat(eye(D),[1,1,K]); % Ideal case where channel unitaries are identities


%% instantiate fully connected 3-user network
q = QLAN([], 'arxiv_users.json');
q.coinc_wind_s = 5e-9;
q.configChans(1:K, U_K, rho_K);

link_AB = q.g.findedge('Alice', 'Bob');
link_AC = q.g.findedge('Alice', 'Charlie');
link_BC = q.g.findedge('Charlie', 'Bob');

alloc1 = uint8([link_AB;link_BC;link_BC;link_BC;link_BC;link_BC;link_BC;link_AC]);
alloc2 = uint8([link_BC;link_BC;link_AB;link_AC;0;0;0;0]);
alloc_vects = [alloc1 alloc2];

alloc_mat1 = double(QLAN.allocVect2Mat(alloc1, q.L));
alloc_mat2 = double(QLAN.allocVect2Mat(alloc2, q.L));
alloc_mats = cat(3, alloc_mat1, alloc_mat2);

% disp('Allocation Matrices')
alloc_mats;
%% Defining Experimental Data
% Calculating total single user counts (scr) from values reported in arXiv paper
arxiv_link_scr1 = [1.2e5 5e3;   % AB scr
    4.7e4 1e4;   % AC scr
    7.5e3 1.4e5];% BC scr

arxiv_link_scr2 = [1.2e5 5.1e3; % AB scr
    8.4e3 1.9e4; % AC scr
    6.4e3 6.7e4];% BC scr

arxiv_scaled_scr1 = nchoosek(2./q.eff_N,2).*(arxiv_link_scr1 - nchoosek(q.d_N,2));
arxiv_scaled_scr2 = nchoosek(2./q.eff_N,2).*(arxiv_link_scr2 - nchoosek(q.d_N,2));

arxiv_scaled_scr = cat(3,arxiv_scaled_scr1,arxiv_scaled_scr2);
%% Deriving Channel Fluxes

% Column 1 of flux vectors correspond to using min; Col 2 to max single
% rate
arxiv_link_fluxes = cat(2,min(arxiv_scaled_scr,[],2),max(arxiv_scaled_scr,[],2));

% Row index = Channel number; Col index = min/max; Pager idx = Alloc 1/2
arxiv_chan_fluxes = zeros(K,size(arxiv_link_fluxes,2),size(alloc_mats,3));
for i = 1:size(arxiv_chan_fluxes,3)
    arxiv_chan_fluxes(:,:,i) = pinv(alloc_mats(2:end,:,i))*arxiv_link_fluxes(:,:,i);
end
%% Evaluating Allocation 1 using Derived Chan Fluxes
Fid = zeros([q.L size(arxiv_chan_fluxes, [2 3]) 2]);
ebr = zeros(size(Fid));
ref_rho = rho_K(:,:,1);

% iterate over allocations to evaluate
for eval_ind = 1:size(alloc_vects,2)
    % iterate over allocations to use in channel flux derivation
    for alloc_ind = 1:size(alloc_mats,3)
        % iterate over which single count data derivation to use
        for scr_ind = 1:size(arxiv_link_fluxes,2)
            %Choosing channel flux
            %             fprintf(1,"-----\nalloc_ind=%d\tscr_ind=%d\teval_ind=%d\n",...
            %                              alloc_ind, scr_ind, eval_ind);
            chan_flux = arxiv_chan_fluxes(:,scr_ind,alloc_ind);
            q.configChans(chan_flux, U_K, rho_K);
            eval_alloc = alloc_vects(:,eval_ind);
            for link_ind = 1:max(eval_alloc)
                iso_vect = eval_alloc;
                iso_vect(iso_vect ~= link_ind) = 0;
                m = q.evalAlloc(iso_vect, ref_rho);
                Fid(link_ind,scr_ind, alloc_ind, eval_ind) = m(link_ind,1);
                ebr(link_ind,scr_ind, alloc_ind, eval_ind) = m(link_ind,2);
                
            end
        end
    end
end
%% Displaying data
for alloc_ind = 1:size(Fid,3)
    for scr_ind = 1:size(Fid,2)
        for eval_ind = 1:size(Fid,4)
            
            x = {'min', 'max'};
            scr_choice = x{scr_ind};
            alloc_choice = sprintf('Allocation %d', alloc_ind);
            eval_choice = sprintf('Allocation %d', eval_ind);
            
            fprintf(1, "Metrics for %s evaluated using %s, %s SCR\n\n",...
                eval_choice, alloc_choice, scr_choice)
            disp('Channel Fluxes')
            disp(arxiv_chan_fluxes(:,scr_ind, alloc_ind)')
            disp('Fidelity')
            disp(Fid(:,scr_ind,alloc_ind,eval_ind)')
            disp('EBR')
            disp(ebr(:,scr_ind,alloc_ind,eval_ind)')
        end
    end
end