QLAN Optimization Test File

Parameters
N = 4; % Max number of users
K = 12; % Number of frequency channel pairs
D = 4; % quantum state dimensionality; dimension of two-qubit system
coinc_wind_s = 1.024e-9; % detector coincidence window in seconds


% list of node pairs constituting network links
links = [1,2;...
         1,3;...
%          1,5;...
         2,3;...
         2,4;...
%          2,5;...
         3,4;...
%          4,5;...
  ];
L = size(links,1); % Number of links in specified network

% Constructing edge table with edge labels
h = @(stem) sprintf("l%s", num2str(stem));
link_names = arrayfun(h, (1:L)', 'UniformOutput', 1);
EdgeTable = table(links, link_names, 'VariableNames', {'EndNodes', 'EdgeNames'});

% Construct table of node labels
h = @(stem) sprintf("n%s", num2str(stem));
user_names = arrayfun(h, (1:N)', 'UniformOutput', 1);
NodeTable = table(user_names, 'VariableNames', {'Name'});

% construct graph object using edge table and node names
g = graph(EdgeTable, NodeTable);

% Plot graph with names
fig = plot(g);
fig.EdgeLabel = g.Edges.EdgeNames;
inc_mat = abs(incidence(g));
disp(full(inc_mat))

User Parameters
% Each user assumed to have two detectors of equal efficiencies
eff_N = [0.85; 0.85; 0.2; 0.1]; % vector of user detector efficiencies
d_N = zeros(size(eff_N))

Channel Parameters
% Each channel has an initial state rho_k
psi_m = [0;1;-1;0]; % EPR state vector |HV> - |VH>
rho_psi_m = psi_m*psi_m'; % construct density matrix for initial entangled states
rho_K = repmat(rho_psi_m, 1, 1, K); % initial channel states are all |HV>-|VH>

% Each channel has unitary describing change in quantum states after
% transmission. U_K(:,:,k) returns the unitary of the kth channel pair
U_K = repmat(eye(D)/D,[1,1,K]); % Ideal case where channel unitaries are identities

% Define spectrum partitioning; central stop band assumed
wss_f0 = 192.3125e12; % WSS center frequency
chan_width = 25e9;
power_spect = gaussPowerSpect(wss_f0,610e9,24e-3);

Calculating channel fluxes
% Each channel has associated photon flux which can be calculated from
% input power spectrum and spectrum partitioning

% Matrices of channel boundaries
lo_freq_bounds = zeros(K,2); % freq bounds of lower spectral slices
hi_freq_bounds = zeros(K,2); % freq bounds of upper spectral slices

for k = 1:K
        lo_freq_bounds(k,:) = [wss_f0 - chan_width*(k+1/2), wss_f0 - chan_width*(k-1/2)];
        hi_freq_bounds(k,:) = [wss_f0 + chan_width*(k-1/2), wss_f0 + chan_width*(k+1/2)];
end


chan_flux = zeros(K, 1);
for k = 1:K
    chan_flux(k) = integral(power_spect, min(lo_freq_bounds(k,:)), max(lo_freq_bounds(k,:)));
    chan_flux(k) = chan_flux(k) + ...
                   integral(power_spect, min(hi_freq_bounds(k,:)), max(hi_freq_bounds(k,:)));
end

plot(freqs, power_spect(freqs))
xline(wss_f0, 'r--')
xline(unique(bound_freqs(:)), 'b:')

Constructing allocation matrix from allocation vector
% Random allocation vector; index with channel number to obtain assigned
% link number
alloc_vect = zeros(K,1);
for k = 1:K
    alloc_vect(k) = mod(k-1, L+1);
end
disp(alloc_vect)
[GC, GR] = groupcounts(alloc_vect)
% Sparse matrix representation of allocation vector where row indices are
% link labels and column indices are channel labels.
% E.g., if alloc_mat(3,5) == 1, then channel 5 is assigned to link 2.
% if alloc_mat(1,2) == 1, then channel 2 is assigned to link 0 (null link)
alloc_mat = double((0:L)' == alloc_vect');

% Link flux = sum(fluxes of allocated channels)
% link_fluxes = total photon fluxes corresponding to each link
link_fluxes = alloc_mat*chan_flux;

% User flux = sum(non-null link fluxes connected to user)
user_fluxes = inc_mat*link_fluxes(2:end);


Calculating count rates
R_N = eff_N/2 .* user_fluxes + d_N % vector of user single count rates
A_L = zeros(L,1); % vector of link accidental rates (excluding null)
C_L = zeros(L,1); % vector of link true coincidence rates (excluding null)
full_inc_mat = logical(full(inc_mat)); % ensure incidence matrix not in sparse form
for i = 1:L
    % for each link, calculate product of single count rates for connected
    % users and coincidence window
    A_L(i) = 8*coinc_wind_s*prod(R_N(full_inc_mat(:,i)));
    
    % for each link calculate product of efficiencies of connceted users
    % and link photon flux
    C_L(i) = link_fluxes(i)*prod(eff_N(full_inc_mat(:,i)));
end

disp(A_L)
disp(C_L)

Calculating Received State
% Noise-entanglement weight per link
lambda_L = 1/(1+A_L./C_L);
lambda_L_3D = reshape(lambda_L,1,1,L);



% Calculating noise state = 1-lamda*eye(D)/D for each link
noise_L = pagemtimes(1-lambda_L_3D, repmat(eye(D)/D, 1, 1, L));
% disp(noise_L)

% rho_L = density matrix shared across links
rho_L = zeros(D,D,L);
for i = 1:L
    U = U_K(:,:,logical(alloc_mat(i,:))); % unitaries of channels allocated to link i
    rho = rho_K(:,:,logical(alloc_mat(i,:))); % initial states of channels allocated to link i
    flux = chan_flux(logical(alloc_mat(i,:)));
    flux = reshape(flux,1,1,numel(flux));
    
    % calculate sum of all U*rho*U' density matrices
    rho_L(:,:,i) = sum(pagemtimes(flux,pagemtimes(pagemtimes(U,rho),'none',U,'ctranspose')),3);
end

% weighting entangled state using lambda and portion of total flux
link_fluxes_3D = reshape(link_fluxes(2:end),1,1,L);
epr_weights = pagemtimes(lambda_L_3D, link_fluxes_3D); % lambda_l / mu_L_l
rho_L = pagemtimes(epr_weights, rho_L); % lambda_l*rho_l

sigma_L = rho_L + noise_L;

disp(noise_L)
disp(rho_L)
disp(sigma_L)
% verifying density matrices still satisfy trace = 1
traces = zeros(L,1);
for i = 1:L
    traces(i) = trace(sigma_L(:,:,i));
end
disp(traces)

Fidelity Calculation
F_L = zeros(L,1);
for i = 1:L
    F_L(i) = psi_m'*sigma_L(:,:,i)*psi_m;    
end
disp(F_L)
