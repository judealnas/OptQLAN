function evalAlloc(Q,alloc_vect,ref_rho, mode)
%EVALALLOC evaluates the given channel allocation for a given QLAN network
%
%   evalAlloc(Q,ALLOC,REF_RHO) evaluates the allocation vector ALLOC by
%   calculating for each link the entangled bit rates and fidelities with
%   respect to density matrix state REF_RHO. If REF_RHO is a DxD matrix,
%   that state is used to calculate fidelity for each link. If REF_RHO is
%   instead a DxDxL matrix, each page of REF_RHO is used to calculate the
%   corresponding link fidelity.
%
%   evalAlloc(__, mode) selects one of three evaluation methods.
%       'iso' (default) - each link is active and evaluated one at a time.
%       This ignores cross-talk from flux for other entanglement links.
%
%       'batch' - The same as 'iso', however, links are evaluated in
%       batches of non-overlapping links instead of one at a time. May
%       provide evaluation speed-up for large networks.
%
%       'simult' - Assumes all links are simultaneously active. Accounts
%       for cross-talk arising from flux for other entanglement links

%% Parameter validations 
narginchk(3,4);

valid_alloc = @(x) isvector(x) && numel(x) == Q.K && ... % Vector of K elements
    max(alloc_vect) <= Q.L && min(alloc_vect) >= 0 && ... % no out of range links
    ~any(mod(x,1)); % all integers

valid_refmat = @(x) (isequal(size(ref_rho, [1 2]), size(Q.U_K, [1 2])) &&...
    size(ref_rho,3) == 1 || size(ref_rho,3) == Q.L);

valid_mode = @(x)any(strcmpi(x,{QLAN.EVALMODE_ISO,QLAN.EVALMODE_SIMULT,QLAN.EVALMODE_BATCH}));

% Check required function arguments
if ~valid_alloc(alloc_vect)
    error('Invalid allocation. Allocation must be vector of integers ranging from 0 to %d', Q.L);
end
if ~valid_refmat(ref_rho)
    error("Reference density matrix incompatible with channel unitaries. Reference matrix must be 1 or %d pages of %d-by-%d matrices.", Q.L, size(Q.U_K,1), size(Q.U_K,2));
end

if nargin == 4 
    if ~valid_mode(mode)
        error("Unrecognized mode. Options are '%s', '%s', or '%s'",...
            QLAN.EVALMODE_ISO, QLAN.EVALMODE_SIMULT, QLAN.EVALMODE_BATCH)
    end
elseif nargin < 4
    mode = QLAN.EVALMODE_ISO;
end
%% Allocation Evaluation
% If single reference state provided, duplicate into 3D matrix
if size(ref_rho,3) == 1
    ref_rho = repmat(ref_rho,1,1,Q.L);
end

% Persistent noise variable
persistent noise_L
if isempty(noise_L)
    noise_L = repmat(eye(Q.D)./Q.D,[1 1 Q.L]);
end

% Save evaluation mode
Q.eval_mode = mode;

% Calculation of link fluxes same regardless of evaluation mode
Q.alloc_vect= alloc_vect(:);
Q.alloc_mat = QLAN.allocVect2Mat(Q.alloc_vect,Q.L);
Q.link_fluxes = Q.calcLinkFluxes();

% 'simult' mode of evaluation is called by 'iso' and 'batch' once the
% isolated allocations are constructed
if strcmpi(mode,QLAN.EVALMODE_SIMULT)
    
    Q.user_fluxes = Q.calcUserFluxes();
    Q.R_N = Q.calcUserSingleRate();
    Q.A_L = Q.calcLinkAccRate();
    Q.C_L = Q.calcLinkCorrRate();
    Q.lambda_L = Q.calcStateWeights();
    Q.rho_L = Q.calcLinkStates();
    Q.sigma_L = Q.calcNoisyLinkStates(noise_L);
    Q.F_L = Q.calcFidelity(ref_rho);
    Q.ebr_L = Q.calcEbitRate();
    
    return
else
    % Preallocated matrices to hold non-simultaneous evaluation and improve
    % performance; User single rates and fluxes have different dimensions
    % depending on evaluation method. In 'iso' or 'batch' an extra
    % dimension of size L is added. Indexing that dimension gives that
    % metric with that link index isolated. E.g. user_fluxes(:,2) gives all
    % user fluxes when only link 2 is active.
    Q.A_L = zeros(Q.L, 1);
    Q.C_L = zeros(Q.L, 1);
    Q.lambda_L = zeros(Q.L, 1);
    Q.F_L = zeros(Q.L, 1);
    Q.ebr_L = zeros(Q.L, 1);
    Q.rho_L = zeros([size(Q.rho_K, [1 2]) Q.L]);
    Q.sigma_L = zeros([size(Q.rho_K, [1 2]) Q.L]);   
    Q.user_fluxes = zeros(Q.N, Q.L);
    Q.R_N = zeros(Q.N, Q.L);
    
    % A copy of the provided QLAN object is used to evaluate the isolated
    % allocations
    qtemp = copy(Q); 
    
    % construct link batches as a cell array of integer arrays
    if strcmpi(mode,QLAN.EVALMODE_ISO)
        % Evaluate links one at a time
        link_batches = num2cell(1:Q.L);
    elseif strcmpi(mode,QLAN.EVALMODE_BATCH)
        % Evaluate as non-adjacent batches
        link_batches = nonAdjEdgeSet(Q.links);
    else
        error("Unrecognized evaluation mode %s", mode)
    end
    
    for i = 1:numel(link_batches)
        % Get isolated links in batch
        iso_links = link_batches{i};
        
        % Construct isolated allocation vector by setting all link indices
        % not in current batch to 0
        iso_alloc = alloc_vect;
        iso_alloc(all(iso_alloc(:) ~= iso_links(:).',2)) = 0;
        
        % Evaluate isolated links
        qtemp.evalAlloc(iso_alloc, ref_rho,QLAN.EVALMODE_SIMULT);
        
        % Copy quantities from temp QLAN object
        for curr_link = iso_links
            user_idx = Q.links(curr_link,:);
            Q.user_fluxes(user_idx,curr_link) = qtemp.user_fluxes(user_idx);
            Q.R_N(user_idx,curr_link) = qtemp.R_N(user_idx);
        end
        Q.A_L(iso_links) = qtemp.A_L(iso_links);
        Q.C_L(iso_links) = qtemp.C_L(iso_links);
        Q.lambda_L(iso_links) = qtemp.lambda_L(iso_links);
        Q.rho_L(:,:,iso_links) = qtemp.rho_L(:,:,iso_links);
        Q.sigma_L(:,:,iso_links) = qtemp.sigma_L(:,:,iso_links);
        Q.F_L(iso_links) = qtemp.F_L(iso_links);
        Q.ebr_L(iso_links) = qtemp.ebr_L(iso_links);
    end
end
end

