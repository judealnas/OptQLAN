function evalAlloc(Q,alloc_vect,ref_rho,varargin)
%EVALALLOC evaluates the given channel allocation for a given QLAN network
%
%   evalAlloc(Q,ALLOC,REF_RHO) evaluates the allocation vector ALLOC by
%   calculating for each link the entangled bit rates and fidelities with
%   respect to density matrix state REF_RHO. If REF_RHO is a DxD matrix,
%   that state is used to calculate fidelity for each link. If REF_RHO is
%   instead a DxDxL matrix, each page of REF_RHO is used to calculate the
%   corresponding link fidelity.
%
%   evalAlloc(__,'EvalMode', mode) selects one of three evaluation methods.
%       'iso' (default) - each link is active and evaluated one at a time.
%       This ignores cross-talk from flux for other entanglement links.
%
%       'batch' - The same as 'iso', however, links are evaluated in
%       batches of non-overlapping links instead of one at a time. May
%       provide evaluation speed-up for large networks.
%
%       'simult' - Assumes all links are simultaneously active. Accounts
%       for cross-talk arising from flux for other entanglement links
valid_alloc = @(x) isvector(x) && numel(x) == Q.K && ... % Vector of K elements
    max(alloc_vect) <= Q.L && min(alloc_vect) >= 0 && ... % no out of range links
    ~any(mod(x,1)); % all integers

valid_refmat = @(x) (isequal(size(ref_rho, [1 2]), size(Q.U_K, [1 2])) &&...
    size(ref_rho,3) == 1 || size(ref_rho,3) == Q.L);

% Use inputParser to get function parameters
p = inputParser;
p.addParameter('EvalMode',QLAN.EVALMODE_ISO,...
    @(x)any(strcmpi(x,{QLAN.EVALMODE_ISO,QLAN.EVALMODE_SIMULT,QLAN.EVALMODE_BATCH})))
p.parse(varargin{:})
mode = lower(p.Results.EvalMode);

% Check required function arguments
if ~valid_alloc(alloc_vect)
    error('Invalid allocation. Allocation must be vector of integers ranging from 0 to %d', Q.L);
end
if ~valid_refmat(ref_rho)
    error("Reference density matrix incompatible with channel unitaries. Reference matrix must be 1 or %d pages of %d-by-%d matrices.", Q.L, size(Q.U_K,1), size(Q.U_K,2));
end

% If single reference state provided, duplicate into 3D matrix
if size(ref_rho,3) == 1
    ref_rho = repmat(ref_rho,1,1,Q.L);
end

% Preallocated matrices to hold evaluation and improve performance; User
% single rates and fluxes have different dimensions depending on evaluation
% method. In 'iso' or 'batch' an extra dimension of size L is added.
% Indexing that dimension gives that metric with that link index isolated.
% E.g. user_fluxes(:,2) gives the user fluxes when only link 2 is active.
Q.link_fluxes = zeros(Q.L, 1);
Q.A_L = zeros(Q.L, 1);
Q.C_L = zeros(Q.L, 1);
Q.lambda_L = zeros(Q.L, 1);
Q.F_L = zeros(Q.L, 1);
Q.ebr_L = zeros(Q.L, 1);
Q.rho_L = zeros([size(Q.rho_K, [1 2]) Q.L]);
Q.sigma_L = zeros([size(Q.rho_K, [1 2]) Q.L]);

lastdim = 1;
if ~strcmpi(mode,QLAN.EVALMODE_SIMULT)
    lastdim = Q.L;
end

Q.user_fluxes = zeros(Q.N, lastdim);
Q.R_N = zeros(Q.N, lastdim);

% Save evaluation mode
Q.eval_mode = mode;

% 'simult' mode of evaluation is called by 'iso' and 'batch' once the
% isolated allocations are constructed
if strcmpi(mode,QLAN.EVALMODE_SIMULT)
    
    noise_L = repmat(eye(Q.D)./Q.D,[1 1 Q.L]);
    
    Q.alloc_vect = alloc_vect(:);
    Q.alloc_mat = QLAN.allocVect2Mat(Q.alloc_vect, Q.L);
    Q.link_fluxes = Q.calcLinkFluxes();
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
    % To accelerate isoalted link evaluation, sets of non-adjacent links
    % are generated. Non-adjacent links are those that do not share users.
    % The links in each set are evaluated. The valid user and link
    % quantities are indexed and stored.
    qtemp = copy(Q); % QLAN copy used for evaluation
    
    
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
        iso_alloc(all(alloc_vect(:) ~= iso_links(:).',2)) = 0;
        
        % Evaluate isolated links
        qtemp.evalAlloc(iso_alloc, ref_rho,'EvalMode',QLAN.EVALMODE_SIMULT);
        
        % Copy quantities from temp QLAN object
        for curr_link = iso_links
            user_idx = Q.links(curr_link,:);
            Q.user_fluxes(user_idx,curr_link) = qtemp.user_fluxes(user_idx);
            Q.R_N(user_idx,curr_link) = qtemp.R_N(user_idx);
        end
        % Storing link quantities
        Q.link_fluxes(iso_links) = qtemp.link_fluxes(iso_links);
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

