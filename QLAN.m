classdef QLAN < handle & matlab.mixin.Copyable
    %QLAN Representation of flex-grid quantum LAN
    %   Q = QLAN(LINKS, FILEPATH) retrieves an array of users from the
    %   JSON-encoded file FILEPATH. LINKS is an Lx2 numeric matrix where
    %   each row contains the indices of users forming a network link.
    %
    %   Q = QLAN(LINKS, EFF_N, D_N) constructs the network using the
    %   corresponding arrays of user efficiencies EFF_N and detector dark
    %   count rates D_N.
    %
    %   Q = QLAN(_,'UserNames', UNAMES) associates names with users of the
    %   network where NAMES is a string vector or cell array of character
    %   vectors withe the same number of elements as users in the network,
    %   i.e. size(EFF_N,1) == numel(UNAMES)
    %
    %   Q = QLAN(_,'LinkNames', LNAMES) associates names with the links of
    %   network. LNAMES is a string vector or cell array of character
    %   vectors with the same number of elements as links in the network,
    %   i.e. size(LINKS,1) == numel(LNAMES)
    
    %% Public Properties
    properties
        alloc_vect (:,1) uint8; % Allocation vector
        coinc_wind_s (1,1) double = 0; % Coincidence window in seconds
        user_names (:,1); % network user labels
        link_labels (:,1); % link names
        inc_mat (:,:) logical; % Network incidence matrix in sparse format
        eff_N (:,1); % Col vector of user detector efficiencies
        d_N (:,1);
        U_K; % DxDxK matrix of channel unitaries
        rho_K; % DxDxK matrix of initial channel states
        chan_fluxes (:,1);
    end
    
    %% Dependent properties
    % Quanitites derived from other properties' set methods
    properties (Dependent=true)
        N; % # of users; calculated as size(eff_N,1);
        K; % # of channel pairs; size(chan_fluxes,1);
        L; % # of non-null network links; size(links,1);
    end
    
    %% Set-protected properties
    % Calculated quantities
    properties (SetAccess=public)
        links (:,2) {};
        user_file (1,:) char = ''; % File from which users were imported
        alloc_mat logical; % matrix expansion of allocation vector
        link_fluxes; % (L+1)x1 vector of link biphoton fluxes
        user_fluxes; % Nx1 vector of user single photon fluxes
        R_N; % Nx1 vector of single photon rates per detector per user
        A_L; % Accidental coincidence rate per link
        C_L; % Correlated coincidence rate per link
        lambda_L; % weight between entangled and noise states
        rho_L; % link states prior to adding noise
        sigma_L; % link states after weighting and noise addition
        F_L; % Per link fidelities
        ebr_L; % Per link ebit rates
        
    end
    
    %% Static Methods
    methods (Static=true)
        function M = allocVect2Mat(V,L)
            %ALLOCVECT2MAT converts an allocation vector to an allocation
            %matrix
            %
            %   M = allocVect2Mat(V,L) converts allocation vector V to
            %   allocation matrix M given L non-null links. Each row of M
            %   corresponds to a link. Each column of M corresponds to a
            %   channel pair. If M(l,k) == 1, then the kth channel is
            %   allocated to the lth link, where the 0th link is the null
            %   link.
            
            if ~(isvector(V) && all(mod(V,1) == 0, 'all'))
                error("V must be a vector of integers >= 0") 
            end
            M = logical((0:L)' == V(:)');
        end
        
        function link_fluxes = calcLinkFluxes_(alloc_mat, chan_fluxes)
            %CALCLINKFLUXES_ is a static method that returns the flux per
            %network link 
            %   link_fluxes = calcLinkFluxes_(ALLOC_MAT, CHAN_FLUX) calculates
            %   the flux through the links as the matrix multiplication of the
            %   allocation  matrix and the channel pair flux vector:
            %       link_fluxes = ALLOC_MAT * CHAN_FLUX
            
            link_fluxes = alloc_mat * chan_fluxes;
        end
        
        function user_fluxes = calcUserFluxes_(incidence, link_fluxes)
            %CALCUSERFLUXES_ is a static method for calculating flux
            %received by each user in network
            %
            %   user_fluxes = calcUserFluxes_(INCID, LINK_FLUX) calculates
            %   the flux received by network users via matrix
            %   multiplication of the incidene matrix and link flux vector.
            %       INCID - NxL incidence matrix
            %       LINK_FLUX - (L+1)x1 vector of L non-null-link fluxes
            %       and 1 null link flux at index 1
            user_fluxes = incidence * link_fluxes(2:end);
        end
        
        function R_N = calcUserSingleRate_(user_fluxes, eff_N, d_N)
            %CALCUSERSINGLERATE_ is a static method calculating the single
            %count rate at a user's detector
            %
            %   R_N = calcUserSingleRate_(USER_FLUX, EFF_N, D_N) returns
            %   the single count rate given the user flux vector USER_FLUX,
            %   user efficiency vector EFF_N, and user dark count rate
            %   vector D_N. All arguments must have same number of
            %   elements.
            R_N = eff_N.*user_fluxes./2 + d_N;
        end
        
        function A_L = calcLinkAccRate_(inc_mat, R_N, coinc_wind_s)
            %CALCLINKACCRATE_ is a static method calculating the accidental
            %coincidence rate of network links
            %
            %   A_L = calcLinkAccRate_(INCID, R_N, TAU) calculates the link
            %   accidental rates A_L given:
            %       INCID - NxL undirected incidence matrix
            %       R_N - Nx1 vector of user single count rates
            %       TAU - scalar double representing coincidence window in
            %       seconds.
            %       A_L - Lx1 vector of accidental coincidence rates for
            %       each link
            L = size(inc_mat,2);
            A_L = zeros(L, 1);
            full_inc_mat = logical(full(inc_mat));
            for i = 1:L
                A_L(i) = ...
                    8*coinc_wind_s.*prod(R_N(full_inc_mat(:,i)));
            end
        end
        
        function C_L = calcLinkCorrRate_(inc_mat, link_fluxes, eff_N)
            %CALCLINKCORRRATE_ is a static method calculating the true,
            %correlated coincidence rate of each link in a network.
            %
            %   C_L = calcLinkCorrRate_(INCID, LINK_FLUX, EFF_N) returns
            %   the link true coincidence rate vector C_L.
            %       INCID - NxL undirected incidence matrix
            %       LINK_FLUX - (L+1)x1 vector of L non-null-link fluxes
            %                   and one null link flux at index 1
            %       EFF_N - Nx1 vector of user detector efficiencies
            %       C_L - Lx1 vector of true coincidence rates for each
            %               link
            L = size(inc_mat,2);
            C_L = zeros(L, 1);
            full_inc_mat = logical(full(inc_mat));
            for i = 1:L
                C_L(i) = ...
                    link_fluxes(i+1)*prod(eff_N(full_inc_mat(:,i)));
            end
        end
        
        function lambda_L = calcStateWeights_(A_L, C_L)
            %CALCSTATEWEIGHTS_ is a static method that calculates the
            %relative weights between the maximally mixed state and true
            %desired state per link
            %
            %   lambda_L = calcStateWeights_(A_L,C_L) calculates the weight
            %   vector as 
            %       lambda_L = 1 ./ (1 + A_L ./ C_L)
            %   where A_L and C_L are the link accidental and true
            %   coincidence rates respectively
            lambda_L = 1./(1 + A_L ./ C_L);
        end
        
        function rho_L = calcLinkStates_(alloc_mat, link_fluxes,...
                                         chan_fluxes, U_K, rho_K)
            %CALCLINKSTATES_ is a static method that calculates the state
            %shared on a network link prior to the introduction of noise
            %
            %   rho_L calcLinkStates_(ALLOC_MAT, LINK_FLUX, CHAN_FLUX,
            %                         U_K, RHO_K)
            %       ALLOC_MAT - LxK allocation matrix
            %       LINK_FLUX - (L+1)x1 vector of L non-null-link fluxes
            %                   and one null link flux at index 1
            %       CHAN_FLUX - Kx1 vector of channel pair fluxes
            %       U_K       - DxDxK 3D matrix of channel unitaries
            %       RHO_K     - DxDxK 3D matrix of initial channel states
            %   The allocation matrix is used to deterine which channels
            %   are allocated to a link. The unitaries U_k of those
            %   channels are applied to their respective intial channel
            %   states rho_k. These states are summed, weighting each
            %   channel state by the channel's flux contribution to the
            %   link flux.
            %   
            L = size(alloc_mat,1) - 1; % get number of non null links
            rho_L = zeros([size(rho_K,[1 2]) L]);
            for i = 2:(L+1) % iterate over non-null links
                alloc_chans = logical(alloc_mat(i,:));
                U = U_K(:,:,alloc_chans); % unitaries of channels allocated to link i-1
                rho = rho_K(:,:,alloc_chans); % initial states of channels allocated to link i-1
                chan_flux = chan_fluxes(alloc_chans); % get channel fluxes allocated to link i-1
                chan_flux = reshape(chan_flux,1,1,numel(chan_flux)); % reshape to 3D vector for quicker math
                
                % Calculate sum of all chan_flux*U*rho*U' density matrices
                % associated with link and normalize by total link flux
                rho_L(:,:,i-1) = sum(pagemtimes(chan_flux,pagemtimes(pagemtimes(U,rho),'none',U,'ctranspose')),3);
                rho_L(:,:,i-1) = rho_L(:,:,i-1)./link_fluxes(i);
            end
        end
        
        function sigma_L = calcNoisyLinkStates_(rho_L, noise_L, lambda_L)
            %CALCNOISYLINKSTATES is a static method that calculates the
            %link states after inclusion of noise.
            %
            %   sigma_L = calcNoisyLinkStates_(RHO, NOISE, LAMBDA) returns
            %   the sum of the corresponding RHO's and NOISE's weighted by
            %   the corresponding LAMBDA and 1-LAMBDA respectively
            %       RHO - DxDxL 3D matrix of link states
            %       NOISE - DxDxL 3D matrix of noise states
            %       LAMBDA - Lx1 vector of weights
            %       sigma_L - DxDxL 3D matrix of noisy link states
            %           sigma_L(:,:,l) = 
            %               lambda_L(l)*rho_L(:,:,l) +
            %               (1-lambda_L(l))*noise_L(:,:,l)
            lambda_L = reshape(lambda_L, 1,1, numel(lambda_L));
            sigma_L = pagemtimes(lambda_L, rho_L) + ...
                      pagemtimes(1-lambda_L, noise_L);
            sigma_L(isnan(sigma_L)) = 0;
        
        end
        
        function F_L = calcFidelity_(rho_L, psi_L)
            %CALCFIDELITY_ is a static method that calculates the
            %fidelities between pages of provided states
            %
            %   F_L = calcFidelity_(RHO, PSI) returns a vector of
            %   fidelities between the corresponding matrices of the 3D
            %   matrices RHO and PSI. RHO and PSI must have the same size
            %   third dimension (i.e. same number of pages) 
            %       RHO - 3D matrix of L pages of matrices
            %       PSI - 3D matrix of L pages of matrices
            %       F_L - Lx1 vector of fidelities between corresponding
            %             RHO's and PSI's
            if size(rho_L,3) ~= size(psi_L,3)
                error("3rd dimension of args must be equal")
            end
            
            L = size(rho_L,3);
            
            F_L = zeros(L, 1);
            for i = 1:L
                F_L(i) = fidelity(rho_L(:,:,i), psi_L(:,:,i));
            end
        end
        
        function ebr_L = calcEbitRate_(A_L, C_L, sigma_L)
            %CALCEBITRATE_ is a static method that calculates the entangled
            %bit rate over a network link
            %   
            %   ebr_L = calcEbitRate_(A_L, C_L, SIGMA_L) calculates the
            %   link entangled bit rate using the link coincidence rates
            %   A_L, C_L and the log negativity of the noisy link states
            %   SIGMA_L
            %       A_L - Lx1 vector of link accidental coincidence rates
            %       C_L - Lx1 vector of link truue coincidence rates
            %       SIGMA_L = DxDxL 3D matrix of noisy link states used for
            %       log negativity calculation
                    
            L = size(sigma_L,3);
            ebr_L = zeros(L,1);
            for i = 1:L
                ebr_L(i) = logneg(sigma_L(:,:,i));
            end
            
            ebr_L = (A_L + C_L).*ebr_L;
        end
    end
    
    %% Instance Methods
    methods
        %% Configuration Methods
        function obj = QLAN(varargin)
            
            % Check if arg is string or character vector
            istext = @(x) isstring(x) || ischar(x);
            
            % Check if arg is string vect or cell string array
            istextvector = @(x) isvector(x) && (isstring(x) || iscellstr(x));
            
            % Check if arg is numeric vector
            isnumcvector = @(x)isvector(x) && isnumeric(x) ;
            
            p = inputParser;
            p.addRequired('links', @(x)isempty(x) || (size(x,2) == 2 && isnumeric(x)));
            p.addRequired('file_eff', @(x) istext(x) || isnumcvector(x));
            p.addOptional('d_N', [], isnumcvector);
            p.addParameter('UserNames',[], istextvector);
            p.addParameter('LinkNames',[], istextvector);
            
            p.parse(varargin{:});
            links = p.Results.links;
            file_eff = p.Results.file_eff;
            d_N = p.Results.d_N;
            user_names = p.Results.UserNames(:);
            link_names = p.Results.LinkNames(:);
                       
            % Getting users from file file_eff arg is text
            if istext(file_eff)
                filename = file_eff;
                [file_user_names, eff_N, d_N] = loadUsersFromFile(filename);
                
                % If no user names provided, use names from file
                if isempty(user_names)
                    user_names = file_user_names;
                end
                
            % Otherwise file_eff is numeric vector of efficiencies
            else
                filename = '';
                eff_N = file_eff(:);
            end
            
            % Save user file, efficiencies, and dark cnts vales
            obj.user_file = filename;
            obj.eff_N = eff_N;
            obj.d_N = d_N;
            
            % Get # of users
            N = size(eff_N,1);
            
            % Generae default user names
            default_user_names = arrayfun(@(x)sprintf('n%d',x), (1:N)',...
                'UniformOutput', false);
            
            % If empty list provided for links, use fully connected net
            if isempty(links)
                links = nchoosek(1:N, 2);
            end
            
            % Check if links include unspecified user
            if max(links, [], 'all') > N
                error("File '%s' contains %d users but provided links assumes %d users.",...
                    file_eff, N, max(links, [], 'all'))
            end
            
            % If user_names empty, use default names 
            if isempty(user_names)
                user_names = default_user_names;
            
            % if insufficient names, pad using default names
            elseif numel(user_names) < N
                warning("%d node names provided when expecting %d. Padding with default names",...
                    numel(user_names), N);
                
                user_names = [cellstr(user_names) ; default_user_names(numel(user_names)+1:end)];
            
            % If excess names, ignore excess
            elseif numel(user_names) > N
                warning("Too many user names specified. Extra names ignored");
                user_names = user_names(1:N);
            end
            obj.user_names = user_names; % Saving user names 
            
            % If no link names provided, use user names to generate names
            if isempty(link_names)
                obj.configLinks(links, 'GenLabels', user_names);
            else
                obj.configLinks(links, link_names);
            end
        end
        
        function obj = configChans(obj, chan_fluxes, chan_unitry, chan_states)
            %CONFIGCHANS(Q, chan_fluxes, chan_unitry, chan_states)
            %configures the following QLAN properties:
            %   chan_fluxes - channel biphoton flux
            %   chan_unitry - (:,:,K) matrix of channel unitaries
            %   chan_states - (:,:,K) matrix of initial channel density matrices
            chan_fluxes = chan_fluxes(:);
            sz1 = [size(chan_fluxes,1), size(chan_unitry,3), size(chan_states,3)];

            if any(sz1 ~= sz1(1))
                error("Dimension Mismatch: Unequal channels")
            end
            
            if any(size(chan_unitry, [1 2]) ~= size(chan_states, [1 2]))
                error("Dimension Mismatch: unitaries and states are incompatible sizes");
            end
            
            obj.U_K = chan_unitry;
            obj.rho_K = chan_states;
            obj.chan_fluxes = chan_fluxes;
        end
        
        function obj = configLinks(obj, links, varargin)
            %CONFIGLINKS re-connects users in a QLAN to form an entanglement
            %   network. 
            %
            %   configLinks(QLAN, LINKS) - connects users in the
            %   QLAN object according to the Lx2 numeric matrix LINKS. Each
            %   row of LINKS contains the indices of users forming the link
            %   endpoints.
            %   
            %   configLinks(QLAN, LINKS, LABELS) uses the string vector or
            %   cell vector of char vectors LABELS as link labels for
            %   plotting
            %
            %   configLinks(QLAN, LINKS, 'GenLabels', USERNAMES) generates
            %   link labels using the matrix LINKS and the string vector or
            %   cell vector of char vectors USERNAMES
                        
            % Check if arg is string or cell string array
            istextvector = @(x)isvector(x) && (isstring(x) || iscellstr(x));
            
            p = inputParser;
%             p.FunctionName = 'configLinks';
            p.addRequired('obj', @(x) isa(x,'QLAN'));
            p.addRequired('links', @(x) size(x,2) == 2 && ...
                                   isnumeric(x) && ...
                                   max(x,[],'all') <= obj.N);
            p.addOptional('labels',[], istextvector);
            p.addParameter('GenLabels', [], istextvector);
            
            p.parse(obj,links,varargin{:});
            obj = p.Results.obj;
            links = p.Results.links;
            labels = p.Results.labels(:);
            unames = p.Results.GenLabels(:);
            
            n = size(links,1); 
            
            % Generating default link labels
            if isempty(unames)
                default_labels = arrayfun(@(x)sprintf('l%d',x), (1:n)',...
                    'UniformOutput', false);
            else
                default_labels = generateLinkLabels(unames, links);
            end
            
            if ~isempty(labels)
                if numel(labels) < n 
                    warning("Insufficient link labels provided");
                    labels = [cellstr(labels)
                              default_labels(numel(labels)+1:end)];
                elseif numel(labels) > n
                    warning("Excess link labels ignored")
                    labels = labels(1:n);
                end
            else
                labels = default_labels;
            end
            
            obj.inc_mat = linksToIncidence(links);
            obj.links = links;
            obj.link_labels = labels;
        end
        
        function tf = isConfiged(obj)
            %ISCONFIGED checks if the necessary QLAN parameters have been
            %set in order to calculate network metrics.
            %
            %   tf = isConfiged(QLAN) returns true if the following
            %   parameters meet the corresponding conditions:
            %       - 
            
            persistent cfg_msgs;
            cfg_msgs = {
                'Invalid coincidence window';
                'Channel unitaries unspecified';
                'Channel states unspecified';
                'Channel fluxes unspecified';
                'Network links unspecified';
                'User efficiencies unspecified';
                'User dark count rates unspecified';
            };
            
            tf_arr = [obj.coinc_wind_s > 0;
            ~isempty(obj.U_K);
            ~isempty(obj.rho_K);
            ~isempty(obj.chan_fluxes);
            ~isempty(obj.inc_mat);
            ~isempty(obj.eff_N);
            ~isempty(obj.d_N)];
            
            tf = all(tf_arr);
            if ~tf
                disp(cfg_msgs{tf_arr})
            end
                
        end
        
        %% Dependent property accessors
        function k = get.K(obj)
            k = numel(obj.chan_fluxes);
        end
        
        function n = get.N(obj)
            n = size(obj.eff_N, 1);
        end
        
        function l = get.L(obj)
            l = size(obj.links,1);
        end

        %% User File Loading/Saving methods
        function obj = saveUsersToFile(obj, filename)
            saveUsersToFile(filename, obj.eff_N, obj.d_N, obj.user_names);
        end
        
        function obj = loadUsersFromFile(obj, filename)
            [obj.user_names, obj.eff_N, obj.d_N] = loadUsersFromFile(filename);
        end
        
        %% Metric Calculation Methods
        function link_fluxes = calcLinkFluxes(obj)
            link_fluxes = obj.calcLinkFluxes_(obj.alloc_mat, obj.chan_fluxes);
            obj.link_fluxes = link_fluxes;
        end
        
        function user_fluxes = calcUserFluxes(obj)
            user_fluxes = obj.calcUserFluxes_(full(obj.inc_mat),obj.link_fluxes);
            obj.user_fluxes = user_fluxes;
        end
        
        function R_N = calcUserSingleRate(obj)
            R_N = obj.calcUserSingleRate_(obj.user_fluxes, obj.eff_N, obj.d_N);
            obj.R_N = R_N;
        end
        
        function A_L = calcLinkAccRate(obj)
            A_L = obj.calcLinkAccRate_(full(obj.inc_mat), obj.R_N, obj.coinc_wind_s); 
            obj.A_L = A_L;
        end
        
        function C_L = calcLinkCorrRate(obj)
            C_L = obj.calcLinkCorrRate_(obj.inc_mat, obj.link_fluxes, obj.eff_N);
            obj.C_L = C_L;
        end
        
        function lambda_L = calcStateWeights(obj)
            lambda_L = obj.calcStateWeights_(obj.A_L, obj.C_L);
            obj.lambda_L = lambda_L;
        end
        
        function rho_L = calcLinkStates(obj)
            rho_L = obj.calcLinkStates_(obj.alloc_mat, obj.link_fluxes,...
                                       obj.chan_fluxes, obj.U_K, obj.rho_K);
            obj.rho_L = rho_L;
        end
        
        function sigma_L = calcNoisyLinkStates(obj)
            D = size(obj.rho_L,1);
            n = size(obj.rho_L,3);
            noise_L = repmat(eye(D)./D, 1,1,n);
            
            sigma_L = obj.calcNoisyLinkStates_(obj.rho_L, noise_L, obj.lambda_L);
            obj.sigma_L = sigma_L;
        end
        
        function F_L = calcFidelity(obj, psi)
            F_L = obj.calcFidelity_(obj.sigma_L, psi);
            obj.F_L = F_L;
        end
        
        function ebr_L = calcEbitRates(obj)
            logneg_L = zeros(obj.L,1);
            for i = 1:obj.L
                logneg_L(i) = logneg(obj.sigma_L(:,:,i));
            end
            
            ebr_L = obj.calcEbitRate_(obj.A_L, obj.C_L, obj.sigma_L);
            obj.ebr_L = ebr_L;
        end
        
        function x = evalAlloc(obj, alloc_vect, ref_rho)
            %EVALALLOC returns QLAN performance metrics given a spectrum
            %allocation
            %
            %   [FID, EBR] = evalAlloc(QLAN, ALLOC_VECT, RHO) returns the
            %   link state fidelities FID with respect to state RHO as well
            %   as the link entangled bit rates EBR given a spectrum
            %   allocation ALLOC_VECT. RHO must have the same dimensions as
            %   the channel unitaries U_k and initial channel states rho_k.
            %   ALLOC_VECT is a Kx1 vector of link indices l where K is the
            %   number of channel pairs and 0 <= l <= L, the number of 
            %   links.
            
            valid_alloc = @(x) isequal(size(x(:)), [obj.K 1]) && ...
                max(alloc_vect) <= obj.L && min(alloc_vect) >= 0;
            
            p = inputParser;
            p.addRequired('obj', @(x) isa(x,'QLAN'));
            p.addRequired('alloc_vect', valid_alloc);
            p.addRequired('ref_rho', @(x) all(size(ref_rho) == size(obj.U_K, [1 2])));
            
            p.parse(obj, alloc_vect, ref_rho)
            obj = p.Results.obj;
            alloc_vect = p.Results.alloc_vect(:);
            ref_rho = p.Results.ref_rho;
            
            obj.alloc_vect = alloc_vect;
            obj.alloc_mat = QLAN.allocVect2Mat(alloc_vect, obj.L);
            obj.calcLinkFluxes();
            obj.calcUserFluxes();
            obj.calcUserSingleRate();
            obj.calcLinkAccRate();
            obj.calcLinkCorrRate();
            obj.calcStateWeights();
            obj.calcLinkStates();
            obj.calcNoisyLinkStates();
            
            ref_rho = repmat(ref_rho,1,1,obj.L);
            obj.calcFidelity(ref_rho);
            obj.calcEbitRates();
            x = [obj.F_L(:), obj.ebr_L(:)];
        end
        
        %% Plotting/Display Methods
        function plotGraph(obj)
            edgetable = table(obj.links, obj.link_labels,...
                              'VariableNames', {'EndNodes', 'EdgeNames'});
            nodetable = table(obj.user_names, 'VariableNames', {'Name'});
            g = graph(edgetable, nodetable);
            plot(g, 'EdgeLabel', g.Edges.EdgeNames)
        end
        
        function dispUsers(obj)
            disp(table(obj.user_names, obj.eff_N*100, obj.d_N,...
                'VariableNames', {'Name', 'Efficiency (%)', 'Dark Count Rate'}))
        end
        
        %% Utility Methods
        function ind = usrIDX(obj, name)
            if ischar(usr) || (isstring(usr) && numel(usr) == 1)
                ind = find(strcmp(obj.user_names, name));
            else
                error("Pass char vect or string as user name")
            end
            
        end
        
        function ind = linkIDX(obj, usr1, usr2)
            isnumber = @(x) isscalar(x) && isnumeric(x);
            istext = @(x) ischar(x) || (isstring(x) && numel(x) == 1);
            
            inds = zeros(2,1);
            if istext(usr1) && istext(usr2)
                inds = [find(strcmp(obj.user_names,usr1));
                        find(strcmp(obj.user_names,usr2))];
            elseif isnumber(usr1) && isnumber(usr2)
                inds = [usr1;
                        usr2];
            else
                error("Arguments must both be integers, char vectors, or strings")
            end
            
            for ind = 1:obj.L
                if all(ismember(inds, obj.links(ind,:)))
                    return
                end
            end
        end
        
        function deg = userDegree(obj, usr)
            if ischar(usr) || (isstring(usr) && numel(usr) == 1)
                usr_ind = find(strcmp(obj.user_names, usr));
            elseif isscalar(usr) && isnumeric(usr)
                usr_ind = usr;
            else
                error("Integer, char vector, or string required")
            end
            
            deg = sum(any(obj.links == usr_ind, 2));
        end
    end
end
