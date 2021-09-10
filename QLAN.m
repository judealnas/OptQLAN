classdef QLAN < handle
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
    
    
    properties
        alloc_vect; % Allocation vector
        g graph = graph(); % network graph object
        coinc_wind_s; % Coincidence window in seconds
        user_names; % network user labels
        link_labels; % link names
        inc_mat; % Network incidence matrix in sparse format
        eff_N (:,1); % Col vector of user detector efficiencies
        d_N (:,1)
        U_K  % DxDxK matrix of channel unitaries
        rho_K % DxDxK matrix of initial channel states
        chan_fluxes;
    end
    
    % Quanitites derived from other properties' set methods
    properties (Dependent=true)
        N; % # of users; calculated as size(eff_N,1);
        K; % # of channel pairs; size(chan_fluxes,1);
        L; % # of non-null network links; size(links,1);
    end
    % Calculated quantities
    properties (SetAccess=protected)
        links (:,2) {};
        user_file; % Filepath from which users were imported
        alloc_mat; % matrix expansion of allocation vector
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
    
    methods
        function obj = QLAN(varargin)
            
            % Check if arg is string or character vector
            istext = @(x) isstring(x) || ischar(x);
            
            % Check if arg is string or cell string array
            istextcolvector = @(x) size(x,2) == 1 && (isstring(x) || iscellstr(x));
            
            % Check if arg is numeric vector
            isnumcolvector = @(x) size(x,2) == 1 && isnumeric(x) ;
            
            p = inputParser;
            p.addRequired('links', @(x) size(x,2) == 2 && isnumeric(x));
            p.addRequired('file_eff', @(x) istext(x) || isnumcolvector(x));
            p.addOptional('d_N', [], isnumcolvector);
            p.addParameter('UserNames',[], istextcolvector);
            p.addParameter('LinkNames',[], istextcolvector);
            
            p.parse(varargin{:});
            links = p.Results.links;
            eff_N = p.Results.file_eff;
            d_N = p.Results.d_N;
            user_names = p.Results.UserNames;
            link_names = p.Results.LinkNames;
            
            % initialize empty graph
            g = graph();
            
            % Getting users from file if necessary
            if istext(eff_N)
                [file_user_names, eff_N, d_N] = loadUsersFromFile(eff_N);
                
                % If no user names provided, use names from file
                if isempty(user_names)
                    user_names = file_user_names;
                end
            end
            obj.eff_N = eff_N;
            obj.d_N = d_N;
            
            N = size(eff_N,1);
            default_user_names = arrayfun(@(x)sprintf('n%d',x), (1:N)',...
                'UniformOutput', false);
            if max(links, [], 'all') > N
                error("LINKS specifies more users than available.")
            end
            
            % If user_names empty, use default names 
            if isempty(user_names)
                user_names = default_user_names;
            
            % if insufficient names, pad using default names
            elseif numel(user_names) < N
                warning("Not enough node names provided. Padding with default names")
                user_names = [cellstr(user_names) default_user_names(numel(user_names)+1:end)];
            
            % If excess names, ignore excess
            elseif numel(user_names) > N
                warning("Too many user names specified. Extra names ignored");
                user_names = user_names(1:N);
            end
            obj.user_names = user_names;
%             % If link labels not provided, generate from user names
%             if isempty(link_names)
%                 link_names = generateLinkLabels(user_names, links);
%             end
            
            % Constructing graph
            g = g.addnode(user_names);
            obj.configLinks(links, 'GenLabels', user_names);
%             g = g.addedge(links(:,1), links(:,2));
%             
%             % Saving instance variables
%             obj.g = g;
%             obj.inc_mat = abs(incidence(obj.g));
%             obj.links = links;
%             obj.link_labels = link_names;
        end
        
        function obj = configChans(obj, chan_fluxes, U_K, rho_K)
            %CONFIGCHANS(Q, chan_fluxes, U_K, rho_K) configures the
            %following QLAN properties:
            %   chan_fluxes - channel biphoton flux
            %   U_K - (:,:,K) matrix of channel unitaries
            %   rho_K - (:,:,K) matrix of initial channel density matrices
            sz1 = [size(chan_fluxes,1), size(U_K,3), size(rho_K,3)];
            if any(sz1 ~= sz1(1))
                error("Dimension Mismatch: Unequal channels")
            end
            
            if any(size(U_K, [1 2]) ~= size(rho_K, [1 2]))
                error("Dimension Mismatch: unitaries and states of different sizes");
            end
            
            obj.U_K = U_K;
            obj.rho_K = rho_K;
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
            istextcolvector = @(x) size(x,2) == 1 && (isstring(x) || iscellstr(x));
            
            p = inputParser;
            p.FunctionName = 'configLinks';
            p.addRequired('obj', @(x) isa(x,'QLAN'));
            p.addRequired('links', @(x) size(x,2) == 2 && ...
                                   isnumeric(x) && ...
                                   max(x,[],'all') <= obj.N);
            p.addOptional('labels',[], istextcolvector);
            p.addParameter('GenLabels', [], istextcolvector);
            
            p.parse(obj,links,varargin{:});
            obj = p.Results.obj;
            links = p.Results.links;
            labels = p.Results.labels;
            user_names = p.Results.GenLabels;
            
            n = size(links,1); 
            
            % Generating default link labels
            if isempty(user_names)
                default_labels = arrayfun(@(x)sprintf('l%d',x), (1:n)',...
                    'UniformOutput', false);
            else
                default_labels = generateLinkLabels(user_names, links);
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
            
            edgetable = table(links, labels,...
                'VariableNames', {'EndNodes', 'EdgeNames'});
            
            obj.g = obj.g.rmedge(1:size(obj.g.Edges,1)); % delete all edges
            obj.g = obj.g.addedge(edgetable);
            
            obj.inc_mat = abs(incidence(obj.g));
            obj.links = links;
            obj.link_labels = labels;
            obj.links = links;
        end
        
        %%Dependent property accessors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function k = get.K(obj)
            k = numel(obj.chan_fluxes);
        end
        
        function n = get.N(obj)
            n = size(obj.eff_N, 1);
        end
        
        function l = get.L(obj)
            l = size(obj.g.Edges,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%User File loading/saving%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = saveUsersToFile(obj, filename)
            saveUsersToFile(filename, obj.eff_N, obj.d_N, obj.user_names);
        end
        
        function obj = loadUsersFromFile(obj, filename)
            [obj.user_names, obj.eff_N, obj.d_N] = loadUsersFromFile(filename);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function link_fluxes = calcLinkFluxes(obj)
            link_fluxes = obj.alloc_mat*obj.chan_fluxes;
            obj.link_fluxes = link_fluxes;
        end
        
        function user_fluxes = calcUserFluxes(obj)
            user_fluxes = full(obj.inc_mat)*obj.link_fluxes(2:end);
            obj.user_fluxes = user_fluxes;
        end
        
        function R_N = calcUserSingleRate(obj)
            R_N = obj.eff_N./2.*obj.user_fluxes + obj.d_N;
            obj.R_N = R_N;
        end
        
        function A_L = calcLinkAccRate(obj)
            A_L = zeros(obj.L, 1);
            full_inc_mat = logical(full(obj.inc_mat));
            for i = 1:obj.L
                A_L(i) = ...
                    8*obj.coinc_wind_s.*prod(obj.R_N(full_inc_mat(:,i)));
            end
            obj.A_L = A_L;
        end
        
        function C_L = calcLinkCorrRate(obj)
            C_L = zeros(obj.L, 1);
            full_inc_mat = logical(full(obj.inc_mat));
            for i = 1:obj.L
                C_L(i) = ...
                    obj.link_fluxes(i+1)*prod(obj.eff_N(full_inc_mat(:,i)));
            end
            obj.C_L = C_L;
        end
        
        function lambda_L = calcStateWeights(obj)
            lambda_L = 1./(1 + obj.A_L./obj.C_L);
            obj.lambda_L = lambda_L;
        end
        
        function rho_L = calcLinkStates(obj)
            rho_L = zeros([size(obj.rho_K,[1 2]) obj.L]);
            for i = 2:(obj.L+1) % iterate over non-null links
                alloc = obj.alloc_mat(i,:);
                U = obj.U_K(:,:,alloc); % unitaries of channels allocated to link i-1
                rho = obj.rho_K(:,:,alloc); % initial states of channels allocated to link i-1
                chan_flux = obj.chan_fluxes(alloc); % get channel fluxes allocated to link i-1
                chan_flux = reshape(chan_flux,1,1,numel(chan_flux)); % reshape to 3D vector for quicker math
                
                % calculate sum of all chan_flux*U*rho*U' density matrices
                % associated with link and normalize by total link flux
                rho_L(:,:,i-1) = sum(pagemtimes(chan_flux,pagemtimes(pagemtimes(U,rho),'none',U,'ctranspose')),3);
                rho_L(:,:,i-1) = rho_L(:,:,i-1)./obj.link_fluxes(i);
            end
            
            obj.rho_L = rho_L;
        end
        
        function sigma_L = calcMixedLinkStates(obj)
            D = size(obj.rho_L,1);
            n = size(obj.rho_L,3);
            noise = repmat(eye(D)./D, 1,1,n);
            
            lambda_L_3D = reshape(obj.lambda_L, 1, 1, obj.L);
            
            sigma_L = pagemtimes(lambda_L_3D, obj.rho_L)...
                + pagemtimes(1-lambda_L_3D, noise);
            obj.sigma_L = sigma_L;
        end
        
        function F_L = calcFidelity(obj, psi)
            if size(psi) ~= size(obj.sigma_L, [1 2])
                error("Reference density matrix psi must have same dimensions as 1st two dims of sigma_L to calculate link fidelities")
            end
            
            F_L = zeros(obj.L,1);
            for i = 1:obj.L
                F_L(i) = fidelity(obj.sigma_L(:,:,i), psi);
            end
            obj.F_L = F_L;
        end
        
        function ebr_L = calcEbitRates(obj)
            logneg_L = zeros(obj.L,1);
            for i = 1:obj.L
                logneg_L(i) = logneg(obj.sigma_L(:,:,i));
            end
            
            ebr_L = (obj.C_L + obj.A_L).*logneg_L;
            obj.ebr_L = ebr_L;
        end
        
        
        function x = evalAlloc(obj, alloc_vect, ref_rho)
            obj.alloc_vect = alloc_vect;
            obj.alloc_mat = logical((0:obj.L)' == alloc_vect(:)');
            obj.calcLinkFluxes();
            obj.calcUserFluxes();
            obj.calcUserSingleRate();
            obj.calcLinkAccRate();
            obj.calcLinkCorrRate();
            obj.calcStateWeights();
            obj.calcLinkStates();
            obj.calcMixedLinkStates();
            x = [obj.calcFidelity(ref_rho) obj.calcEbitRates()];
        end
        
        function plotGraph(obj)
            plot(obj.g, 'EdgeLabel', obj.g.Edges.EdgeNames)
        end
    end
end

