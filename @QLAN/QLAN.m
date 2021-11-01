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
        D; % Dimension of total quantum system; usually 4
    end
    
    %% Set-protected properties
    % Calculated quantities
    properties (SetAccess=private)
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
        F_Lmax; % Max possible link fidelities
        mu_L_Fmax; % Link fluxes for max fidelity
        ebr_L; % Per link ebit rates
        ebr_Lmax; % Max possible entangled bit rates
        mu_L_ebrmax; % Link fluxes for max EBR
        eval_mode; % evaluation method chosen for current state
    end
    
    properties (Constant)
        EVALMODE_ISO = 'iso';
        EVALMODE_SIMULT = 'simult';
        EVALMODE_BATCH = 'batch';
    end
    
    %% Static Methods
    methods (Static=true)
        M = allocVect2Mat(V,L)
        
        function obj = saveUsersToFile(obj, filename)
            saveUsersToFile(filename, obj.eff_N, obj.d_N, obj.user_names);
        end
        
        function obj = loadUsersFromFile(obj, filename)
            [obj.user_names, obj.eff_N, obj.d_N] = loadUsersFromFile(filename);
        end
    end
    
    %% Instance Methods
    methods
        %% Configuration Methods
        function obj = QLAN(varargin)
            % Check if arg is string vect or cell string array
            istextvector = @(x) isvector(x) && (isstring(x) || iscellstr(x));
            
            % Check if arg is numeric vector
            isnumcvector = @(x)isvector(x) && isnumeric(x) ;
            
            p = inputParser;
            p.addRequired('links', @(x)isempty(x) || (size(x,2) == 2 && isnumeric(x)));
            p.addRequired('eff_N', isnumcvector);
            p.addRequired('d_N', isnumcvector);
            p.addRequired('chan_fluxes', isnumcvector);
            p.addRequired('chan_unitaries', @(x)isnumeric(x));
            p.addRequired('chan_states', @(x)isnumeric(x));
            p.addRequired('coinc_wind_s', @(x) isnumeric(x) && isscalar(x));
            p.addParameter('UserNames',[], istextvector);
            p.addParameter('LinkNames',[], istextvector);
            
            % Parsing input
            p.parse(varargin{:});
            links = p.Results.links;
            eff_N = p.Results.eff_N;
            d_N = p.Results.d_N;
            chan_fluxes = p.Results.chan_fluxes;
            chan_unitaries = p.Results.chan_unitaries;
            chan_states = p.Results.chan_states;
            coinc_wind_s = p.Results.coinc_wind_s;
            user_names = p.Results.UserNames(:);
            link_names = p.Results.LinkNames(:);
            
            % Verify states and unitaries are compatible
            if ~isequal(size(chan_unitaries),size(chan_states))
                error("Channel unitary and state dimensions incompatible. Unitaries of size %dx%dx%d while states of size %dx%dx%d",...
                    size(chan_unitaries,1),size(chan_unitaries,2),size(chan_unitaries,3),...
                    size(chan_states,1), size(chan_states,2),size(chan_states,3))
            end
            
            % Verify channel quantites have compatible dimensions
            if any(numel(chan_fluxes) ~= [size(chan_unitaries,3) size(chan_states,3)])
                error("Channel quantity dimensions not compatible. %d channel fluxes, %d channel unitaries, and %d channel states provided.",...
                    numel(chan_fluxes),size(chan_unitaries,3),size(chan_states,3))
            end
            
            % Saveing quantities to object
            obj.eff_N = eff_N;
            obj.d_N = d_N;
            obj.chan_fluxes = chan_fluxes;
            obj.U_K = chan_unitaries;
            obj.rho_K = chan_states;
            obj.coinc_wind_s = coinc_wind_s;
            
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
            
            % Getting user names
            if isempty(user_names) % If user_names empty, use default names
                user_names = default_user_names;
            elseif numel(user_names) < N % if insufficient names, pad using default names
                warning("%d node names provided when expecting %d. Padding with default names",...
                    numel(user_names), N);
                user_names = [cellstr(user_names) ; default_user_names(numel(user_names)+1:end)];
                
            elseif numel(user_names) > N % If excess names, ignore excess
                warning("Too many user names specified. Extra names ignored");
                user_names = user_names(1:N);
            end
            obj.user_names = user_names; % Saving user names
            
            % Constructing links
            if isempty(link_names) % If no link names provided, use user names to generate names
                obj.configLinks(links, 'GenLabels', user_names);
            else
                obj.configLinks(links, link_names);
            end
            
            % Calculating theoretical maxima
            usr_noise_params = obj.coinc_wind_s.*obj.d_N./obj.eff_N;
            for i = 1:obj.L
                usrs = obj.links(i,:);
                model = QLANLinkModel(usr_noise_params(usrs),obj.D);
                
                [obj.F_Lmax, xF] = model.fidelityMaximum();
                obj.mu_L_Fmax = xF./obj.coinc_wind_s;
                
                [ebrmax, xE] = model.EBRMaximum();
                obj.ebr_Lmax = prod(obj.eff_N(usrs))/obj.coinc_wind_s*ebrmax;
                obj.mu_L_ebrmax = xE./obj.coinc_wind_s;
            end
        end
        
        function obj = configChans(obj, chan_fluxes, chan_unitaries, chan_states)
            %CONFIGCHANS(Q, chan_fluxes, chan_unitry, chan_states)
            %configures the following QLAN properties:
            %   chan_fluxes - channel biphoton flux vector
            %   chan_unitaries - (:,:,K) matrix of channel unitaries
            %   chan_states - (:,:,K) matrix of initial channel density matrices
            if ~(isvector(chan_fluxes) && isnumeric(chan_fluxes))
                error("Channel fluxes must be a numeric vector")
            end
            
            chan_fluxes = chan_fluxes(:);
            
            if any(numel(chan_fluxes) ~= [size(chan_unitaries,3), size(chan_states,3)])
                error("Dimension Mismatch: Unequal channels")
            end
            
            if any(size(chan_unitaries, [1 2]) ~= size(chan_states, [1 2]))
                error("Dimension Mismatch: unitaries and states are incompatible sizes");
            end
            
            obj.U_K = chan_unitaries;
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
            p.FunctionName = 'configLinks';
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
        
        function D = get.D(obj)
            D = size(obj.rho_K,1);
        end
        
        %% Metric Calculation Methods
        % Signatures of methods in other files
        link_fluxes = calcLinkFluxes(obj, alloc_mat, chan_fluxes)        
        user_fluxes = calcUserFluxes(obj, inc_mat, link_fluxes)
        R_N = calcUserSingleRate(obj, user_fluxes, eff_N, d_N)
        A_L = calcLinkAccRate(obj, inc_mat, R_N, coinc_wind_s)
        C_L = calcLinkCorrRate(obj, inc_mat, link_fluxes, eff_N)
        lambda_L = calcStateWeights(obj, A_L, C_L)
        rho_L = calcLinkStates(obj, alloc_mat, link_fluxes, chan_fluxes,...
            U_K, rho_K)
        sigma_L = calcNoisyLinkStates(obj, noise_L, rho_L, lambda_L)
        F_L = calcFidelity(obj, rho_L, psi_L)
        ebr_L = calcEbitRate(obj, A_L, C_L, sigma_L)
        evalAlloc(obj, alloc_vect, ref_rho, varargin)
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
