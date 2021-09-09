classdef QLAN < handle
    %QLAN Representation of flex-grid quantum LAN
    %   Q = QLAN(FILEPATH) constructs a fully connected network using the
    %   user properties contained in the JSON-encoded file FILEPATH
    %
    %   Q = QLAN(EFF_N, D_N) constructs a fully connected network using the
    %   vectors of user efficiencies EFF_N and dark count rates D_N
    %
    %   Q = QLAN(_,'UserNames', NAMES) associates names with users of the
    %   network where NAMES is a string vector or cell vector of character
    %   vectors corresponding with the order of specified users
    %
    %   Q = QLAN(_,'Links', LINKS) connects the users according to LINKS.
    %   LINKS can be:
    %       - An Lx2 matrix where rows specify the index of endpoint users
    %       forming the link
    %       - A table with the variable "EndNodes". Rows of this table are
    %       the user indices or user names forming the link. The individual
    %       links can be given labels by including in the table a cell
    %       vector of char vectors with the variable name "EdgeNames".
   
    
    properties
        alloc_vect; % Allocation vector
        g graph = graph(); % network graph object
        coinc_wind_s; % Coincidence window in seconds
        user_names; % network user labels
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
            string_or_char = @(x) isstring(x) || ischar(x);
            
            % Check if arg is string or cell string array
            cellstr_or_string = @(x) isstring(x) || iscellstr(x);
            
            % Check if arg is numeric vector
            num_vect = @(x) isvector(x) && isnumeric(x);
            
            % Check if argument is Nx2 matrix or table with the variable
            % 'EndNodes'
            mat_or_table = @(x) size(x,2) == 2 ||...
                (istable(x) && any(ismember(x.Properties.VariableNames,{'EndNodes'})));
            
            p = inputParser;
            p.addRequired('file_eff', @(x) string_or_char(x) || num_vect(x));
            p.addOptional('d_N', [], num_vect);
            p.addParameter('UserNames',[], cellstr_or_string);
            p.addParameter('Links',[], mat_or_table);
            
            p.parse(varargin{:});
            file_eff = p.Results.file_eff;
            d_N = p.Results.d_N;
            names = p.Results.UserNames;
            links = p.Results.Links;
            
            % If string or char vect input, interpret as filepath and load
            % users from file
            if nargin == 1 && string_or_char(varargin{1})
                [obj.user_names, obj.eff_N, obj.d_N] = ...
                                            loadUsersFromFile(varargin{1});
                                        
            % Else if two numeric vetors passed, interpret as corresponding
            % efficiency and dark count rate vectors
            elseif nargin == 2 && num_vect(varargin{1}) && num_vect(varargin{2})
                eff_N = varargin{1};
                d_N = varargin{2};
                if size(eff_N(:)) ~= size(d_N(:))
                    error("Efficiency and dark count vectors must have same size")
                end
                
                obj.eff_N = eff_N;
                obj.d_N = d_N;
            end
            
            % Handling link parameter
            links = p.Results.Links;
            if ~isempty(links)
                if isnumeric(links)
                    obj.links = links;
                    g = graph(links(:,1), links(:,2));
                else % Else links is table with variable 'EndNodes'
                    g = graph(links)
                end
                
            else
                links = nchoosek(1:obj.L, 2);
                obj.links = links;
                g = graph(links(:,1), links(:,2));
            end
        end
        
        function obj = configChans(obj, chan_fluxes, U_K, rho_K)
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
        
        function obj = createLinks(links, labels)
            obj.links = links;
            obj.L = size(links,2);
            obj.g = graph(links(:,1),links(:,2));
            obj.inc_mat = abs(incidence(obj.g));
            
            if nargin > 1
                obj.g.Edges = [obj.g.Edges table(labels, 'VariableNames', {'EdgeNames'})];
            end    
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
            if size(psi_L) ~= size(obj.sigma_L, [1 2])
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
            
            obj.ebr_L = (obj.C_L + obj.A_L).*logneg_L;
        end
        
        
        function x = evalAlloc(obj, alloc_vect)
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
            x = [obj.calcFidelity() obj.calcEbitRates()];
        end
    end
end

