classdef QLAN < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        D = 4; % System dimensionality; D = 4 for 2-qubit system
        g graph; % network graph object
        user_names; % network user labels
        inc_mat; % Network incidence matrix in sparse format
        eff_N (:,1); % Col vector of user detector efficiencies
        d_N (:,1)
        U_K  % DxDxK matrix of channel unitaries
        rho_K % DxDxK matrix of initial channel states
        chan_fluxes;
        links (:,2);
    end
    
    % Quanitites derived from other properties' set methods
    properties (Dependent=true)
        N; % # of users; calculated as size(eff_N,1);
        K; % # of channel pairs; size(chan_fluxes,1);
        L; % # of non-null network links; size(links,1);
    end
    % Calculated quantities
    properties (SetAccess=protected)
        
        link_fluxes; % (L+1)x1 vector of link biphoton fluxes
        user_fluxes; % Nx1 vector of user single photon fluxes
        R_N; % Nx1 vector of single photon rates per detector per user
        A_L; % Accidental coincidence rate per link
        C_L; % Correlated coincidence rate per link
        F_L; % Per link fidelities
        EN_L; % Per link ebit rates

    end
    
    methods
        function obj = QLAN(varargin)
            % If string or char vect input, interpret as filepath and load
            % users from file
            if nargin == 1 && (isstring(varargin{1}) || ischar(varargin{1}))
                [obj.user_names, obj.eff_N, obj.d_N] = ...
                                            loadUsersFromFile(varargin{1});
            else 
                error("Currently only supports initialization from file");
            end
        end
        
        function obj = createLinks(links)
            obj.links = links;
            obj.L = size(links,2);
            obj.link_names = generateLinkLabels(obj.user_names, links);
            obj.g = graph(links(:,1),links(:,2));
        end
        
        function obj = set.chan_fluxes(obj,fluxes)
            obj.K = numel(fluxes);
        end
        
        function obj = saveUsersToFile(obj, filename)
            saveUsersToFile(filename, obj.eff_N, obj.d_N, obj.user_names);
        end
        
        function obj = loadUsersFromFile(obj, filename)
            [obj.user_names, obj.eff_N, obj.d_N] = loadUsersFromFile(filename);
        end
        
        function user_fluxes = calcUserFluxes(obj)
        end
        
        function link_fluxes = calcLinkFluxes(obj)
        end
        
        function R_N = calcUserSingleRate(obj)
        end
        
        function A_L = calcLinkAccRate(obj)
        end
        
        function C_L = calcLinkCorrRate(obj)
        end
        
        function rho_L = calcLinkStates(obj)
        end
        
        function sigma_L = calcMixedLinkStates(obj)
        end
           
        function F_L = calcFidelity()
        end
        
        function ebr_L = calcEbitRates()
        end
        
        
        function x = evalAlloc(obj, alloc_vect)
            alloc_mat = double((0:obj.L)' == alloc_vect(:)');
        end
    end
end

