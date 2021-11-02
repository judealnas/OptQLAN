classdef QLANLinkModel
    %QLANLINKMODEL encapsulates the dimensionless equations modelling an
    %ideal, one-link network. Channel unitaries are assumed identity.
    %
    %   model = QLANLinkModel(Y_VECT,D) Returns a link model object.
    %       Y_VECT - a two-element vector containing the noise parameters
    %               of the link endpoint users. A user's noise parameter is
    %               defined as: y = (coincidence window)*(dark count
    %               rate)/(efficiency)
    %       D - dimension of the quantum system shared across the link; if
    %           unprovided, defaults to 4, the size of a two-qubit system
    %
    % -----------------------------------------------------------------
    %   Static Properties:
    %
    %   QLANLinkModel.general_fidelity_func - returns a function handle
    %       with signature
    %              fid(x,yA,yB,D)
    %       which calculates the fidelity of a link state shared between
    %       two users with noise parameters yA,yB with normalized flux x
    %       and state dimension D. Flux is normalized by the following
    %       equation:
    %               x = (coincidence window)*(biphoton flux)
    %
    %   QLANLinkModel.general_scaledEBR_func - returns a function handle
    %       with signature
    %           scaledEBR(x,yA,yB,D)
    %       which calculates the scaled EBR of a link receiving normalized
    %       flux of x with user noise parameters yA,yB and system dimension
    %       D. The scaled EBR can be converted to units of ebits/s by the
    %       following equation:
    %           EBR = [etaA*etaB/tau^2]*scaledEBR
    %       where etaN is the efficiency of user N and tau is the
    %       coincidence window in seconds
    %
    %   QLANLinkModel.isnzEBR(yA,yB,D) - returns a boolean value. If false,
    %       then a link with noise parameters yA, yB and state dimension D
    %       will have a maximum EBR of 0, i.e. no communication possible.
    
    properties
        y (1,2) double; % Vector of user noise parameters
        D = 4; % Dimension of shared link state
    end
    
    properties (Constant)
        % Fidelity as a function of normalized flux x, noise parameters yA and yB, and shared state dimension D
        general_fidelity_func = @(x,yA,yB,D)1/D*...
            (1+(D-1).*x./(x.^2+2.*(yA+yB+0.5).*x+4*yA.*yB)) 
        
        % Entangled bit rate as a function of normalized flux x, noise parameters yA and yB, and shared state dimension D

        general_scaledEBR_func = @(x,yA,yB,D)(x.^2+2.*(yA+yB+0.5).*x+4*yA.*yB).*...
            log2(sqrt(D).*QLANLinkModel.general_fidelity_func(x,yA,yB,D))
        
        general_DscaledEBr_func = @(x,yA,yB,D)...
            (2*(x+yA+yB)+1)*...
            log2(sqrt(D).*QLANLinkModel.general_fidelity_func(x,yA,yB,D))+...
            (D-1).*(4*yA.*yB-x.^2)./(log(2).*(x.^2+(2*(yA+yB)+D).*x+4*yA.*yB))
        
        general_bound_nzEBR = @(y,D) y-sqrt(QLANLinkModel.alpha(D).*y*2)+QLANLinkModel.alpha(D)./2;
        
        isEBRnz = @(yA,yB,D)(yA-yB).^2 - QLANLinkModel.alpha(D).*(yA+yB) + QLANLinkModel.alpha(D).^2/4 > 0
    end
    
    properties (Constant,Hidden)
        alpha = @(D) (D-sqrt(D))./(sqrt(D)-1);
    end
    
    methods
        function obj = QLANLinkModel(y,D)
            narginchk(1,2);
            if (numel(y) ~= 2) || any(y<0)
                error("y must be a two-element vector of non-negative real numbers")
            end
            
            obj.y = y(:);
            if nargin == 2
                obj.D = D;
            end
            
            if ~QLANLinkModel.isEBRnz(obj.y(1),obj.y(2),obj.D)
                warning("Noise parameters %.3f and %.3f result in a max EBR of 0",y(1),y(2))
            end
        end
        
        function func = boundNZEBR(obj)
            %%BOUNDNZEBR Returns a function handle for the boundary line
            %%marking where the maximum EBR is 0.
            %
            %   func = boundNZEBR(obj) returns a function handle for the
            %   boundary line marking where the maximum EBR is 0 with the
            %   value of D stored in OBJ substituted in.
            
            func = @(y)obj.general_bound_nzEBR(y,obj.D);
        end
        function func = fidelityFunc(obj)
            %FIDELITYFUNC Returns a function handle calculating link state
            %fidelity
            %
            %   func = fidelityFunc(obj) returns a function handle for link
            %   state fidelity with the noise parameters YA, YB and
            %   dimension D stored in OBJ substitued in.
            func = @(x)obj.general_fidelity_func(x,obj.y(1),obj.y(2),obj.D);
        end
        
        function func = scaledEBRFunc(obj)
            %%SCALEDEBRFUNC Returns a function handle for the scaled EBR.
            %
            %   func = scaledEBRFunc(OBJ) returns a function handle for the
            %   scaled EBR with the noise parameters YA, YB and dimension D
            %   substitued in.
            func = @(x)obj.general_scaledEBR_func(x,obj.y(1),obj.y(2),obj.D);
        end
        
        function func = scaledDiffEBRFunc(obj)
            %%SCALEDEBRFUNC Returns a function handle for the derivative of
            %%the scaled EBR with resepct to normalized flux
            %
            %   func = scaledEBRFunc(OBJ) returns a function handle for the
            %   derivative with respect to normalized flux of scaled EBR
            %   with the noise parameters YA, YB and dimension D substitued
            %   in.
            func = @(x)obj.general_DscaledEBr_func(x,obj.y(1),obj.y(2),obj.D);
        end
        
        function [Fmax, xFmax] = fidelityMaximum(obj)
            %%FIDELITYMAXIMUM Calculates the maximum possible fidelity of a
            %%link and the flux required to achieve that fidelity
            %
            %   Fmax = fidelityMaximum(OBJ) returns the maximum fidelity
            %   possible for a link described by the noise parameters and
            %   dimension size contained in OBJ
            %
            %   [Fmax, xFmax] = fidelityMaximum(obj) also returns the
            %   normalized link flux requried to achieve the maximum
            %   fidelity
            Dim = obj.D;
            Fmax = 1/Dim*(1+(Dim-1)./(4*sqrt(prod(obj.y))+2*sum(obj.y)+1));            
            xFmax = 2*sqrt(prod(obj.y));
        end
        
        function [scaledEBRmax, xEBRmax] = EBRMaximum(obj)
            %%EBRMAXIMUM calculates the maximum possible scaled entangled
            %%bit rate of a link and the flux required to achieve that rate
            %
            %   scaledEBRmax = EBRMaximum(OBJ) returns the maximum scaled
            %   EBR possible for a link described by the noise
            %   parameters and dimension size contained in OBJ
            %
            %   [Fmax, xFmax] = fidelityMaximum(obj) also returns the
            %   normalized link flux requried to achieve the maximum
            %   rate
            %
            % Scaled EBR can be converted to true EBR via:
            %       (true EBR) = (scaled EBR)*(tau)/(etaA*etaB)
            % where tau is the coincidence window and the etaN's are the
            % user detector efficiencies.
            yA = obj.y(1);
            yB = obj.y(2);
            
            % Check if non-zero EBR is possible; immediately return 0 if
            % not
            if ~obj.isEBRnz(yA,yB,obj.D)
                scaledEBRmax = 0;
                xEBRmax = NaN;
                return
            end
            
            % Otherwise, find
            ebr_func = obj.scaledEBRFunc();
            [xEBRmax, scaledEBRmax] = fminbnd(@(x)-ebr_func(x),0,2);
            scaledEBRmax = abs(scaledEBRmax);
        end
    end
end

