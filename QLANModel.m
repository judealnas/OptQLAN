classdef QLANModel
    %QLANMODEL encapsulates the symbolic functions defining the experimental
    %model of a quantum flex-grid local access network
    
    properties (Constant=true)
        mu = sym('mu',{'positive','real'});
        etaA = sym('etaA',{'positive','real'});
        etaB = sym('etaB',{'positive','real'});
        dA = sym('dA',{'positive','real'});
        dB = sym('dB',{'positive','real'});
        tau = sym('tau',{'positive','real'});
        D = sym('D',{'positive','integer'});
        
        RA  = QLANModel.etaA*QLANModel.mu/2+QLANModel.dA;
        RB  = QLANModel.etaB*QLANModel.mu/2+QLANModel.dB;
        A = 8*QLANModel.tau*QLANModel.RA*QLANModel.RB;
        C = QLANModel.etaA*QLANModel.etaB*QLANModel.mu;
        lambda = 1/(1+QLANModel.A/QLANModel.C);
    end
    
    methods
        function obj = QLANModel()            
        end
    end
    
    methods (Static=true)
        function f = subsParamVars(func, vals)
            sub_vars = QLANModel.getParamVars();
            if numel(vals) ~= numel(sub_vars)
                error("Insufficient values provided");
            elseif numel(vals) > numel(sub_vars)
                warning("Excessive substitution values provided")
                vals = vals(1:numel(sub_vars));
            end
            
            f = subs(func, sub_vars, vals);
        end
        
        function vars = getParamVars()
            vars = [QLANModel.etaA,...
                QLANModel.etaB,...
                QLANModel.dA,...
                QLANModel.dB,...
                QLANModel.tau,...
                QLANModel.D];
        end
        
        function F = fidUniformState()
            F = ((QLANModel.D-1)*QLANModel.lambda+1)/QLANModel.D;
        end
        
        function RE = ebrUniformState()
            RE = (QLANModel.A + QLANModel.C)*...
                log2(sqrt(QLANModel.D)*QLANModel.fidUniformState());
        end
    end
end

