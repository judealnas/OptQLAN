classdef QLANModel
    %QLANMODEL contains the symbolic functions modelling entanglement
    %distribution between two users of a QLAN neglecting channel effects
    
    properties (Constant=true)
        sym_mu = sym('mu',{'positive','real'});
        sym_etaA = sym('etaA',{'positive','real'});
        sym_etaB = sym('etaB',{'positive','real'});
        sym_dA = sym('dA',{'positive','real'});
        sym_dB = sym('dB',{'positive','real'});
        sym_tau = sym('tau',{'positive','real'});
        sym_D = sym('D',{'positive','integer'});
        
        RA  = QLANModel.sym_etaA*QLANModel.sym_mu/2+QLANModel.sym_dA;
        RB  = QLANModel.sym_etaB*QLANModel.sym_mu/2+QLANModel.sym_dB;
        A = 8*QLANModel.sym_tau*QLANModel.RA*QLANModel.RB;
        C = QLANModel.sym_etaA*QLANModel.sym_etaB*QLANModel.sym_mu;
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
            vars = [QLANModel.sym_etaA,...
                QLANModel.sym_etaB,...
                QLANModel.sym_dA,...
                QLANModel.sym_dB,...
                QLANModel.sym_tau,...
                QLANModel.sym_D];
        end
        
        function F = fidUniformState()
            F = ((QLANModel.sym_D-1)*QLANModel.lambda+1)/QLANModel.sym_D;
        end
        
        function RE = ebrUniformState()
            RE = (QLANModel.A + QLANModel.C)*...
                log2(sqrt(QLANModel.sym_D)*QLANModel.fidUniformState());
        end
    end
    
    methods
        function testMethod(obj, arg1, arg2)
            narginchk(1,3)
            if nargin == 1 && isa(obj,'QLANModel')
                disp('Good')
                arg1 = 9;
                arg2 = 0;
            end
            
            disp(arg1 + arg2);
        end
    end
end

