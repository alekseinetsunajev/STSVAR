
classdef data < specification
    properties (SetAccess = private)
        Y
        T
        Z
        Theta0
        U0
        B0
        
    end
    properties
        options
        Lambda0
    end
    
    methods
        function obj = data()
                    obj.Y = xlsread('Y_CPI_NonEnergy_RealS_FF197001-200706.xls');
 %                   obj.Y = xlsread('Y_CPI_NonEnergy_RealS_FF198301-200706.xls');
                    obj.T = size(obj.Y);
                    obj.Z = calc_Z(obj.Y, obj.lags);
                    [obj.Theta0, obj.U0, ~, ~] = Var(obj.Y', obj.lags);
                    obj.U0 = obj.U0';
                    obj.B0 =  0.5*(1/obj.T(1,1)*(obj.U0'*obj.U0))^0.5 + 1/100000*randn(obj.T(1,2),obj.T(1,2));
                    %obj.Lambda0 = obj.LM *eye(obj.T(1,2));
        end
    end
end