classdef glucoseAbsorptionDistKgutMpcCL < glucoseAbsorptionCL
    methods
        function self = glucoseAbsorptionDistKgutMpcCL(subjectNum)
            self@glucoseAbsorptionCL(subjectNum);
        end
        function x_dot = mpcDynamics(self,~,x,~,dist)
            
            k_gri = self.k_gri;
            k_abs = self.k_abs;
            
            Q_sto1 = x(1);
            Q_sto2 = x(2);
            Q_gut = x(3);
            Q_sto = Q_sto1 + Q_sto2;
            
            d_in = dist(1);
            k_gut = dist(2);
            
            x_dot = [-k_gri*Q_sto1+d_in
                -k_gut*Q_sto2+k_gri*Q_sto1
                -k_abs*Q_gut+k_gut*Q_sto2
                ];
        end            
    end
end
