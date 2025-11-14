classdef glucoseAbsorptionMpcCL < glucoseAbsorptionCL
    methods
        function self = glucoseAbsorptionMpcCL(subjectNum)
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
            D = dist(2);

            k_gut = self.getK_gut(Q_sto, D);
            
            x_dot = [-k_gri*Q_sto1+d_in
                -k_gut*Q_sto2+k_gri*Q_sto1
                -k_abs*Q_gut+k_gut*Q_sto2
                ];
        end        
        
        function k_gut = getK_gut(self,Q_sto, D)
            k_min = self.k_min;
            k_max = self.k_max;
            b = self.b;
            d = self.d;
            
            alpha_gut = 5 / (2 * D * (1 - b));
            beta_gut = 5 / (2 * D * d); 
            k_gut = k_min + (k_max - k_min) * (tanh(alpha_gut * (Q_sto - b * D))-tanh(beta_gut * (Q_sto - d * D))+2)/2;
        end      
    end
end
