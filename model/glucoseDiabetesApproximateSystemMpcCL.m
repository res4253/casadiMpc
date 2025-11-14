classdef glucoseDiabetesApproximateSystemMpcCL < glucoseDiabetesSystemMpcCL
    properties 
        MMchangeRate, MMintercept
    end
    methods
        function self = glucoseDiabetesApproximateSystemMpcCL(subjectNum)
            self@glucoseDiabetesSystemMpcCL(subjectNum);

            [changeRate, intercept] = self.diterminMichaelis; %todo
            
            self.MMchangeRate = changeRate;
            self.MMintercept = intercept;
        end

        function U_id = getU_id(self, G_t, X)
            V_mx = self.V_mx;
            V_m0 = self.V_m0;
            K_m0 = self.K_m0;
            G_tb = self.G_tb;
            MMchangeRate = self.MMchangeRate;

            V_m = V_m0 + V_mx*X;
            U_idb =(V_m * G_tb) / (K_m0 + G_tb);
            U_id = U_idb + (MMchangeRate(1)*X + MMchangeRate(2))*(G_t - G_tb);
        end

        function E = getE(~,~)
            E = 0;
        end

        function [changeRate, intercept] = diterminMichaelis(self)
            V_m0 = self.V_m0; 
            V_mx = self.V_mx;
            K_m0 = self.K_m0;
            G_b = self.G_b;
            G_tb = self.G_tb;
            
            syms G_t_ X_ changeRate_
            
            assume(G_t_, 'positive')
            assume(X_, 'positive')
            assume(changeRate_, 'positive')
            
            % F_cns = 1;
            V_m = V_m0 + V_mx*X_;
            % G_tb = (F_cns-EGP_b+k_1*G_b*V_G)/k_2;
            U_idb =(V_m * G_tb) / (K_m0 + G_tb);
            
            U_id_ = (V_m0 + V_mx*X_)*G_t_ / (K_m0 + G_t_);
            approxiLine_ = U_idb + changeRate_ * (G_t_ - G_tb);
            
            J_ = (U_id_ - approxiLine_)^2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            intFrom = 20 * G_tb/G_b;
            intTo = 250* G_tb/G_b;
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            squareSum_ = int(J_,G_t_,intFrom,intTo);
            diffSquareSum_ = diff(squareSum_, changeRate_);
            
            changeRate_ = solve(diffSquareSum_ == 0, changeRate_);
            
            intercept_ = (-(changeRate_*G_tb)) + U_idb;
            
            %todo%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %不明な誤差が発生, おそらくsymbolic式と数値との間で
            changeRate = sym2poly(changeRate_);
            intercept = sym2poly(intercept_);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
