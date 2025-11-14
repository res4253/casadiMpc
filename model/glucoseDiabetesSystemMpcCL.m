classdef glucoseDiabetesSystemMpcCL < glucoseDiabetesSystemCL
    methods
        function self = glucoseDiabetesSystemMpcCL(subjectNum)        
            self@glucoseDiabetesSystemCL(subjectNum)
        end
        
        % 不連続な関数を近似してオーバーライド, getEはgrantParentのmethod
        function E = getE(self, G_p)
            k_e2 = self.k_e2;
            k_e1 = self.k_e1;
            
            E = k_e1*(G_p - k_e2)/(1+exp(-(G_p-k_e2)));
        end
        
        function EGP = getEGP(self, G_p, I_d, ~)
            k_p1 = self.k_p1;
            k_p2 = self.k_p2;
            k_p3 = self.k_p3;
            EGP = k_p1 - k_p2 * G_p - k_p3 * I_d;
        end
    end
end