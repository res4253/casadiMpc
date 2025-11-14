classdef glucoseDiabetesSystemCL < glucoseSystemCL
    methods
        function self = glucoseDiabetesSystemCL(subjectNum)
            self@glucoseSystemCL(subjectNum)
        end
        % 抽象メソッドを実体化
        function EGP = getEGP(self, G_p, I_d, ~)
            k_p1 = self.k_p1;
            k_p2 = self.k_p2;
            k_p3 = self.k_p3;
            EGP = max(k_p1 - k_p2 * G_p - k_p3 * I_d, 0);
        end
    end
end
