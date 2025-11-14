classdef glucoseAbsorptionK_minCL < glucoseAbsorptionMpcCL
    methods 
        function self = glucoseAbsorptionK_minCL(subjectNum)
            self@glucoseAbsorptionMpcCL(subjectNum)
        end
        function k_gut = getK_gut(self,~, ~)
            k_min = self.k_min;
            k_max = self.k_max;

            k_gut = k_min;
        end   
    end
end
