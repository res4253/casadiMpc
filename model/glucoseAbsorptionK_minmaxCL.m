classdef glucoseAbsorptionK_minmaxCL < glucoseAbsorptionMpcCL
    methods 
        function self = glucoseAbsorptionK_minmaxCL(subjectNum)
            self@glucoseAbsorptionMpcCL(subjectNum)
        end
        function k_gut = getK_gut(self,~, ~)
            k_min = self.k_min;
            k_max = self.k_max;

            k_gut = (k_max-k_min)/2;
%             k_gut = (k_max+k_min)/2;
%             k_gut = k_max;
        end   
    end
end
