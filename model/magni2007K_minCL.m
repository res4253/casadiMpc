% k_gut を, k_minに変更

classdef magni2007K_minCL < magni2007CL
    methods
        function self = magni2007K_minCL(~)
            self@magni2007CL('mpc')
            self.glucoseAbsorption = glucoseAbsorptionK_minCL(self.subjectNum);
        end
    end
end