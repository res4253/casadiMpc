% k_gut を, k_minに変更

classdef magni2007K_maxCL < magni2007CL
    methods
        function self = magni2007K_maxCL(~)
            self@magni2007CL('mpc')
            self.glucoseAbsorption = glucoseAbsorptionK_maxCL(self.subjectNum);
        end
    end
end