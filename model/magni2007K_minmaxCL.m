% k_gut を, k_minに変更

classdef magni2007K_minmaxCL < magni2007CL
    methods
        function self = magni2007K_minmaxCL(~)
            self@magni2007CL('mpc')
            self.glucoseAbsorption = glucoseAbsorptionK_minmaxCL(self.subjectNum);
        end
    end
end