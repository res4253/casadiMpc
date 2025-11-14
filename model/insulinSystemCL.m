classdef (Abstract = true) insulinSystemCL 
    properties
        m_1, m_2, m_4, HE_b, x_0, V_I, I_b, I_lb, I_pb
    end
    
    methods
        function self = insulinSystemCL(subjectNum)
            self.m_1 = readParam(subjectNum, 'm_1');
            self.m_2 = readParam(subjectNum, 'm_2');
            self.m_4 = readParam(subjectNum, 'm_4');
            self.HE_b = readParam(subjectNum, 'HE_b');
            self.V_I = readParam(subjectNum, 'V_I');
            self.I_b = readParam(subjectNum, 'I_b');        
        end

        function I = getI(self, I_p)
            I = I_p/self.V_I;
        end
    end
end