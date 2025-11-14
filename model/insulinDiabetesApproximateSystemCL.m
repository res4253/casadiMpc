classdef insulinDiabetesApproximateSystemCL < insulinDiabetesSystemCL
    properties
        m_7
    end
    methods 
        function self = insulinDiabetesApproximateSystemCL(subjectNum)
            self@insulinDiabetesSystemCL(subjectNum);
            self.x_0 = self.I_pb;
            self.m_7 = self.ins_b_perMinKg/(self.I_b*self.V_I);
        end

        function x_dot = dynamics(self,~,x,~,d)
            m_7 = self.m_7;

            I_p = x(1);    
            S = d(1);

            x_dot = -(m_7)*I_p+S;
        end
    end
end