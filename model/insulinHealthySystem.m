classdef insulinHealthySystem < insulinSystemCL
    properties
        m_5, m_6, 
    end
    methods
        function self = insulinHealthySystem(subjectNum)
            self@insulinSystemCL(subjectNum)
            self.m_5 = readParam(subjectNum, 'm_5'); 
            self.m_6 = readParam(subjectNum, 'm_6');
            
            self.I_pb = self.I_b*self.V_I;
            self.I_lb = (self.m_2+self.m_4)*self.I_pb/self.m_1;
            
            self.x_0 = [self.I_lb; self.I_pb];
        end
        
        function x_dot = dynamics(self, ~, x, ~, d)
            
            I_l = x(1);
            I_p = x(2);
            S = d(1);
            
            m_1 = self.m_1;
            m_2 = self.m_2;
            m_4 = self.m_4;
            m_5 = self.m_5;
            m_6 = self.m_6;
            
            HE = -m_5 * S + m_6;
            m_3 = HE * m_1/(1-HE);  
            
            
            x_dot = [-(m_1+m_3)*I_l+m_2*I_p+S
                    -(m_2+m_4)*I_p+m_1*I_l
                    ];
        end
    end
end