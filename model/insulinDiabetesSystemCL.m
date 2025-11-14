classdef insulinDiabetesSystemCL < insulinSystemCL
    properties
        ins_b_perMinKg
    end
    methods 
        function self = insulinDiabetesSystemCL(subjectNum)
            self@insulinSystemCL(subjectNum)

            m_3b = self.getM_3(0);
            self.I_lb = self.m_2*self.I_b*self.V_I/(self.m_1+m_3b);
            self.I_pb = self.I_b*self.V_I;
            self.x_0 = [self.I_lb; self.I_pb];

            self.ins_b_perMinKg = ...
                ((self.m_2+self.m_4)*self.I_b*self.V_I-self.m_1*self.m_2*self.I_b*self.V_I/(self.m_1+m_3b)); %pmol/min/kg
        end
        
        function x_dot = dynamics(self, t, x, u, d)
            m_1 = self.m_1;
            m_2 = self.m_2;
            m_4 = self.m_4;

            I_l = x(1);
            I_p = x(2);          
            S = d(1);
            
            m_3 = self.getM_3(S);
            
            x_dot = [-(m_1+m_3)*I_l+m_2*I_p
                    -(m_2+m_4)*I_p+m_1*I_l+S
                    ];
        end
        
        function m_3 = getM_3(self, ~)
            m_1 = self.m_1;
            HE_b = self.HE_b;
            m_3 = m_1*(HE_b/(1-HE_b));
        end
    end
end
