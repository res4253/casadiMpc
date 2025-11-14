classdef glucoseSystemCL
    properties
        F_cns, k_1, k_2, k_p1, k_p2, k_p3, k_e1, k_e2, ...
            K_m0, V_m0, K_mx, V_mx, EGP_b, G_b, V_G, x_0, I_b, ...
            G_tb, G_pb,
    end
    
    methods 
        function self = glucoseSystemCL(subjectNum)
            self.F_cns = readParam(subjectNum,'F_cns');
            self.k_1 = readParam(subjectNum,'k_1');
            self.k_2 = readParam(subjectNum,'k_2');
            self.k_p2 = readParam(subjectNum,'k_p2');
            self.k_p3 = readParam(subjectNum,'k_p3');
            self.k_e1 = readParam(subjectNum,'k_e1');
            self.k_e2 = readParam(subjectNum,'k_e2');
            self.K_m0 = readParam(subjectNum,'K_m0');
            self.K_mx = readParam(subjectNum,'K_mx');
            self.V_mx = readParam(subjectNum,'V_mx');
            self.EGP_b = readParam(subjectNum,'EGP_b');
            self.G_b = readParam(subjectNum,'G_b');
            self.I_b = readParam(subjectNum,'I_b');
            self.V_G = readParam(subjectNum,'V_G');
            self.k_p1 = self.EGP_b+self.k_p2*self.G_b*self.V_G+self.k_p3*self.I_b;
            self.G_tb = (self.F_cns-self.EGP_b+self.k_1*self.G_b*self.V_G)/self.k_2;
            self.V_m0 = (self.EGP_b-self.F_cns)*(self.K_m0+self.G_tb)/self.G_tb;
            self.G_pb = self.G_b*self.V_G;
            self.x_0 = [self.G_pb; self.G_tb];
        end
        
        function x_dot = dynamics(self,~, x, ~, d) 
            k_1 = self.k_1;
            k_2 = self.k_2;
            U_ii = self.F_cns;
            
            G_p = x(1);
            G_t = x(2);
            R_a = d(1);
            X   = d(2);
            I_d = d(3);
            
            I_po = [];
            
            if length(d) > 3
                I_po = d(4);
            end
            
            EGP = self.getEGP(G_p, I_d, I_po);
            U_id = self.getU_id(G_t, X);
            E = self.getE(G_p);
            
            x_dot = [EGP+R_a-U_ii-k_1*G_p+k_2*G_t-E
                    -U_id+k_1*G_p-k_2*G_t
                    ];
        end
            
        function E = getE(self, G_p)
            k_e2 = self.k_e2;
            k_e1 = self.k_e1;
            if G_p > k_e2
                E = k_e1*(G_p - k_e2);
            else
                E = 0;
            end
        end
        
        function U_id = getU_id(self, G_t, X)
            K_m0 = self.K_m0;
            K_mx = self.K_mx;
            V_m0 = self.V_m0;
            V_mx = self.V_mx;
            
            K_m = K_m0 + K_mx*X;
            V_m = V_m0 + V_mx*X;
            U_id =(V_m * G_t) / (K_m + G_t);
        end
        
        function G = getG(self, G_p)
            G = G_p/self.V_G;
        end
    end
    
    methods (Abstract = true)
        EGP = getEGP(self, G_p, I_d, I_po);
    end
end