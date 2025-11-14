classdef glucoseHealthySystemCL < glucoseSystemCL
    properties 
        k_p4, I_pob
    end
    methods
        function self = glucoseHealthySystemCL(subjectNum)
            self@glucoseSystemCL(subjectNum)
            self.k_p4 = readParam(subjectNum,'k_p4');
            self.I_pob = readParam(subjectNum, 'I_pob');
            self.k_p1 = self.EGP_b+self.k_p2*self.G_b*self.V_G+self.k_p3*self.I_b+self.k_p4*self.I_pob;

            self.G_pb = self.G_b*self.V_G;
            self.x_0 = [self.G_pb; self.G_tb];
        end
        
        function EGP = getEGP(self, G_p, I_d, I_po)
            k_p1 = self.k_p1;
            k_p2 = self.k_p2;
            k_p3 = self.k_p3;
            k_p4 = self.k_p4;
            EGP = max(k_p1 - k_p2 * G_p - k_p3 * I_d - k_p4*I_po, 0);            
        end 
    end
end
