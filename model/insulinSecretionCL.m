classdef insulinSecretionCL
    properties
        S_b, alpha, beta, gamma, K, h, x_0, I_pob, Y_b
        
    end
    methods
        function self = insulinSecretionCL(subjectNum)
            self.S_b = readParam(subjectNum, 'S_b');
            self.alpha = readParam(subjectNum, 'alpha');
            self.beta = readParam(subjectNum, 'beta');
            self.gamma = readParam(subjectNum, 'gamma');
            self.h = readParam(subjectNum, 'h');
            self.K = readParam(subjectNum, 'K');            
            self.I_pob = self.S_b/self.gamma;
            self.Y_b = 0;

            
            self.x_0 = [self.I_pob; self.Y_b];
        end
        
    
        function x_dot = dynamics(self,~,x,~,d)   
            S_b = self.S_b;
            alpha = self.alpha;
            beta = self.beta;
            gamma = self.gamma;
            h = self.h;
            K = self.K;
            
            I_po = x(1);
            Y = x(2);
            G_dot = d(1);
            G = d(2);
            
            if G_dot > 0
                S_po = Y + K*G_dot + S_b;
            else
                S_po = Y + S_b;
            end
            
            if beta*(G-h) >= -S_b
                YIF = -alpha*(Y-beta*(G-h));
            else
                YIF = -alpha*Y-alpha*S_b;
            end   
            x_dot =  [-gamma*I_po+S_po
                     YIF
                     ];
        end
            
        function S = getS(self,I_po)
            gamma = self.gamma;
            S = gamma*I_po;
        end  
    end
end
