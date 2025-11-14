classdef man2007CL
    properties
        attribute, subjectNum, glucoseSystem, insulinSystem, insulinSecretion, insulinSignal, ...
            stomachD, glucoseAbsorption, meal, insulin, y_0, x_0, d_0, u_0, nx, nu, nd, ny, V_G,        
    end
    methods
        function self = man2007CL(~)
            self.attribute = 'normal';
            self.subjectNum = 0; % normal パラメータはsubjectNum = 0
            self.V_G = readParam(self.subjectNum, 'V_G');
            
            self.insulinSystem = insulinHealthySystem(self.subjectNum);
            self.insulinSignal = insulinSignalCL(self.subjectNum);
            self.insulinSecretion = insulinSecretionCL(self.subjectNum);
            self.glucoseSystem = glucoseHealthySystemCL(self.subjectNum);
            self.glucoseAbsorption = glucoseAbsorptionCL(self.subjectNum);
            
            self.meal = mealCL;
            self.stomachD = stomachDCL(self.glucoseAbsorption, self.meal);
             
            self.nx = 12;
            self.nu = 0;
            self.nd = 2;
            self.ny = 1;
            
            self.x_0 = [self.glucoseSystem.x_0
                        self.insulinSystem.x_0
                        self.insulinSignal.x_0
                        self.insulinSecretion.x_0
                        self.glucoseAbsorption.x_0
                        ];
                    
            self.y_0 = self.x_0(1)/self.V_G;
            self.d_0 = [self.meal.d_0
                        self.stomachD.D_0
                       ];
        end
        
        function x_dot = dynamics(self, t, x, ~, d)
            G_p = x(1);
            G_t = x(2);
            I_l = x(3);
            I_p = x(4);
            I_1 = x(5);
            I_d = x(6);
            X = x(7);
            I_po = x(8);
            Y = x(9);
            Q_sto1 = x(10);
            Q_sto2 = x(11);
            Q_gut = x(12);    
            
            d_in = d(1);
            D = d(2);      
            
            Ra = self.glucoseAbsorption.getR_a(Q_gut);
            I = self.insulinSystem.getI(I_p); 
            G = self.glucoseSystem.getG(G_p);
            S = self.insulinSecretion.getS(I_po);
            
            G_dot = self.glucoseSystem.dynamics(t,[G_p; G_t], 0, [Ra; X; I_d; I_po]);
            I_dot = self.insulinSystem.dynamics(t,[I_l; I_p], 0, S);
            Is_dot = self.insulinSignal.dynamics(t, [I_1; I_d; X], 0, I);
            S_dot = self.insulinSecretion.dynamics(t, [I_po; Y], 0, [G_dot(1)/self.V_G; G]);
            Q_dot = self.glucoseAbsorption.dynamics(t,[Q_sto1; Q_sto2; Q_gut],0,[d_in; D]);
            
            x_dot = [G_dot; I_dot; Is_dot; S_dot; Q_dot];
        end
        
        function y = output(self,~,x,~,~)
            y = x(1)/self.V_G;
        end    
    end
end

            
            
            
            
            