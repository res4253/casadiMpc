classdef magni2007withoutGlucoseAbsorptionCL < magni2007CL
    properties
        R_aArray
    end
    methods
        function self = magni2007withoutGlucoseAbsorptionCL(~)
            self@magni2007CL('mpc')

            glucoseRa = glucoseRaCL(self.subjectNum, self.meal, self.stomachD);
            self.R_aArray = glucoseRa.getR_aArray();

            self.nx = 9;
            self.nu = 1;
            self.nd = 1;
            self.ny = 1;
                       
            self.x_0 = [self.glucoseSystem.x_0
                        self.insulinSystem.x_0
                        self.insulinSignal.x_0
                        self.insulinInjection.x_0
                       ];
            self.y_0 = self.x_0(1)/self.V_G;
            self.u_0 = self.insulinInjection.ins_b;
            self.d_0 = 0;
        end

        function x_dot = dynamics(self, t, x, u, d)
            G_p = x(1);
            G_t = x(2);
            I_l = x(3);
            I_p = x(4);
            I_1 = x(5);
            I_d = x(6);
            X = x(7);
            S_1 = x(8);
            S_2 = x(9);
            
            insInput = u(1);
            Ra = d(1);

            I = self.insulinSystem.getI(I_p);   
            S = self.insulinInjection.getS([S_1;S_2]);

            G_dot = self.glucoseSystem.dynamics(t,[G_p;G_t],0,[Ra; X; I_d]);
            I_dot = self.insulinSystem.dynamics(t, [I_l;I_p], 0, S);
            Is_dot = self.insulinSignal.dynamics(t, [I_1; I_d; X], 0, I);
            S_dot = self.insulinInjection.dynamics(t, [S_1;S_2],insInput,0);

            x_dot = [G_dot; I_dot; Is_dot; S_dot];
        end    
        
        function y = output(self,~,x,~,~) 
             y = [x(1); x(2); x(3); x(4); x(5); x(6); x(7); x(8); x(9)];
%              y = x(1)/self.V_G;
        end        
    end
end