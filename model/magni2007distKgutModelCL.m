classdef magni2007distKgutModelCL < magni2007CL
    properties
        attribute, subjectNum, glucoseSystem, insulinSystem, insulinInjection, insulinSignal, ...
            stomachD, glucoseAbsorption, meal, insulin, y_0, x_0, d_0, u_0, nx, nu, nd, ny, V_G, CR,
    end
    methods 
        function self = magni2007distKgutModelCL(use)
            self@magni2007CL(use);
            self.glucoseAbsorption = glucoseAbsorptionDistKgutMpcCL(self.subjectNum);
        end

        function x_dot = dynamics(self,t,x,u,d)
            G_p = x(1);
            G_t = x(2);
            I_l = x(3);
            I_p = x(4);
            I_1 = x(5);
            I_d = x(6);
            X = x(7);
            S_1 = x(8);
            S_2 = x(9);
            Q_sto1 = x(10);
            Q_sto2 = x(11);
            Q_gut = x(12);
            
            insInput = u(1);
            d_in = d(1);
            k_gut = d(2);

            Ra = self.glucoseAbsorption.getR_a(Q_gut);
            I = self.insulinSystem.getI(I_p);   
            S = self.insulinInjection.getS([S_1;S_2]);

            G_dot = self.glucoseSystem.dynamics(t,[G_p;G_t],0,[Ra; X; I_d]);
            I_dot = self.insulinSystem.dynamics(t, [I_l;I_p], 0, S);
            Is_dot = self.insulinSignal.dynamics(t, [I_1; I_d; X], 0, I);
            S_dot = self.insulinInjection.dynamics(t, [S_1;S_2],insInput,0);
            Q_dot = self.glucoseAbsorption.dynamics(t,[Q_sto1; Q_sto2; Q_gut],0,[d_in; k_gut]);

            x_dot = [G_dot; I_dot; Is_dot; S_dot; Q_dot];
        end
    end
end