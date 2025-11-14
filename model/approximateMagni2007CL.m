classdef approximateMagni2007CL
    properties
        attribute, subjectNum, glucoseSystem, insulinSystem, insulinInjection, insulinSignal, ...
            stomachD, glucoseAbsorption, meal, insulin, y_0, x_0, d_0, u_0, nx, nu, nd, ny, V_G, CR,
    end
    methods 
        function self = approximateMagni2007CL(use)
            self.attribute = 'diabetes';
            self.subjectNum = readTable('simulationConfig','subjectNum');
            self.V_G = readParam(self.subjectNum, 'V_G');
            self.CR = readParam(self.subjectNum, 'CHOr');

            self.insulinSystem = insulinDiabetesApproximateSystemCL(self.subjectNum);
            self.insulinSignal = insulinSignalCL(self.subjectNum);
            self.insulinInjection = insulinInjectionCL(self.subjectNum, self.insulinSystem.ins_b_perMinKg);

            if strcmp(use,'mpc')
                self.glucoseSystem = glucoseDiabetesApproximateSystemMpcCL(self.subjectNum);
                self.glucoseAbsorption = approximateGlucoseAbsorptionCL(self.subjectNum);
            else
                error('normal is not exist with this model')
            end

            self.meal = mealCL;
            self.insulin = insulinCL(self.insulinInjection.ins_b, self.CR);

            self.nx = 9;
            self.nu = 1;
            self.nd = 1;
            self.ny = 1;
                       
            self.x_0 = [self.glucoseSystem.x_0
                        self.insulinSystem.x_0
                        self.insulinSignal.x_0
                        self.insulinInjection.x_0
                        self.glucoseAbsorption.x_0
                       ];
            self.y_0 = self.x_0(1)/self.V_G;
            self.u_0 = self.insulinInjection.ins_b;
            self.d_0 = self.meal.d_0;          
        end

        function x_dot = dynamics(self,t,x,u,d)
            G_p = x(1);
            G_t = x(2);
            I_p = x(3);
            I_1 = x(4);
            I_d = x(5);
            X = x(6);
            S_1 = x(7);
            S_2 = x(8);
            Q_gut = x(9);
            
            insInput = u(1);
            d_in = d(1);

            Ra = self.glucoseAbsorption.getR_a(Q_gut);
            I = self.insulinSystem.getI(I_p);   
            S = self.insulinInjection.getS([S_1;S_2]);

            G_dot = self.glucoseSystem.dynamics(t,[G_p;G_t],0,[Ra; X; I_d]);
            I_dot = self.insulinSystem.dynamics(t, I_p, 0, S);
            Is_dot = self.insulinSignal.dynamics(t, [I_1; I_d; X], 0, I);
            S_dot = self.insulinInjection.dynamics(t, [S_1;S_2],insInput,0);
            Q_dot = self.glucoseAbsorption.dynamics(t,Q_gut,0,d_in);

            x_dot = [G_dot; I_dot; Is_dot; S_dot; Q_dot];
        end
        
        % 出力変数, 次元を変えたらny, y_0もかえること
        function y = output(self,~,x,~,~) 
%              y = [x(1); x(2); x(3); x(4); x(5); x(6); x(7); x(8); x(9); x(10); x(11); x(12)];
             y = x(1)/self.V_G;
        end
    end
end