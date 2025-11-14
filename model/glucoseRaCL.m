classdef glucoseRaCL < glucoseAbsorptionCL
    properties
        stomachD, meal, 
    end

    methods 
        function self = glucoseRaCL(subjectNum, meal, stomachD)
            self@glucoseAbsorptionCL(subjectNum)
            self.stomachD = stomachD;
            self.meal = meal;
        end

        function R_aArray = getR_aArray(self)
            t_end = readTable('simulationConfig','t_end');
            N = readTable('mpcConfig','N');
            Ts = readTable('mpcConfig','Ts');
            meal = self.meal;
            stomachD = self.stomachD;

            options = odeset('MaxStep', 0.5, 'RelTol',1e-12);
            [t, q] = ode45(@(t,x)self.dynamics(t,x,0,[meal.getMealInput(t); stomachD.getDInput(t)]), 0:1:(t_end+Ts*N), self.x_0, options);
            
            
            %縦が状態次元, 横が時間
            q = q';
            
            R_aArray = self.getR_a(q(3,:));
        end
    end
end

