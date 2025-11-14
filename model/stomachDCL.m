classdef stomachDCL < timeDeterminedInputCL
    properties
        D, D_0, beforePhaseNum, glucoseAbsorption, q_0, stomachDArray, meal, 
    end
    methods
        function self = stomachDCL(glucoseAbsorption, meal)
            self@timeDeterminedInputCL();
            self.inputArray = self.makeInputArray('meal');
            self.glucoseAbsorption = glucoseAbsorption;
            self.q_0 = self.glucoseAbsorption.x_0;
            self.D_0 = 45000;
            self.meal = meal;
            self.stomachDArray = self.makeStomachDArray;
        end
        
        %
        function stomachDArray = makeStomachDArray(self, ~)
            inputArray = self.inputArray;
            q_0 = self.q_0;
            
            stomachDArray = [];
            for phaseNum = 1:size(inputArray,1)
                D = q_0(1)+q_0(2) + inputArray(phaseNum,3); % todo Dの下限値は設定すべきな気がする
 
                if phaseNum == size(inputArray,1)
                   stomachDArray = [stomachDArray; 
                                 inputArray(phaseNum,1), self.t_end-inputArray(phaseNum,1),D];                    
                    break
                else
                   stomachDArray = [stomachDArray; 
                                 inputArray(phaseNum,1), inputArray(phaseNum+1,1)-inputArray(phaseNum,1),D];
                end
                
                t1 = inputArray(phaseNum,1);
                t2 = inputArray(phaseNum+1,1);
                options = odeset('MaxStep', 0.5, 'RelTol',1e-12);
                [t, q] = ode45(@(t,q)self.glucoseAbsorption.dynamics(t,q,0,[self.meal.getMealInput(t),D]),[t1,t2],q_0, options);
                q_0 = q(end,:);
            end
        end
        
        % timeDetermindInputCLの getDeterminedInput をoveride
        % stomachD用の getDeterminedInput
        % inputArray(i, 2)で割らない
        function input = getDeterminedInput(self, t) 
            inputArray = self.stomachDArray;
            input = self.D_0;
            for i = 1:size(inputArray,1)
                if t >=  inputArray(i,1) && t < inputArray(i,1) + inputArray(i,2)
                    input = inputArray(i,3);
                    break
                end
            end
        end

        function D = getDInput(self, t)
            D = self.getDeterminedInput(t);
        end
            
    end
end