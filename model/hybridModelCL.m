classdef hybridModelCL
    properties
        insulin, meal, stomachD, x_0, dynamics1, dynamics2, sw1, sw2
    end
    
    methods
        function self = hybridModelCL(model1, model2, sw1, sw2)
            self.dynamics1 = @model1.dynamics1;
            self.dynamics2 = @model2.dynamics2;
            self.x_0 = model1.x_0;
            self.meal = model1.meal;
            self.stomachD = model1.stomachD;
            self.insulin = model1.insulin;
            self.sw1 = sw1;
            self.sw2 = sw2;
        end

        function x_dot = dynamics(t,x,u,d)
            sw1 = self.sw1;
            sw2 = self.sw2;
            dynamics1 = self.dynamics1;
            dynamics2 = self.synamics2;
            if x(10)+x(11) > sw2 || x(10)+x(11) < sw1
                x_dot = dynamics1(t,x,u,d);
            elseif x(10)+x(11) >= sw1 && x(10)+x(11) <= sw2
                x_dot = dynamics2(t,x,u,d);
            end
        end
    end
end

            
    