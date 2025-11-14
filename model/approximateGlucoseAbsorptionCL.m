classdef approximateGlucoseAbsorptionCL < glucoseAbsorptionMpcCL
    methods
        function self = approximateGlucoseAbsorptionCL(subjectNum)
            self@glucoseAbsorptionMpcCL(subjectNum)

%             k_min = self.k_min;
%             k_max = self.k_max;
            
            self.k_abs = readParam(subjectNum, 'k');
            
%             self.k_abs =(k_max*0.0004)/k_min + k_min - 0.0022;
            self.x_0 = 0;
        end

        function x_dot = dynamics(self, ~, x, ~, d)
            k_abs = self.k_abs;
            %% todo
%             global k_abs
            
            Q_gut = x(1);
            d_in = d(1);
            
            x_dot = -k_abs*Q_gut+d_in;
        end
        
        function R_a = getR_a(self, Q_gut)
%             global k_abs % todo
%             
            BW = self.BW;
            f = self.f;
            k_abs = self.k_abs;
            R_a = f .* k_abs .* Q_gut ./ BW;
        end
    end
end
