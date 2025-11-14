classdef insulinInjectionCL
    properties
        k_a1, k_a2, k_d, BW, x_0, S_1b, S_2b, ins_b;
    end
    methods
        function self = insulinInjectionCL(subjectNum, ins_b_perMinKg)
            self.k_a1 = readParam(subjectNum, 'k_a1');
            self.k_a2 = readParam(subjectNum, 'k_a2');
            self.k_d = readParam(subjectNum, 'k_d');
            self.BW = readParam(subjectNum, 'BW');
            self.ins_b = ins_b_perMinKg*self.BW; %pmol/min
            self.S_1b = self.ins_b/(self.BW*(self.k_a1+self.k_d)); %pmol/kg
            self.S_2b = self.k_d*self.S_1b/self.k_a2;            
            self.x_0 = [self.S_1b; self.S_2b];
        end

        function x_dot = dynamics(self, ~, x, u, ~)

            S_1 = x(1);
            S_2 = x(2);
            ins = u(1);

            k_a1 = self.k_a1;
            k_a2 = self.k_a2;
            BW = self.BW;
            k_d = self.k_d;

            x_dot = [-(k_a1+k_d)*S_1+ins/BW
                    k_d*S_1-k_a2*S_2
                    ];
        end

        function S = getS(self, x)
            k_a1 = self.k_a1;
            k_a2 = self.k_a2;
            S_1 = x(1);
            S_2 = x(2);
            S = k_a1*S_1+k_a2*S_2;
        end
    end
end