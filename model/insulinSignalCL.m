classdef insulinSignalCL
    properties
        k_i, p_2U, x_0, I_b, I_1b, I_db, X_b
    end
    methods
        function self = insulinSignalCL(subjectNum)
            self.k_i = readParam(subjectNum, 'k_i');
            self.p_2U = readParam(subjectNum, 'p_2U');
            self.I_b = readParam(subjectNum, 'I_b');

            self.I_1b = self.I_b;
            self.I_db = self.I_b;
            self.X_b = 0;

            self.x_0 = [self.I_1b; self.I_db; self.X_b];
        end

        function x_dot = dynamics(self, ~, x, ~, d)

            k_i = self.k_i;
            p_2U = self.p_2U;
            I_b = self.I_b;

            I_1 = x(1);
            I_d = x(2);
            X = x(3);
            I = d(1);

            x_dot = [ -k_i*(I_1-I)
                    -k_i*(I_d-I_1)
                    -p_2U*X+p_2U*(I-I_b)
                    ];
        end
    end
end
