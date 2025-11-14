classdef linearModelCL
    properties
        attribute, subjectNum, nonlinearModel, meal, insulin, stomachD, ...
        Ac, Bc, Cc, Dc, Buc, Bdc, Duc, Ddc, nx, nu, nd, ny, y_0, x_0, u_0, d_0, V_G, R_aArray
    end
    methods 
        function self = linearModelCL(model)
            self.nonlinearModel = model;
            self = self.setABCD;
            self.meal = model.meal;
            self.insulin = model.insulin;
            self.stomachD = model.stomachD;
            self.subjectNum = model.subjectNum;
            self.attribute = model.attribute;
            self.y_0 = model.y_0;
            self.x_0 = model.x_0;
            self.u_0 = model.u_0;
            self.d_0 = model.d_0;
            self.nx = model.nx; %状態の次元
            self.nu = model.nu; %入力の次元 
            self.nd = model.nd; %外乱の次元
            self.ny = model.ny;
            self.V_G = model.V_G;
        end
                
        function self = setABCD(self)
            [Ac, Buc, Bdc, Cc, Duc, Ddc] = linearization(self.nonlinearModel);
            
            self.Ac = double(Ac);
            self.Buc = double(Buc);
            self.Bdc = double(Bdc);
            self.Bc = [self.Buc self.Bdc];
            self.Cc = double(Cc);
            self.Duc = double(Duc);
            self.Ddc = double(Ddc);
            self.Dc = [self.Duc self.Ddc];
        end  
        
        function x_dot = dynamics(self, ~, x, u, d)
            Ac = self.Ac;
            Buc = self.Buc;
            Bdc = self.Bdc;
                
            x_tilda = (x - self.x_0);
            u_tilda = (u - self.u_0);
            d_tilda = (d - self.d_0);
            
            x_dot = Ac*x_tilda + Buc*u_tilda + Bdc*d_tilda;
        end
        
        function y = output(self, t, x, u, d)
            nonlinearModel = self.nonlinearModel;
            y = nonlinearModel.output(t, x, u, d);
        end
    end
end

        
        