function [A_eq, Bu_eq, Bd_eq, C_eq, Du_eq, Dd_eq] = linearization(nonlinearModel)
nx = nonlinearModel.nx;
nu = nonlinearModel.nu;
nd = nonlinearModel.nd;

x_0 = nonlinearModel.x_0;
u_0 = nonlinearModel.u_0;
d_0 = nonlinearModel.d_0;

syms x [nx, 1];
syms u [nu, 1];
syms d [nd, 1];

% Dは結局tanhの中でQsto/Dという形ででてくるので, 並行点ではQsto=0より, Dの偏微分は0
x_dot = nonlinearModel.dynamics(0, x, u, d);
y = nonlinearModel.output(0,x,u,d);

% jacobian
A = jacobian(x_dot, x);
Bu = jacobian(x_dot, u);
Bd = jacobian(x_dot, d);
C = jacobian(y, x);
Du = jacobian(y, u);
Dd = jacobian(y, d);

%%
% equilibrium point
A_eq = subs(A, [x;u;d], [x_0;u_0;d_0]);
Bu_eq = subs(Bu, [x;u;d], [x_0;u_0;d_0]);
Bd_eq = subs(Bd, [x;u;d], [x_0;u_0;d_0]);
C_eq = subs(C, [x;u;d], [x_0;u_0;d_0]);
Du_eq = subs(Du, [x;u;d], [x_0;u_0;d_0]);
Dd_eq = subs(Dd, [x;u;d], [x_0;u_0;d_0]);
end










