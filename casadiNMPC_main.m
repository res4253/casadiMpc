%% 
% EXCELの設定を確認 !!
% EXCELの設定を確認 !!
% EXCELの設定を確認 !!
%%
clear
addpath('/csmnt/macos1014/applications/cfw/share/matlab/casadi/:');
addpath('/home/tanabe/casadi-3/')

global subjectNum
for subjectNum = [53]% 53 54 55 56 58 59 60]

import casadi.* 

Path = ['./model/:',...
        './database/:'...
        './fig/:'];
addpath(Path);

magniModel = magni2007CL('mpc');
manModel = man2007CL('mpc');


%% 設定 %%%%%%%%%%%%%
Ts = readTable('mpcConfig','Ts');
% t_end = readTable('simulationConfig','t_end');
t_end = 3000;
N = readTable('mpcConfig','N');
% subjectNum = readTable('simulationConfig','subjectNum');
%%%%%%%%%%%%%%%%%%%%

mpcModel = magniModel;

nx = mpcModel.nx;
nu = mpcModel.nu;
nd = mpcModel.nd;
ny = mpcModel.ny;


%% make reference %%%%%%
options = odeset('MaxStep', 0.5, 'RelTol',1e-12);
% [t_refs,x_refs] = ode45(@(t,x)magniModel.dynamics(t,x,...
%                       magniModel.insulin.getInsulinInput(t),...
%                       [magniModel.meal.getMealInput(t);
%                       magniModel.stomachD.getDInput(t)]),0:Ts:t_end+N*Ts ,magniModel.x_0,options);   

[t_refs,x_refs] = ode45(@(t,x)manModel.dynamics(t,x,0,...
                      [manModel.meal.getMealInput(t);
                      manModel.stomachD.getDInput(t)]),0:Ts:t_end+N*Ts ,manModel.x_0,options);   

x_refs_block = x_refs;                  
x_refs = reshape(x_refs', 1, [])'; 

%% make dist %%%%%%
dists = [];
for t = 0:Ts:t_end+(N-1)*Ts
    meal = [magniModel.meal.getMealInput(t);...
            magniModel.stomachD.getDInput(t)];
    dists = [dists; meal];
end

%% mpc %%%%%%%%%%%%%%%%%

Path = [Path, './casadi/:'];
addpath(Path)

nmpc = casadiNMpcCL(mpcModel);

xs = [];
ts = [];
us = [];
Xopts = [];
Uopts = [];
Yopts = [];
cts = [];
x = magniModel.x_0; 


% dist = repmat([0;75000],N,1); % 予測に用いる外乱
ub_x = inf(nx*(N+1),1); % 状態の上限
lb_x = -inf(nx*(N+1),1); % 状態の下限
ub_u = 30000*ones(nu*N,1); % 入力の上限
lb_u = 0*ones(nu*N,1); % 入力の下限
u_ref = magniModel.u_0*ones(nu*N,1); % 入力の参照値
% x_ref = repmat(magniModel.x_0,N+1,1);

step = 0;
for t = 0:Ts:t_end
    disp(['t = ', num2str(t)])
 
    x_ref = x_refs((step*nx)+1 : (step*nx)+nx*(N+1)); % 各時刻における時変目標値を x_refs から取り出す
    dist = dists((step*nd)+1 : (step*nd)+nd*N);
   
    nmpc = nmpc.solveNLP(t, x, dist, x_ref, u_ref, ub_x, lb_x, ub_u, lb_u);
    
    [time, x_hist] = ode45(@(t,x)magniModel.dynamics(t,x,...
                                nmpc.u,...
                                [magniModel.meal.getMealInput(t);...
                                magniModel.stomachD.getDInput(t)]),...
                            t:t+Ts,...
                            x,...
                            odeset('MaxStep', 0.5));
                        
    x = x_hist(end,:)';
    step = step + 1;
    
    xs = [xs; x_hist];
    ts = [ts; time];
    us = [us, nmpc.u];
    Xopts = [Xopts, nmpc.Xopt]; 
    Uopts = [Uopts, nmpc.Uopt];
    Yopts = [Yopts, nmpc.Uopt];
    cts = [cts, nmpc.compute_t];
end


figure
plot(t_refs, x_refs_block(:,1)./magniModel.V_G);
hold on
plot(ts,xs(:,1)./magniModel.V_G);
yline(mpcModel.x_0(1)./magniModel.V_G)

for i = 1:12
    figure
    plot(t_refs,x_refs_block(:,i))
    hold on
    plot(ts,xs(:,i))
    yline(mpcModel.x_0(i))
end

t_grid = 0:Ts:t_end;
figure
stairs(t_grid,us)

rmpath(Path);

end






