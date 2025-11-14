classdef casadiAdaptiveQpCL < casadiMpcCL
    properties
        nlp, solver,
    end
    methods
        function self = casadiAdaptiveQpCL(mpcModel)
           addpath('../casadi');
           import casadi.*
        
           self@casadiMpcCL(mpcModel); 
           self.attribute = 'adaptive';
           self.max_iter = readTable('mpcConfig','max_iter');
           self.mu_min = readTable('mpcConfig','mu_min');
           
           % システム行列の離散化
           self = self.discritize;
           
           % 変数定義
           self = self.makeVariableAndConst(); 

           % 非線形計画問題の定式化
           self.nlp = struct('f',self.J,'x',self.w,'g',self.g,'p',self.p);
           self.solver = nlpsol('solver', 'ipopt', self.nlp, ...
               struct('calc_lam_p', true, 'calc_lam_x', true, 'print_time', false, ...
               'ipopt', struct('max_iter', self.max_iter, 'mu_min', self.mu_min, 'warm_start_init_point', 'yes', 'print_level', 0, 'print_timing_statistics', 'no')));    
        end
        
        
        function Xk_next = getDynamics(self, ~, Xk, Uk, Dk, Ad, Bd)
            Xk_next = Ad*Xk + Bd*[Uk; Dk];
        end
        
        % 実際に解く関数, ここれ各状態制約や参照値を毎ステップ指定.
        % t : このステップの時刻 1×1
        % x0 : このステップにおける状態の初期値 nx × 1
        % dist : 決定外乱, nd*N × 1
        % y_ref : 出力参照, ny*(N+1) × 1
        % u_ref : 入力参照, nu*N × 1
        % lbx, ubx   : 状態制約, nx*(N+1) × 1
        % lbu, ubu   : 入力制約, nu*N × 1  
        % Ad, Bd : 適応したシステム行列 nx × nx, nx × nu
        function  self = solveAdaptiveQP(self, t, x0, dist, y_ref, u_ref, ubx, lbx, ubu, lbu, Ad, Bd)  
            
            Ad_elem = reshape(Ad,[],1);
            Bd_elem = reshape(Bd,[],1);
            
            % パラメータ配列を作成  
            if ~isempty(dist)
               self.dist = dist; 
            end            
            if ~isempty(y_ref)
               self.y_ref = y_ref; 
            end
            if ~isempty(u_ref)
               self.u_ref = u_ref; 
            end
            if ~isempty(ubx)
                self.ubx = ubx;
            end
            if ~isempty(lbx)
                self.lbx = lbx;
            end
            if ~isempty(ubu)
                self.ubu = ubu;
            end            
            if ~isempty(lbu)
                self.lbu = lbu;
            end
            
            self.ubw = [self.ubx; self.ubu];
            self.lbw = [self.lbx; self.lbu];
            self.p = [t; self.dist; self.y_ref; self.u_ref; Ad_elem; Bd_elem];

            % 初期状態についての制約を設定
            self.lbw(1:self.nx) = x0;
            self.ubw(1:self.nx) = x0;
            
            
            
            % primal variables (x) と dual variables（ラグランジュ乗数）の初期推定解も与えつつ solve（warm start）
            % 'x0': 推定状態変数解, N_wx+N_wu × 1　(上のx0とは違うので注意)
            % 'p' : 時変パラメータ, 1+N_wd+N_Yref+N_Uref × 1
            % 'lbx' : 状態制約　N_wx+N_wu × 1
            % 'ubl' : 状態制約　N_wx+N_wu × 1
            % 'lbg' : 制約式の上限　N_g × 1
            % 'ubg' : 制約式の下限　N_g × 1
            % 'lam_x0' : 状態変数の推定ラグランジュ乗数 N_wx+N_wu × 1
            % 'lam_g0' : 制約式の推定ラグランジュ乗数 N_g × 1
            
            sol = self.solver('x0', self.w0, 'p', self.p, 'lbx', self.lbw, 'ubx', self.ubw, 'lbg', self.lbg, 'ubg', self.ubg, 'lam_x0', self.lam_w, 'lam_g0', self.lam_g);
            % sol.x : 最適変数, N_wx+N_wu × 1
            
            % 今回の解を保存
            self.w0     = full(sol.x);     
            self.lam_w = full(sol.lam_x);
            self.lam_g = full(sol.lam_g);
            
            
            N = self.N;
            nx = self.nx;
            nu = self.nu;
            nd = self.nd;
            
            % 最適入力や最適変数などを取り出しておく
            % U : 最適入力, nu × 1;
            % Uopt : 最適入力系列 nu*N × 1
            % Xopt : 最適状態系列 nx*(N+1) × 1
            % Yopt : 最適出力系列 ny*(N+1) × 1
            self.u     = full(sol.x( nx*(N+1)+1 : nx*(N+1)+nu )); 
            self.Uopt = full(sol.x( nx*(N+1)+1 : end )); 
            self.Xopt = full(sol.x( 1 : nx*(N+1) )); 

            
            %%%% 出力変数の計算, 最適化変数に入れてしまうのは違うきがするので.
            Xopt_block = reshape(self.Xopt, nx, []);
            Uopt_block = reshape(self.Uopt, nu, []);
            D_block = reshape(self.dist, nd, []);
            
            Yopt = [];
            mpcModel = self.mpcModel;
            Ts = self.Ts;
            
            for k = 1:N
                y_step =  mpcModel.output(t+Ts*(k-1), Xopt_block(:,k), Uopt_block(:,k), D_block(:,k));
                Yopt = [Yopt; y_step];
            end
            % N+1の点ではUとDは一つ前のステップのものを使う
            k = N+1;
            y_step =  mpcModel.output(t+Ts*(k-1), Xopt_block(:,k), Uopt_block(:,k-1), D_block(:,k-1));
            Yopt = [Yopt; y_step];             
            self.Yopt = Yopt;
            %%%%
        end   
    end
end
