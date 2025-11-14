classdef casadiMpcCL
    % casadiでのMPC用の変数定義クラス
    % mpcModelオブジェクトをもちいてインスタンス
    % mpcModelに必要な属性 -> nx, nu, nd, ny, dynamics, output
    properties
        attribute, mpcModel, nx, nu, nd, ny, N, Ts, evalFunction, J, w, g, p, ...
            w0, lam_w, lam_g, lbg, ubg, ubu, lbu, u, dist, y_ref, u_ref, ubw, lbw, ubx, lbx,...
            Yopt, Uopt, Xopt, compute_t
    end
    methods
        function self = casadiMpcCL(mpcModel)
            self.mpcModel = mpcModel;
            self.nx = mpcModel.nx; % 各変数の次元
            self.nu = mpcModel.nu;
            self.nd = mpcModel.nd;
            self.ny = mpcModel.ny;
            self.N = readTable('mpcConfig','N');
            self.Ts = readTable('mpcConfig','Ts');
            self.evalFunction = evalFunctionCL(self.N);
        end
        
        function self = makeVariableAndConst(self)
            import casadi.* 
            
            N = self.N;
            Ts = self.Ts;
            nx = self.nx;
            nu = self.nu;
            nd = self.nd;
            ny = self.ny;
            w_x      = [];
            w_y      = [];
            w_u      = [];
            w_d      = [];
            g        = [];
            ubg      = [];
            lbg      = [];
            p        = [];
            
            % qp の場合は, 状態ベクトルも作成する(時変システムの場合に備えて)
            if strcmp(self.attribute,'qp') 
                Ad = SX.sym('Ad',[nx nx]);
                Bd = SX.sym('Bd',[nx nd+nu]);
                Ad_elem = reshape(Ad,nx*nx,1);
                Bd_elem = reshape(Bd,nx*(nd+nu),1);
            end
            
            % 状態変数, 入力変数, 外乱変数, 時変パラメータ変数をSX.Symで定義
            % モデルのダイナミクスを用いて拘束条件も定義する
            % 外乱や, 参照値は時変にできるようにここでは全て変数として定義
            % その後, casadiNlpCLのsolveNLPにて改めて指定
            
            % 時変パラメータtを作成 システムの時変パラメータは, 外乱に追加されたし
            t = SX.sym('t',1);   
            
            Xk = SX.sym('X_0', nx);               % 初期時刻の状態ベクトル x0
            w_x      = [w_x; Xk];                 % 最適化変数の配列に追加

 
            for k = 0:N-1
                t_f = t + k*Ts;
                
                Uk = SX.sym(['U_', num2str(k)], nu);       % 時間ステージ k の制御入力 uk を表す変数
                w_u  = [w_u; Uk];                          % Uk を外乱配列に追加 
                Dk = SX.sym(['D_', num2str(k)], nd);       % 時間ステージ k の外乱 dk を表す変数
                w_d  = [w_d; Dk];                          % dk を外乱配列に追加
                
                Yk = self.mpcModel.output(t,Xk,Uk,Dk);          % 出力変数
                w_y      = [w_y; Yk];     
                
                Xk1 = SX.sym(['X_', num2str(k+1)], nx); % 時間ステージ k+1 の状態を表す変数
                w_x      = [w_x; Xk1];                  % xk1 を最適化変数配列(w_x)に追加                                
                
                % 状態方程式による拘束条件   
                if strcmp(self.attribute, 'nonlinear')
                    Xk_next = self.getDynamics(t_f, Xk, Uk, Dk);
                elseif strcmp(self.attribute, 'qp')
                    Xk_next = self.getDynamics(t_f, Xk, Uk, Dk, Ad, Bd);
                end   

                g      = [g; Xk_next-Xk1];
                ubg    = [ubg; zeros(nx,1)];       % 等式制約は lower-bound と upper-bound を同じ値にすることで設定
                lbg    = [lbg; zeros(nx,1)];       % 等式制約は lower-bound と upper-bound を同じ値にすることで設定
 
                % その他の制約
                [formula, const] = casadiConstraint(t, Xk, Uk); % その他の制約があればここで指定
                if ~isempty(const)
                    g = [g; formula];
                    ubg = [ubg; const(:,2)];
                    lbg = [lbg; const(:,1)];
                end  
                Xk = Xk1;
            end
            
            % 最後のステップの出力変数
            % 直達項に入る入力と外乱は, 前のステップの入力と外乱をそのまま使っている.
            % 本当にそれでいいのかわからない(きっと良くない)ので, 保留
            % uの長さをN+1にするほうがいいのか？<-それでもよさそう
            Yk = self.mpcModel.output(t+N*Ts,Xk1,Uk,Dk);
            w_y      = [w_y; Yk];
            
            % 出力の参照値, 入力の参照値を変数として定義
            Y_ref = SX.sym('Y_ref',[(N+1)*ny,1]);
            U_ref = SX.sym('U_ref',[N*nu,1]);
            
            % 評価関数
            % w_y　と y_ref の次元を一致させること
            % w_y は　w_x　の関数, 最適化変数として動かすのはw_xとw_uである.
            J = self.evalFunction.getEvalFunction(w_y, w_u, Y_ref, U_ref);
            
            self.J = J;
            %　各変数ベクトルの整形, すべて列ベクトル
            self.w = vertcat(w_x,w_u);
            self.g = vertcat(g{:});
            self.p = vertcat(t, w_d, Y_ref, U_ref); 
            
            % adaptiveの場合は, 状態ベクトル変数も追加
            if strcmp(self.attribute,'qp') 
                self.p = vertcat(t, w_d, Y_ref, U_ref, Ad_elem, Bd_elem);
            end
                
            % 拘束条件式ベクトル
            self.ubg = ubg;
            self.lbg = lbg;
            
            % 他の必要な配列を作成
            self.w0 = zeros(length(self.w),1);
            self.lam_w = zeros(length(self.w),1);
            self.lam_g = zeros(length(self.g),1);
            self.ubx = Inf(length(w_x),1);
            self.lbx = -Inf(length(w_x),1);
            self.ubu = Inf(length(w_u),1);
            self.lbu = -Inf(length(w_u),1);
            self.dist = zeros(size(w_d,1),1);
            self.y_ref = zeros(length(w_y),1); % y_ref は w_y と次元一致
            self.u_ref = zeros(length(w_u),1);
        end
    end
    
    methods(Abstract = true)
        Xk_next = getDynamics(self, t_f, Xk, Uk, Dk);
    end
    
end

%% note
% 出力参照にした.　入力変数の長さをどうするか
