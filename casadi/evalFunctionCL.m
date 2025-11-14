classdef evalFunctionCL
    properties
        Q_stage, R_stage, Q_terminal, N
    end
    
    methods
        function self = evalFunctionCL(N)
            % Q_stage : ステージコスト, ny × 1
            % R_stage : ステージコスト, nu × 1
            % Q_terminal : ステージコスト, ny × 1
            Q_stage = readTable('mpcConfig','Q_stage');
            R_stage = readTable('mpcConfig','R_stage');
            Q_terminal = readTable('mpcConfig','Q_terminal');
            
            self.Q_stage = Q_stage;
            self.R_stage = R_stage;
            self.Q_terminal = Q_terminal;
            self.N = N;
        end
        
        function J = getEvalFunction(self, y, u, y_ref, u_ref)
            Q_stage = self.Q_stage;
            R_stage = self.R_stage;
            Q_tarminal = self.Q_terminal;
            N = self.N;

            Q_stage    = diag(Q_stage);        % 状態への重み(ステージコスト)
            R_stage    = diag(R_stage);         % 制御入力への重み(ステージコスト)
            Q_tarminal = diag(Q_tarminal);        % 状態への重み(終端コスト)
            
            
            %% Q,Rの行列を作成
            Q = [];
            R = [];
            for i = 0:N-1
                Q = blkdiag(Q,Q_stage);
            end
            Q = blkdiag(Q,Q_tarminal);
            
            for i = 0:N-1
                R = blkdiag(R,R_stage);
            end
            
            % 評価関数  
            J = 0.5 * (((y-y_ref)' * Q * (y-y_ref)) +  ((u-u_ref)' * R * (u-u_ref)));
        end
    end
end