% 時間依存の入力の設定, 
% insulinArray = [[入力開始時刻, 入力時間, 入力量]; ... [];];
% 重複する場合は, 前の入力が優先される
% 入力時間等はサンプリング時刻に一致させること % todo
% インスリン量の単位は[pmol]
% inputの単位は[pmol/min]

classdef insulinCL < timeDeterminedInputCL
    properties
        ins_0
    end
    methods
        function self = insulinCL(basalInsulin, CR)
            self@timeDeterminedInputCL()
            self.ins_0 = basalInsulin;
            CRInsulinTF = readTable('inputConfig','CRInsulinTF');
            
            % CRInsulinTFがtrueの場合は, インスリン量は全てCR比によって決定される
            if CRInsulinTF == true
                mealInputArray = self.makeInputArray('meal');
                % CR (g/U) で割り, 1U = 6000 をかけ、pmol単位でinputArrayを作成
                %mealinputArrayが空なら, inputArrayも空にする
                if ~isempty(mealInputArray)   
                    inputArray = ...
                        [mealInputArray(:,[1 2]), mealInputArray(:,3).*6000./(1000*CR)];
                elseif isempty(mealInputArray)
                    inputArray = [];
                end
                self.inputArray = inputArray;
            else
                self.inputArray = self.makeInputArray('insulin');
            end
        end
        
        function insulinInput = getInsulinInput(self, t)
            determinedInsulinInput= self.getDeterminedInput(t);
            insulinInput = determinedInsulinInput + self.ins_0;
        end
    end
end
