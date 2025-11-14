% 時間依存の入力の設定, 
% insulinArray = [[入力開始時刻, 入力時間, 入力量]; ... [];];
% 重複する場合は, 前の入力が優先される
% 入力時間等はサンプリング時刻に一致させること % todo
% 食事量の単位は[mg]
% inputの単位は[mg/min]

classdef mealCL < timeDeterminedInputCL
    properties
        d_0
    end
    methods
        function self = mealCL()
            self@timeDeterminedInputCL()
            self.d_0 = 0;
            self.inputArray = self.makeInputArray('meal');
        end
        
        function mealInput = getMealInput(self, t)
            determinedMealInput = self.getDeterminedInput(t);
            mealInput = determinedMealInput;
        end
    end
end
