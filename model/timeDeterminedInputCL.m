% 時間依存の入力の設定,
% inputArray = [[入力開始時刻, 入力時間, 入力量]; ... [];];
% 重複する場合は, 前の入力が優先される
% 入力時間等はサンプリング時刻に一致させること % todo

classdef timeDeterminedInputCL
    properties
        inputArray, t_end
    end
    methods
        function self = timeDeterminedInputCL()
            self.t_end = readTable('simulationConfig','t_end');
        end
        function input = getDeterminedInput(self, t)
            inputArray = self.inputArray;
            input = 0;
            for i = 1:size(inputArray,1)
                if t >=  inputArray(i,1) && t < inputArray(i,1) + inputArray(i,2)
                    input = inputArray(i,3) / inputArray(i,2);
                    break
                end
            end
        end

        % makeInputArray
        % kind : meal, insulin, (leptin, glucagon ...)
        % regularTF, irregularTF : irregularやregularのどちらを使うか

        function inputArray = makeInputArray(self, kind)

            regularTF = readTable('inputConfig','regularTF');
            irregularTF = readTable('inputConfig','irregularTF');

            regularInputArray = readInputTable(['regular',kind]);
            regularInputArray(:,1) = regularInputArray(:,1)*60;
            irregularInputArray = readInputTable(['irregular',kind]);
            irregularInputArray(:,1) = irregularInputArray(:,1)*60;
            
            
            % インプットが0, の場合ははじく．
            if ~isempty(regularInputArray)
            zeroInputIndex = find(regularInputArray(:,3)==0);
            regularInputArray(zeroInputIndex,:) = [];
            end
            
            if ~isempty(irregularInputArray)
            zeroInputIndex = find(irregularInputArray(:,3)==0);
            irregularInputArray(zeroInputIndex,:) = [];
            end
            
            inputArray = [];
            for j = 1:size(regularInputArray,1)
                simulationDays = floor(self.t_end/(60*24));
                last_t_end = self.t_end-simulationDays*(60*24);

                for i = 1:simulationDays
                    inputArray = [inputArray;
                        regularInputArray(j,1)+24*60*(i-1), regularInputArray(j,2:3)];
                end

                if regularInputArray(j,1) < last_t_end
                    inputArray = [inputArray;
                        regularInputArray(j,1)+24*60*simulationDays, regularInputArray(j,2:3)];
                    
                end
            end
            
            index = find(irregularInputArray(:,1) > self.t_end);
            irregularInputArray(index,:) = [];

            if irregularTF == true && regularTF == true
                inputArray = [inputArray; irregularInputArray];
            elseif irregularTF == false && regularTF == true
                inputArray = [inputArray];
            elseif irregularTF == true && regularTF == false
                inputArray = [irregularInputArray];
            else
                error('regularTFかirregularTFのどちらかはtrueにしてください')
            end
            inputArray = sortrows(inputArray);

            for i = 1:size(inputArray,1)-1
                if inputArray(i,1)+inputArray(i,2) > inputArray(i+1,1)
                    error('入力の範囲が重なっているようです.それぞれ独立にしてください')
                else
                    if inputArray(end,1)+inputArray(end,2) > self.t_end
                        inputArray(end,2) = self.t_end - inputArray(end,1);
                    end
                end
            end
        end
    end
end
