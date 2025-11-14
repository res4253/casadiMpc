function parameter = readParam(subjectNum, parameterName)
%         global m_2;
% 
%         if subjectNum ~= 9001 && subjectNum ~= 9002
%             if strcmp(parameterName,'m_2')
%                 if isempty(m_2)
%                     error('パラメータ が定義されていません, それか, globalに変更した部分をなおしてください')
%                 end
%                 parameter = m_2;
%                 return
%             end
%         end

        value = readtable('T1DMparameter.xlsx');
        eval(['parameterArray = value.s',num2str(subjectNum),';']);
        parameterNameArray = value.Parameter; 
        paramIndex = find(strcmp(parameterName, parameterNameArray));
        parameter = parameterArray(paramIndex);
        if isempty(parameter)
            error([parameterName,'という名前のパラメータは存在しません'])
        end
end

