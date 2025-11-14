function array = readInputTable(sheet)
%     % todo
%     global period
%     if isempty(period)
%         error('perodが定義されていません, それか, globalに変更した部分をなおしてください')
%     end
%     %%%%
    

    value = readtable('nominalInput.xlsx','sheet',sheet);
    tableValue = table2array(value);
    array = tableValue(:,2:4);
    
%     %%%
%     array(:,2) = period*ones(length(array(:,2)),1);
%     %%%
end