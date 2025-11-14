function value = readTable(sheet, label)
%% todo
%     global subjectNum;
% 
%     if isempty(subjectNum)
%         error('subjectNum が定義されていません, それか, globalに変更した部分をなおしてください')
%     end
%     if strcmp(label,'subjectNum')
%         value = subjectNum;
%         return
%     end
%     
%     
%     global Ts;
% 
%     if isempty(Ts)
%         error('Tsが定義されていません, それか, globalに変更した部分をなおしてください')
%     end
%     if strcmp(label,'Ts')
%         value = Ts;
%         return
%     end
%%

    sheetValue = readtable('nominalInput.xlsx','sheet',sheet);
    eval(['value = sheetValue.',label,';']);
    if iscell(value)
        value = cell2mat(value);
    else
        findValue = find(~isnan(value));
        value = value(findValue);
    end
end