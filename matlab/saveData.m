function [dataStruct] = saveData(data, var, exp)
    fileName = './../../vars/histoneCorrResults.mat';
    varName = strcat(string(var), var, exp);
    if exist(fileName, 'file')
        save(fileName, varName, '-append')
    else
        save(fileName, varName)
    end
end