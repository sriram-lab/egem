% Check metabolic_sensitivity module

fileList = dir('./../metabolic_sensitivity');
for file = 1:length(fileList)
    disp(fileList(file))
end
    checkcode fileList(file)
end