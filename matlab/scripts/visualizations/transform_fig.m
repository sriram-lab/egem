%% Convert all fig files to tif
function transform_fig(path, old_ext, new_ext)

%path = './../figures/new-model/';
str = strcat(path, "*", old_ext);
files = dir(str);
for fil = 1:length(files)
    file = files(fil);
    old_filename = file.name;
    new_filename = regexprep(old_filename, old_ext, new_ext);
    %new_filename = regexprep(old_filename, '.fig', '.tif');
    tmp = openfig([path, old_filename]);
    set(tmp,...
        'Units', 'Inches',...
        'Position', [0, 0, 15, 10],...
        'PaperUnits', 'Inches',...
        'PaperSize', [8, 11],...
        'Resize', 'on');
    saveas(tmp, strcat(path, new_filename));
end
    
end