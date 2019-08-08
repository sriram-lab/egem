%% @author: Scott Campit 
function var = read_txt(path, txt)
% read_txt takes in a path and a file name containing new lines for each
% element and outputs a matlab variable of either character vector or num
% vector. 
var = textread([path, txt], '%s', 'whitespace', '\n');
end
