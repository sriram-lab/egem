function var = read_txt(txtfile)

nlines = 897;
var = cell(nlines, 1);

fileID = fopen('ccle_names.txt');
for idx = 1:nlines
    var(idx) = {fgetl(fileID)};
end
fclose(fileID);

end

h3_relval = table2array(h3relval)
h3_marks = table2array(h3marks)
h3_media = table2array(cclemedia)