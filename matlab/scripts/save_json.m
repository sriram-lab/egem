function save_json(data, filename)

fid = fopen(strcat(filename, '.json'), 'w');

if fid == -1
    error('Cannot create JSON file'); 
end

fwrite(fid, data, 'char');
fclose(fid);

end
