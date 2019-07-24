%% read_json: function to read json files
function dict = read_json(fname)
    % Read in file
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid);
    
    % Decode .json
    val = jsondecode(str);
    
    % extract useful information
    field = fieldnames(val);
    values = struct2cell(val);
    dict = [field, values];
    
    %dict = dict';
    %dict = horzcat(dict{:});
    %dict = dict';