function merged_structure = MergeStructs(struct1, struct2)

if isempty(struct1)
    merged_structure = struct2;
    return
end
if isempty(struct2)
    merged_structure = struct1;
    return
end

merged_structure = struct1;

size1 = length(fieldnames(merged_structure));
for i=1:length(struct2)
    fields2 = fieldnames(struct2);
    for j=1:length(fields2)
        merged_structure(size1+i).(fields2{j}) = struct2.(fields2{j});
    end
end

end