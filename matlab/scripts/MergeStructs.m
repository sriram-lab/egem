function merged_structure = MergeStructs(struct1, struct2)

merged_structure = struct1;
size1 = length(merged_structure);
fields2 = fieldnames(struct2);

for i=1:length(struct2)
    for j=1:length(fields2)
        merged_structure(size1+i).fields2({k}) = struct2.fields2({j});
    end
end

end