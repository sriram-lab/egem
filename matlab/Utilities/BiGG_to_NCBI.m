function newmodel = BiGG_to_NCBI(model)
    model.genes = strrep(model.genes, '_AT', '.');
    newmodel = model;
end

