function [diffExp_genes] = find_diffexp_genes(model, listOfCellLines)
    load ./../../vars/ccle_geneExpression_vars

    diffExp_genes.name = [];
    
    for cell = 1:length(listOfCellLines)
        cellLine = listOfCellLines(cell);
        cell_line_positions = find(contains(celllinenames_ccle1, cellLine));
        
        if length(cell_line_positions) < 2
            cell_gene_exp = ccle_expression_metz(:, cell_line_positions);
            
            cell_line_ON = ccleids_met(cell_gene_exp >= 2);
            cell_line_OFF = ccleids_met(cell_gene_exp <= -2);
            
            model_ON = model.genes(ismember(model.genes, cell_line_ON));
            model_OFF = model.genes(ismember(model.genes, cell_line_OFF));
        
        elseif cell_line_positions > 2
            cell_gene_exp = ccle_expression_metz(:, cell_line_positions);
            cell_gene_exp = mean(cell_gene_exp, 2);
            
            cell_line_ON = ccleids_met(cell_gene_exp >= 2);
            cell_line_OFF = ccleids_met(cell_gene_exp <= -2);
            
            model_ON = model.genes(ismember(model.genes, cell_line_ON));
            model_OFF = model.genes(ismember(model.genes, cell_line_OFF));
        end
        
        ON_fieldname = string(strcat('ON_', cellLine));
        OFF_fieldname = string(strcat('OFF_', cellLine));
        
        diffExp_genes.(ON_fieldname) = model_ON;
        diffExp_genes.(OFF_fieldname) = model_OFF;
    
    end
end