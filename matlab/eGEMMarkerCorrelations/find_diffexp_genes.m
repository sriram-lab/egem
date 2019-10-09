function [diffExp_genes] = find_diffexp_genes(model, celllinenames_ccle1, locGeneExp)
    
    load ccle_geneExpression_vars;
    diffExp_genes.name = [];
    celllinenames_ccle1 = celllinenames_ccle1(locGeneExp);
    cell_gene_exp = ccle_expression_metz(:, locGeneExp);

    
    for cell = 1:length(celllinenames_ccle1)
        cellLine = celllinenames_ccle1(cell);
        cell_line_positions = find(ismember(string(celllinenames_ccle1), ...
            string(cellLine)));
        tmpExpression = cell_gene_exp(:, cell_line_positions);
        
        if length(cell_line_positions) < 2
            
            cell_line_ON = ccleids_met(tmpExpression >= 2);
            cell_line_OFF = ccleids_met(tmpExpression <= -2);
            
            model_ON = model.genes(ismember(model.genes, cell_line_ON));
            model_OFF = model.genes(ismember(model.genes, cell_line_OFF));

        elseif cell_line_positions > 2
            tmpExpression = ccle_expression_metz(:, cell_line_positions);
            tmpExpression = mean(tmpExpression, 2);
            
            cell_line_ON = ccleids_met(tmpExpression >= 2);
            cell_line_OFF = ccleids_met(tmpExpression <= -2);
            
            model_ON = model.genes(ismember(model.genes, cell_line_ON));
            model_OFF = model.genes(ismember(model.genes, cell_line_OFF));
        end
        
        ON_fieldname = string(strcat('ON_', cellLine));
        OFF_fieldname = string(strcat('OFF_', cellLine));
        
        diffExp_genes.(ON_fieldname) = model_ON;
        diffExp_genes.(OFF_fieldname) = model_OFF;
    
    end
end