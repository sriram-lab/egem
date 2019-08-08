function environmental_perturbation(model, condition)
    switch condition 
        case 'hypoxic'
            [~, pos] = ismember({'EX_o2(e)'}, model.rxns);
            model.lb(pos) = 0;
            model.ub(pos) = 0;
            disp('Hypoxic model')
        case 'normoxic'
            disp('Normoxic model')
    end
end