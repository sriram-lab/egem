function new_model = addBiomassRxn(model)
    [vals, rxns] = read_txt('./../../../biomass_rxn.xlsx');
    stoichiometry = vals;
    bigg_id = rxns(2:end, 2);
    rxn_names = rxns(2:end, 3);
    
    
    