%% @Author: Scott Campit
function densityplot(fil)
% densityplot.m displays the demand reactions and excess or depletion of
% medium as a function of epsilon2. The values were obtained from excel
% files generated in make_heatmap.m

% INPUTS:
    % xlsx: the Excel file containing the values we want to plot
    
% OUTPUTS:
    % density plots
    
%% densityplot.m
filename = ['./../tables/', fil];
[~, sheets] = xlsfinfo(filename);
for sheet = 1:length(sheets)
    % Metabolic flux or growth rate
    excess_flux = xlsread(filename, sheet, 'B2:U51');
    depletion_flux = xlsread(filename, sheet, 'X2:AQ51');
    
    % Reduced costs (log)
    excess_redcost = xlsread(filename, sheet, 'B54:U103');
    excess_redcost = log(excess_redcost);
    depletion_redcost = xlsread(filename, sheet, 'X54:AQ103');
    depletion_redcost = log(depletion_redcost);
    
    % Shadow price (log)
    excess_shadow = xlsread(filename, sheet, 'B106:U155');
    excess_shadow = log(excess_shadow);
    depletion_shadow = xlsread(filename, sheet, 'X106:AQ155');
    depletion_shadow = log(depletion_shadow);
    
    % test
    surf(excess_flux, sheet)
    
    
    
    
    
    
end