%% @Author: Scott Campit
function A = densityplot(fil)
% densityplot.m displays the demand reactions and excess or depletion of
% medium as a function of epsilon2. The values were obtained from excel
% files generated in make_heatmap.m

% INPUTS:
    % xlsx: the Excel file containing the values we want to plot
    
% OUTPUTS:
    % density plots
    
%% densityplot.m
filename = ['C:\Users\scampit\Desktop\MeGEM\matlab\tables\', fil];
[~, sheets] = xlsfinfo(filename);

for sheet = 2:length(sheets)
    % Metabolic flux or growth rate
    excess_flux = xlsread(filename, sheets(sheet), 'B2:U51');
    depletion_flux = xlsread(filename, sheets(sheet), 'X2:AQ51');
    
    % Reduced costs (log)
    excess_redcost = xlsread(filename, sheets(sheet), 'B54:U103');
    excess_redcost = log(excess_redcost);
    depletion_redcost = xlsread(filename, sheets(sheet), 'X54:AQ103');
    depletion_redcost = log(depletion_redcost);
    
    % Shadow price (log)
    excess_shadow = xlsread(filename, sheets(sheet), 'B106:U155');
    excess_shadow = log(excess_shadow);
    depletion_shadow = xlsread(filename, sheets(sheet), 'X106:AQ155');
    depletion_shadow = log(depletion_shadow);
    
    eflux(:,:,sheet) = excess_flux;
    dflux(:,:,sheet) = depletion_flux;
end
    
end