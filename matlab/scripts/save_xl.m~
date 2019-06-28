function save_xl(filename, rowname, columnname, struct, epsilon2, exp)
% save_xl will save the matlab structure into an excel sheet. Each sheet
% name is the epsilon value. 

% INPUT:
    % For SRA and FBA:
        % row = medium components (50x1)
        % column = reactions (1x20)
        % body = data (50x20)
            % Growth rate
            % Metabolic flux
            % Shadow price
            % Reduced cost
    % For FVA:
        % row = medium components (50x1)
        % column = reactions (1x20)
        % body = data (50x20)
            % Excess min/max flux
            % Depletion min/max flux

switch exp
    case ('sra' | 'fba')
        %% Growth rate
        excess_grate = struct.excess_grate;
        depletion_grate = struct.depletion_grate;

        % Excess
        writecell({'Excess Growth Rate'}, filename, 'Sheet', epsilon2, 'Range', 'A1');
        writematrix(rowname, filename, 'Sheet', epsilon2, 'Range', 'B1:U1');
        writematrix(columnname, filename, 'Sheet', epsilon2, 'Range', 'A2:A51');
        writematrix(excess_grate, filename, 'Sheet', epsilon2, 'Range', 'B2:U51');

        % Depletion
        writecell({'Depletion Growth Rate'}, filename, 'Sheet', epsilon2, 'Range', 'W1');
        writematrix(rowname, filename, 'Sheet', epsilon2, 'Range', 'X1:AQ1');
        writematrix(columnname, filename, 'Sheet', epsilon2, 'Range', 'W2:W51');
        writematrix(depletion_grate, filename, 'Sheet', epsilon2, 'Range', 'X2:AQ51');

        %% Metabolic flux
        excess_flux = struct.excess_flux;
        depletion_flux = struct.depletion_flux;
        
        % Excess
        writecell({'Excess Flux'}, filename, 'Sheet', epsilon2, 'Range', 'A53');
        writematrix(rowname, filename, 'Sheet', epsilon2, 'Range', 'B53:U53');
        writematrix(columnname, filename, 'Sheet', epsilon2, 'Range', 'A54:A103');
        writematrix(excess_flux, filename, 'Sheet', epsilon2, 'Range', 'B53:U103');

        % Depletion
        writecell({'Depletion Flux'}, filename, 'Sheet', epsilon2, 'Range', 'W53');
        writematrix(rowname, filename, 'Sheet', epsilon2, 'Range', 'X53:AQ53');
        writematrix(columnname, filename, 'Sheet', epsilon2, 'Range', 'W54:W103');
        writematrix(depletion_flux, filename, 'Sheet', epsilon2, 'Range', 'X53:AQ103');
        
        %% Shadow price
        excess_sp = struct.excess_flux_sp;
        depletion_sp = struct.depletion_flux_sp;
        
        % Excess
        writecell({'Excess Shadow Price'}, filename, 'Sheet', epsilon2, 'Range', 'A105');
        writematrix(rowname, filename, 'Sheet', epsilon2, 'Range', 'B105:U105');
        writematrix(columnname, filename, 'Sheet', epsilon2, 'Range', 'A106:A155');
        writematrix(excess_sp, filename, 'Sheet', epsilon2, 'Range', 'B106:U155');

        % Depletion
        writecell({'Depletion Shadow Price'}, filename, 'Sheet', epsilon2, 'Range', 'W105');
        writematrix(rowname, filename, 'Sheet', epsilon2, 'Range', 'X105:AQ105');
        writematrix(columnname, filename, 'Sheet', epsilon2, 'Range', 'W106:W155');
        writematrix(depletion_sp, filename, 'Sheet', epsilon2, 'Range', 'X106:AQ155');

        %% Reduced cost
        excess_rc = struct.excess_flux_rc;
        depletion_rc = struct.depletion_flux_rc;
        
        % Excess
        writecell({'Excess Reduced Cost'}, filename, 'Sheet', epsilon2, 'Range', 'A157');
        writematrix(rowname, filename, 'Sheet', epsilon2, 'Range', 'B157:U157');
        writematrix(columnname, filename, 'Sheet', epsilon2, 'Range', 'A158:A207');
        writematrix(excess_rc, filename, 'Sheet', epsilon2, 'Range', 'B158:U207');

        % Depletion
        writecell({'Depletion Reduced Cost'}, filename, 'Sheet', epsilon2, 'Range', 'W157');
        writematrix(rowname, filename, 'Sheet', epsilon2, 'Range', 'X157:AQ157');
        writematrix(columnname, filename, 'Sheet', epsilon2, 'Range', 'W158:W207');
        writematrix(depletion_rc, filename, 'Sheet', epsilon2, 'Range', 'X158:AQ207');
    
    case 'fva'
        %% Max Flux
        excess_maxflux = struct.excess_maxflux;
        depletion_maxflux = struct.depletion_maxflux;

        % Excess
        writecell({'Excess Max Flux'}, filename, 'Sheet', epsilon2, 'Range', 'A1');
        writematrix(rowname, filename, 'Sheet', epsilon2, 'Range', 'B1:U1');
        writematrix(columnname, filename, 'Sheet', epsilon2, 'Range', 'A2:A51');
        writematrix(excess_maxflux, filename, 'Sheet', epsilon2, 'Range', 'B2:U51');

        % Depletion
        writecell({'Depletion Max Flux'}, filename, 'Sheet', epsilon2, 'Range', 'W1');
        writematrix(rowname, filename, 'Sheet', epsilon2, 'Range', 'X1:AQ1');
        writematrix(columnname, filename, 'Sheet', epsilon2, 'Range', 'W2:W51');
        writematrix(depletion_maxflux, filename, 'Sheet', epsilon2, 'Range', 'X2:AQ51');

        %% Metabolic flux
        excess_minflux = struct.excess_minflux;
        depletion_minflux = struct.depletion_minflux;
        
        % Excess
        writecell({'Excess Min Flux'}, filename, 'Sheet', epsilon2, 'Range', 'A53');
        writematrix(rowname, filename, 'Sheet', epsilon2, 'Range', 'B53:U53');
        writematrix(columnname, filename, 'Sheet', epsilon2, 'Range', 'A54:A103');
        writematrix(excess_minflux, filename, 'Sheet', epsilon2, 'Range', 'B53:U103');

        % Depletion
        writecell({'Depletion Min Flux'}, filename, 'Sheet', epsilon2, 'Range', 'W53');
        writematrix(rowname, filename, 'Sheet', epsilon2, 'Range', 'X53:AQ53');
        writematrix(columnname, filename, 'Sheet', epsilon2, 'Range', 'W54:W103');
        writematrix(depletion_minflux, filename, 'Sheet', epsilon2, 'Range', 'X53:AQ103');
end
end    


