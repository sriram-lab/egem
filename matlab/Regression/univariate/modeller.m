function table = modeller(data, type, experiment)
array = table2array(data);
array = meanCenterPredictors(array);
    
switch type
        
        case 'CCLE'
                
                switch experiment
                        
                        case 'all'
                                X0 = array(:, [1:201, 244:247]);
                                Y = array(:, 202:243);

                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end

                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(data.Properties.VariableNames(202:243))', ...
                                        rho, pvalue, R2, RMSE];

                        case 'geneExpOnly'
                                X0 = array(:, 1:201);
                                Y = array(:, 202:243);

                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end
                                
                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(data.Properties.VariableNames(202:243))', ...
                                        rho, pvalue, R2, RMSE];

                        case 'fluxOnly'
                                X0 = array(:, 244:247);
                                Y = array(:, 202:243);

                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end

                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(data.Properties.VariableNames(202:243))', ...
                                        rho, pvalue, R2, RMSE];

                        case 'Me1'
                                markers = string(data.Properties.VariableNames);
                                array = table2cell(data);
                                Me1 = array(:, ~cellfun(@isempty, regexp(markers, 'me1')));
                                Me1Markers = data.Properties.VariableNames(~cellfun(@isempty, ...
                                        regexp(markers, 'me1')));
                                X0 = cell2mat(array(:, [1:201, 244:247]));
                                Y = cell2mat(Me1);

                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end

                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(Me1Markers)', ...
                                        rho, pvalue, R2, RMSE];

                        case 'Me2'
                                markers = string(data.Properties.VariableNames);
                                array = table2cell(data);
                                Me2 = array(:, ~cellfun(@isempty, regexp(markers, 'me2')));
                                Me2Markers = data.Properties.VariableNames(~cellfun(@isempty, ...
                                        regexp(markers, 'me2')));
                                X0 = cell2mat(array(:, [1:201, 244:247]));
                                Y = cell2mat(Me2);
                                
                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end
                                
                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(Me2Markers)', ...
                                        rho, pvalue, R2, RMSE];
                                
                        case 'Me3'
                                markers = string(data.Properties.VariableNames);
                                array = table2cell(data);
                                Me3 = array(:, ~cellfun(@isempty, regexp(markers, 'me3')));
                                Me3Markers = data.Properties.VariableNames(~cellfun(@isempty, ...
                                        regexp(markers, 'me3')));
                                X0 = cell2mat(array(:, [1:201, 244:247]));
                                Y = cell2mat(Me3);
                                
                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end
                                
                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(Me3Markers)', ...
                                        rho, pvalue, R2, RMSE];
                                
                        case 'Ac'
                                markers = string(data.Properties.VariableNames);
                                array = table2cell(data);
                                Ac = array(:, ~cellfun(@isempty, regexp(markers, 'ac')));
                                AcMarkers = data.Properties.VariableNames(~cellfun(@isempty, ...
                                        regexp(markers, 'ac')));
                                X0 = cell2mat(array(:, [1:201, 244:247]));
                                Y = cell2mat(Ac);
                                
                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end
                                
                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(AcMarkers)', ...
                                        rho, pvalue, R2, RMSE];
                end
                
        case 'LeRoy'
                
                switch experiment
                        
                        case 'all'
                                X0 = array(:, [1:201, 239:242]);
                                Y = array(:, 202:239);
                                
                                for markers = 1:size(Y, 2)
                                        [rho(markers, :), pvalue(markers, :), R2(markers, :), ...
                                                RMSE(markers, :), MAE(markers, :)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end
                                
                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(data.Properties.VariableNames(202:239))', ...
                                        rho, pvalue, R2, RMSE];
                                
                        case 'geneExpOnly'
                                X0 = array(:, 1:201);
                                Y = array(:, 202:239);
                                
                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end
                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(data.Properties.VariableNames(202:239))', ...
                                        rho, pvalue, R2, RMSE];
                                
                        case 'fluxOnly'
                                X0 = array(:, 239:242);
                                Y = array(:, 202:239);
                                
                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end
                                
                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(data.Properties.VariableNames(202:239))', ...
                                        rho, pvalue, R2, RMSE];
                                
                        case 'Me1'
                                markers = string(data.Properties.VariableNames);
                                array = table2cell(data);
                                Me1 = array(:, ~cellfun(@isempty, regexp(markers, 'me1')));
                                Me1Markers = data.Properties.VariableNames(~cellfun(@isempty, ...
                                        regexp(markers, 'me1')));
                                X0 = cell2mat(array(:, [1:201, 239:242]));
                                Y = cell2mat(Me1);
                                
                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end
                                
                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(Me1Markers)', ...
                                        rho, pvalue, R2, RMSE];
                                
                        case 'Me2'
                                markers = string(data.Properties.VariableNames);
                                array = table2cell(data);
                                Me2 = array(:, ~cellfun(@isempty, regexp(markers, 'me2')));
                                Me2Markers = data.Properties.VariableNames(~cellfun(@isempty, ...
                                        regexp(markers, 'me2')));
                                X0 = cell2mat(array(:, [1:201, 239:242]));
                                Y = cell2mat(Me2);
                                
                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end
                                
                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(Me2Markers)', ...
                                        rho, pvalue, R2, RMSE];
                                
                        case 'Me3'
                                markers = string(data.Properties.VariableNames);
                                array = table2cell(data);
                                Me3 = array(:, ~cellfun(@isempty, regexp(markers, 'me3')));
                                Me3Markers = data.Properties.VariableNames(~cellfun(@isempty, ...
                                        regexp(markers, 'me3')));
                                X0 = cell2mat(array(:, [1:201, 239:242]));
                                Y = cell2mat(Me3);
                                
                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end
                                
                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(Me3Markers)', ...
                                        rho, pvalue, R2, RMSE];
                                
                        case 'Ac'
                                markers = string(data.Properties.VariableNames);
                                array = table2cell(data);
                                Ac = array(:, ~cellfun(@isempty, regexp(markers, 'ac')));
                                AcMarkers = data.Properties.VariableNames(~cellfun(@isempty, ...
                                        regexp(markers, 'ac')));
                                X0 = cell2mat(array(:, [1:201, 239:242]));
                                Y = cell2mat(Ac);
                                
                                for markers = 1:size(Y, 2)
                                        [rho(markers,:), pvalue(markers,:), R2(markers,:), ...
                                                RMSE(markers,:), MAE(markers,:)] ...
                                                = crossValidation(X0, Y(:, markers), 10, 0.8);
                                end
                                
                                colNames = ["Markers", "R", "P-value", "R2", "RMSE"];
                                table = [colNames; ...
                                        string(AcMarkers)', ...
                                        rho, pvalue, R2, RMSE];
                end
end

end
