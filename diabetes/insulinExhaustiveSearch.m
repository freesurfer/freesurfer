Subject = 4;
Visits = 2;
gpuIndex = 1; % gpuIndex = 1 (or index of the GPU) will use the GPU. gpuIndex = 0 will use the CPU.
showPlot = true; % Disable to gain speed
nDim = 6; % Use matrix operations in the search in the first 'nDim' dimensions.

load SubjectData
param = param(Subject,Visits);

varNames = {'k1',  'k3',  'k45', 'dG0',          'SC',  'SI', 'b1',     'b2',     'bTime',            'cTime'};
varDelta = [0.001, 0.003, 0.001, 1,              .5,    1,    .025 /60, .025 /60, 30,                 -1];
varMin =   [0.005, 0.001, 0.005, 0*varDelta(4),  1.5,   1,    .3 /60,   .3 /60,   0,                  15];
varMax =   [0.05,  0.01,  0.04,  0*varDelta(4),  15,    100,  1.5 /60,  1.5 /60, param(1).time(end),  5];
nVisits = size(param,2);

if gpuIndex
    gpu = gpuDevice(gpuIndex);
    disp(['Using the ' gpu.Name ' GPU...'])
    for v = 1:nVisits
        paramFieldNames = fieldnames(param(v));
        for n = 1:length(paramFieldNames)
            param(v).(paramFieldNames{n}) = gpuArray(param(v).(paramFieldNames{n}));
        end
    end
end
nVars = length(varNames);
varRange = cell(1,nVars);
for n = 1:nVars
    varRange{n} = varMin(n):varDelta(n):varMax(n);
    varRange{n} = permute(varRange{n}', circshift(1:nVars, n-1));
    if gpuIndex
        varRange{n} = gpuArray(varRange{n});
    end
end
L = cellfun(@length, varRange);
for v = 1:nVisits
    param(v).time = permute(param(v).time', circshift(1:(nVars+1), nVars));
    param(v).BG   = permute(param(v).BG',   circshift(1:(nVars+1), nVars));
end
sub = cell(1, nVars - nDim);
subMin = cell(1, nDim);
totalMin = inf;
nInd = prod(L(nDim+1:end));
xl = gather([param(1).time(1) param(1).time(end)]);
tic
for ind = 1:nInd
    [sub{:}] = ind2sub(L(nDim+1:end), ind);
    varSingle = cellfun(@(c1,c2) c1(c2), varRange(nDim+1:end), sub, 'UniformOutput', false);
    E = insulinFitError(varRange{1:nDim}, varSingle{:}, param(1));
    for v = 2:nVisits
        E = E + insulinFitError(varRange{1:nDim}, varSingle{:}, param(v));
    end
    [m, indMin] = min(E(:));
    if m < totalMin
        totalMin = m;
        [subMin{:}] = ind2sub(L(1:nDim), indMin);
        totalMinVar = [cellfun(@(c1,c2) c1(c2), varRange(1:nDim), subMin, 'UniformOutput', false) varSingle];
        strCell = cellfun(@(c1, c2) [c1 ': ' num2str(c2) ', '], varNames, totalMinVar, 'UniformOutput', false);
        disp(['New minimum after ' num2str(toc/60) ' minutes, at index ' num2str(ind) ' out of ' num2str(nInd) ':'])
        str = [strCell{:} 'Error: ' num2str(totalMin)];
        disp(str)
        if showPlot
            for v = 1:nVisits
                [~, BG] = insulinFitError(totalMinVar{:}, param(v));
                subplot(1,nVisits,v)
                plot(param(v).time(:), BG(:), param(v).time(:), param(v).BG(:))
                legend Estimated Measured
                xlabel('Time (min)')
                ylabel('Blood Glucose (mg/dL)')
                title(['Visit ' num2str(v)])
                xlim(xl)
                yl(v,:) = ylim;
            end
            ylAll = [min(yl(:,1)) max(yl(:,2))];
            for v = 1:nVisits
                subplot(1,nVisits,v)
                ylim(ylAll)
            end
            drawnow
            set(gcf, 'Name', str)
        end
    end
end
if gpuIndex
    totalMinVar = cellfun(@gather, totalMinVar, 'UniformOutput', false);
end
save(['insulinFit_Subject' num2str(Subject) '_Visits' strrep(num2str(Visits), '  ', '-')], 'totalMinVar', 'varNames')

disp(['Insulin sensitivity: ' num2str(totalMinVar{find(strcmp(varNames, 'SI'))}) ' mg/dL / u'])
disp(['Insulin-to-carb ratio: ' num2str(totalMinVar{find(strcmp(varNames, 'SI'))} / totalMinVar{find(strcmp(varNames, 'SC'))}) ' g/u'])
disp(['(Carb sensitivity: ' num2str(totalMinVar{find(strcmp(varNames, 'SC'))}) ' mg/dL / g)'])
disp(['Basal: ' num2str(totalMinVar{find(strcmp(varNames, 'b1'))}*60) ' u/hr before <' num2str(totalMinVar{find(strcmp(varNames, 'bTime'))}) ' min>, ' num2str(totalMinVar{find(strcmp(varNames, 'b2'))}*60) ' u/hr after that.'])
disp(['Carb-ingestion delay: ' num2str(totalMinVar{find(strcmp(varNames, 'cTime'))}) ' min'])
