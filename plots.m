%% Load data
% Struct: tests
% Fields (Performed experiments): data(1), data(2), ...
% Subfields: t,         Time series [s],
%            i,         Current series [A],
%            v,         Voltage series [V],
%            label,     'Charge' or 'Discharge' behaviour
%            mode,      0 - Charge, 1 - Discharge
%            V0,        Initial internal supercapacitor voltage [V]
%            Vs,        Source voltage [V]
%
load('experiments.mat');

% Set other parameters
legendNames = cell(1,6);
lineColours = lines(6);%{'k','r','b','m','c','g'};
lineStyles  = {'-','-','--','-','--'};

%% Plot experimental data
eSort = [1 2 3 4 5 6];

% Plot each current from tests structure
fiEx1 = Graphics('ExpMeasurements');

% Normalised data
normalisedData = struct();
for dataId = eSort
    
    % Read the corresponding data
    testData = data(dataId);
    
    % Get the voltage measures
    voltage = testData.v;
    current = testData.i;
    time = testData.t;
    
    % Normalise the time measurement
    normalisedData(dataId).maxTime = max(time(:));
    normalisedData(dataId).minTime = min(time(:));
    normalisedData(dataId).time = (time - normalisedData(dataId).minTime)/...
        (normalisedData(dataId).maxTime - normalisedData(dataId).minTime);
    
    % Normalise the voltage measurements
    normalisedData(dataId).maxVoltage = max(voltage(:));
    normalisedData(dataId).minVoltage = min(voltage(:));
    normalisedData(dataId).voltage = (voltage - normalisedData(dataId).minVoltage)/...
        (normalisedData(dataId).maxVoltage - normalisedData(dataId).minVoltage);
    
    % Normalise the current measurements
    normalisedData(dataId).maxCurrent = max(current(:));
    normalisedData(dataId).minCurrent = min(current(:));
    normalisedData(dataId).current = (current - normalisedData(dataId).minCurrent)/...
        (normalisedData(dataId).maxCurrent - normalisedData(dataId).minCurrent);
    
    % Plot the measures
    plot(normalisedData(dataId).time,normalisedData(dataId).voltage,...
        'Color',lineColours(dataId,:),'LineWidth',1,...
        'DisplayName',sprintf('Experiment %d',dataId)); hold on,
    
end
hold off, box on,

% Assign labels and other unimportant things
%xlim([0,3500])
xlabel('Dimensionless time, $$\tilde{t}$$~[s/s]','Interpreter','LaTeX'),
ylabel('Dimensionless voltage, $$\tilde{v}(t)$$~[V/V]','Interpreter','LaTeX'),
leg1 = legend('show'); set(leg1,'Interpreter','LaTeX','box','off')
set(gca,'LineWidth',1,'TickLabelInterpreter','LaTeX');
%axis([0 3500 -0.1 1.1]);
setall(fiEx1,1,[2,1],12,1)

% Plot each current from tests structure
fiEx2 = Graphics('ExpMeasurementsLOG');

for dataId = eSort
    % Select an experiment
    nData = normalisedData(dataId);
    
    % Plot the measures
    semilogy(nData.time,(nData.voltage + 1),'Color',lineColours(dataId,:),...
        'LineWidth',1,'DisplayName',sprintf('Experiment %d',dataId)); hold on,
    
end
hold off, box on,

% Assign labels and other unimportant things
%xlim([0,3500]), ylim([1e-4 1])
xlabel('Dimensionless time, $$\tilde{t}$$~[s/s]','Interpreter','LaTeX'),
ylabel('Dimensionless voltage, $$\log(\tilde{v}(t)+1)$$~[V/V]','Interpreter','LaTeX'),
%legend(legendNames,'Interpreter','LaTeX','box','off')
set(gca,'LineWidth',1,'TickLabelInterpreter','LaTeX');
%axis([0 3500 -0.1 1.1]);
setall(fiEx2,1,[2,1],12,1)

%% Load results
load('processedData.mat');

% Models
Models = {'Caputo-Dzhrbashyan','Caputo-Fabrizio','Atangana-Baleanu-Caputo',...
    'Conformable','Traditional'};
ModelsSh = {'CD','CF','ABC','Conformable','Traditional'};

% Let Tbest from processedData.mat be T
T = Tbest;

% Selection function
selrow     = @(EXPERIMENT, MODEL) T.Experiment == eSort(EXPERIMENT) & ...
    strcmp(T.Model,Models(MODEL));

%% Plot result data

for dataId = eSort
    % Select an experiment
    nData = normalisedData(dataId);
    testData = data(dataId);
    
    % Create the figure
    fiEy = Graphics(sprintf('Models_vsEx%d',dataId));
    
    semilogy(nData.time,(nData.voltage + 1),'--','Color','k','LineWidth',1,...
        'DisplayName','Experimental'); hold on,
    fprintf('\\multirow{5}{*}{%d} ',dataId);
    for modelId  = 1 : numel(Models)
        % Select data
        selection = selrow(dataId,modelId);
        xx = T(selection,3:end);
        
        % Prepare the objective function
        estValues = [xx.C, xx.R_s, xx.R_x, xx.rho, xx.alpha];
        [errValues,vcValues] = voltageModels(estValues,...
            testData,Models{modelId});
        
        % Normalise determined voltage values
        normalisedVoltages = (vcValues - min(vcValues))/...
            (max(vcValues) - min(vcValues));
        
        % Plot the voltage
        semilogy(nData.time,normalisedVoltages + 1,...
            'Color',lineColours(modelId,:),'LineWidth',1,...
            'DisplayName',ModelsSh{modelId},'LineStyle',lineStyles{modelId}); hold on,
        
        fprintf('& %s & %.2f & %.4f & %.4f & %.4f &%.4f & %.4f\\\\ \n',...
            ModelsSh{modelId},xx.Error,xx.C,xx.R_s,xx.R_x,(xx.rho),xx.alpha);
    end
    fprintf('\\midrule\n');
    
    hold off, box on,
    
    % Assign labels and other unimportant things
    %xlim([0,max(nData)])
    xlabel('Dimensionless time, $$\tilde{t}$$','Interpreter','LaTeX'),
    ylabel('Dimensionless voltage, $$\log(\tilde{v}(t)+1)~[V/V]$$','Interpreter','LaTeX'),
    lg = legend(); set(lg,'Interpreter','LaTeX','box','off','location','best');
    set(gca,'LineWidth',1,'TickLabelInterpreter','LaTeX');
    setall(fiEy,1,[2,1],12,1)
    
end

%% Errors

M = reshape(T.Error,5,[])';
M = M(eSort(:),:);

% Plot each current from tests structure
fiEz1 = Graphics('ModelsErrors');

Error = struct('means',mean(M),'stdvs',std(M),'mins',min(M),'max',max(M));

ModelsSh{end-1} = 'Conf.';
ModelsSh{end} = 'Trad.';
boxplot((M),ModelsSh,'Widths',0.5,'PlotStyle','traditional');
% errorbar(1:5,Error.means,Error.stdvs)
% set(gca,'YScale','log');
xlabel('Models','Interpreter','LaTeX'),
ylabel('$$Error=1-R^2$$','Interpreter','LaTeX'),
% ylim([0 0.1])

setall(fiEz1,2,[1,1],12,1)

%%

fiEz2 = Graphics('ExperimentsErrors');
bp2 = bar(M);
for kk = 1 : 5, bp2(kk).DisplayName = ModelsSh{kk}; end

set(gca,'YScale','log');
xlabel('Experiments','Interpreter','LaTeX'),
ylabel('$$Error=1-R^2$$','Interpreter','LaTeX'),
lg = legend(); set(lg,'Interpreter','LaTeX','box','off',...
    'orientation','horizontal','location','north','NumColumns',5);
% ylim([0 105])
setall(fiEz2,1,[2,1],12,1)