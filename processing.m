%% Load results
name_file_results = 'results_220619.txt'; 

T = readtable(name_file_results,'Format','%u%u%s%f%f%f%f%f%f%f');
varNames = {'Repetition','Experiment','Model','Time','Error','C','R_s','R_x','rho','alpha'};
T.Properties.VariableNames = varNames;
T.Properties.VariableDescriptions = ...
    {'Rep.','Experiment','Model','$$t$$','$$Error$$','$$C$$','$$R_s$$','$$R_x$$','$$\rho$$','$$\alpha$$'};
T.Properties.VariableUnits = ...
    {'','','','s','','F','$$\Omega$$','$$\Omega$$','',''};

% Models
Models = {'Caputo-Dzhrbashyan','Caputo-Fabrizio','Atangana-Baleanu-Caputo','Conformable','Traditional'};
ModelsSh = {'CD','CF','ABC','Conformable','Traditional'};

% Selection function
selrow     = @(EXPERIMENT, MODEL) T.Experiment == EXPERIMENT & ...
    strcmp(T.Model,Models{MODEL});

% Replace NaNs with Ones
T.alpha(isnan(T.alpha)) = 1;

% Get experiments and models
Experiments = unique(T.Experiment);
% Find stats
mods4T = {}; exps4T = [];
means4T = []; stds4T = []; mins4T = []; maxs4T = []; medians4T = []; bests4T = [];
for iiE = 1 : numel(Experiments)
    for iiM = 1 : numel(Models)
        miniT = T(selrow(Experiments(iiE),iiM),4:end);
        
        chunkData = table2array(miniT);
        
        mods4T(end+1,:)     = Models(iiM);
        exps4T(end+1,:)     = Experiments(iiE);
        means4T(end+1,:)    = mean(chunkData,1);
        stds4T(end+1,:)     = std(chunkData,[],1);
        mins4T(end+1,:)     = min(chunkData,[],1);
        maxs4T(end+1,:)     = max(chunkData,[],1);
        medians4T(end+1,:)  = median(chunkData,1);
        
        [~,idmin] = min(miniT.Error);
        bests4T(end+1,:) = table2array(miniT(idmin,:));
        %disp([sprintf('%2d\t%30s\t',Experiments(iiE),Models{iiM}),...
%             num2str(bests4T(end,:))]);
    end
end

%% The best table

Tbest = table(exps4T,mods4T, ...
    bests4T(:,1), ...
    bests4T(:,2), ...
    bests4T(:,3), ...
    bests4T(:,4), ...
    bests4T(:,5), ...
    bests4T(:,6), ...
    bests4T(:,7));
Tbest.Properties.VariableNames = varNames(2:end);

disp(Tbest)

%% The mean and stdv tables
%means4T(:,1), stds4T(:,1), mins4T(:,1), maxs4T(:,1), ...
Tadd = table(exps4T,mods4T, ...    
    means4T(:,2), stds4T(:,2), mins4T(:,2), maxs4T(:,2), ...
    means4T(:,3), stds4T(:,3), mins4T(:,3), maxs4T(:,3), ...
    means4T(:,4), stds4T(:,4), mins4T(:,1), maxs4T(:,4),...
    means4T(:,5), stds4T(:,5), mins4T(:,5), maxs4T(:,5), ...
    means4T(:,6), stds4T(:,6), mins4T(:,6), maxs4T(:,6));
    %'Time_avg', 'Time_std', 'Time_min', 'Time_max', ...
Tadd.Properties.VariableNames = {'Experiment','Model', ... 
    'Error_avg',  'Error_std',  'Error_min',  'Error_max', ...
    'C_avg',    'C_std',    'C_min',    'C_max', ...
    'R_s_avg',  'R_s_std',  'R_s_min',  'R_s_max', ...
    'rho_avg',  'rho_std',  'rho_min',   'rho_max', ...
    'alpha_avg','alpha_std','alpha_min', 'alpha_max'};

%% Save results
save('processedData.mat','Tbest','Tadd');

%%
% figure, 
% for iiE = 1 : numel(Experiments)
%         miniT = Tadd(Tadd.Experiment == Experiments(iiE), 2:end);
%         errorbar(miniT.C_avg,miniT.C_std), hold on,
% end
