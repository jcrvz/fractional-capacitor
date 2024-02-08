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

%% Boudaries selection
%{
    From maxwell.com/images/documents/datasheet_16v_small_cell_module.pdf
    Capacitance is about 58-70 F (350 F/cell with 6 cells in parallel) and
    ESR (Rs) is about 22 mOhm. But, after 10 years of use at ideal conditons
    the capacitance is reduced 20 % from its minimum initial value and the
    ESR is increased 100% from its maximum initial value. Then,

        0.8 * 58 F <= C <= 70 F, and 22 mOhm <= Rs <= 2*22 mOhm.

    However, considering a practical case a bit worse than the ideal, we
    could say that

        0.5 * 58 F <= C <= 70 F, and 22 mOhm <= Rs <= 10*22 mOhm.

    Maximum Rs value changed a lot because it can include wire resistances.

    Otherwise, rho = 1 + Rs/Rp has a well-defined lower limit equals 1,
    which corresponds to the ideal case when Rp tends infinity. The rho
    upper limit can be found for the closest worst case when Rp = Rs, so
    rho = 2.

    Rx is the real source's Thévenin equivalent resistance, which is
    assumed in a range of

                        1 Ohm <= Rx <= 1000 Ohm.

    Finally, alpha must be between 0 and 1. Therefore,
    
                   1 <= rho <= 2 and 0 < alpha <= 1.
%}

%    Boundaries     = [C,       Rs,         Rx,     rho,    alpha]
lowerBoundaries     = [0.001*58,  22e-3,      1e-3,   1,      0.5];
upperBoundaries     = [1.00*70,  2*22e-3,    100e6,    2,      1.0];

boundaries          = [lowerBoundaries',upperBoundaries'];

%% Create data files
datafilename = sprintf('results_%s.txt',datestr(now,'ddmmyy'));
fid = fopen(datafilename,'a+');
fprintf(fid,'%s\n-----\n',datestr(now));              % Print date
fprintf(fid,'%s\n',sprintf('%20s\t', 'Rep', 'Experiment', 'Model', ...
    'Time', 'Error_Value', 'X_values'));

%% Initialise some features
Models = {'Caputo-Dzhrbashyan','Caputo-Fabrizio','Atangana-Baleanu-Caputo', ...
    'Conformable','Traditional'};

%% Run tests

% Define some additional parameters
options = optimoptions('fmincon','Display','none','Algorithm',...
    'interior-point','OptimalityTolerance',1e-10,'UseParallel',false);

for dataId = 1 : 6
    % Select an experiment
    testData    = data(dataId);
    
    for modelId  = 1 : numel(Models)
        
        % Prepare the objective function
        objectiveFunction = @(x) voltageModels(x,testData,Models{modelId});
        
        % Repeat the fitting procedure
        for repetitions = 1 : 1
            
            % If it is the traditional model, alpha parameter is disregarded.
%             xStart = [58 22e-3 100 1.1 0.9];
            if strcmp(Models{modelId},'Traditional'), bnd = boundaries(1:end-1,:);
            else, bnd = boundaries; end
            
            % Optimisation procedure 1: Find the initial guess
            tic,
%             [xSol1,fSol1] = ga(objectiveFunction,size(bnd,1),[],[],[],[],...
%                 bnd(:,1)',bnd(:,2)',[],[]);
%             xSol1 = inf; fSol1 = inf;
            [xSol1,fSol1,details] = CSOA(objectiveFunction,bnd);%,struct('visualMode',true));
            time1 = details.elapsedTime;
%             time1 = toc;
            
            % Print found initial guess solution
            x1ValuesString = sprintf('%10.6f\t',xSol1);
            fprintf('St%3d (Ex-%d, %s) :: x_opt = [%s], f_opt = %.4g, time = %.2f s\n', ...
                repetitions,dataId,Models{modelId},x1ValuesString,fSol1,time1);
            
            % Try to refine this solution
            tic,
            [xSol2,fSol2] = fmincon(objectiveFunction,xSol1,[],[],[],[],...
                bnd(:,1)',bnd(:,2)',[],options);
            time2 = toc;
            
            % Print found refined solution
            x2ValuesString  = sprintf('%10.6f\t',xSol2);
            fprintf('Dt%3d (Ex-%d, %s) :: x_opt = [%s], f_opt = %.4g, time = %.2f s\n', ...
                repetitions,dataId,Models{modelId},x2ValuesString,fSol2,time2);
            
            % Check which one is the best found solution
            if fSol2 < fSol1
                fOpt = fSol2;
                xOptString = x2ValuesString;
            else
                fOpt = fSol1;
                xOptString = x1ValuesString;
            end
            
            % Save results
            fprintf(fid,'%3d\t%10d\t%30s\t%10.2f\t%10.4g\t%50s\n',...
                repetitions,dataId,Models{modelId},time1+time2,fOpt,xOptString);
        end
    end
end

fclose(fid);


