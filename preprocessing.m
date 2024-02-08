%% Load raw data
load Datosn

% Define data modes: % 0 - Charge, 1 - Discharge
modes       = [0 0 0 1 1 1];

% Starting time and ending time
t0 = [1.09 3.3 1.25 2.59 2.14 3.69];
t1 = [81.52 103.10 158.20 55.37 162.20 32.81];

% Initial voltages
ViValues = [0 0 0 16 16 16];

% Final voltages
VfValues = [16 16 16 0 0 0];

% Current values
currentValues = [3 2 1 5 1 10];

% Store charge data
tests(1)    = test4;
tests(2)    = test5;
tests(3)    = test6;

% Store discharge data
tests(4)    = test1;        
tests(5)    = test2;        
tests(6)    = test3;

%% Trim raw data

% Store data
decimationFactor = 100;
data = struct();
for id = 1 : 6
    test = tests(id);
    
    % Measured data
    time    = test.t;    
    voltage = test.v;
    current = test.i;
    
    if modes(id), label = 'Discharge'; else, label = 'Charge'; end
    
    startingTimeIndex  = find(time == t0(id));
    endingTimeIndex    = find(time == t1(id)); % 
%     endingTimeIndex    = numel(time); % default
    
    Vi  = ViValues(id); %mean(voltage(1 : startingTimeIndex));
    Vf  = VfValues(id); %mean(voltage(time > t1(id)));
    
    % Ideal current source
    currentSource   = ((-1)^modes(id))*currentValues(id)*double(time >= t0(id) & time <= t1(id));
    
    % Store data
    data(id).label  = label;
%     data(id).mode   = modes(id);
    data(id).Imax   = ((-1)^modes(id))*currentValues(id);
    data(id).Vi     = Vi;
    data(id).Vf     = Vf;
    data(id).t      = time(startingTimeIndex : endingTimeIndex) - t0(id);    % Trimmed time
    data(id).v      = voltage(startingTimeIndex : endingTimeIndex); % Trimmed voltage
    data(id).i      = ((-1)^modes(id))*(current(startingTimeIndex : endingTimeIndex)); % Trimmed current
    
    % Ideal current source
    data(id).is     = currentSource(startingTimeIndex : endingTimeIndex);
    
    % Decimate data
    data(id).t      = data(id).t(1 : decimationFactor : end);
    data(id).i      = data(id).i(1 : decimationFactor : end);
    data(id).v      = data(id).v(1 : decimationFactor : end);
    data(id).is     = data(id).is(1 : decimationFactor : end);
    
    % Just to visualise
    figure,
    subplot(2,1,1)
    plot(time,current,'b','LineWidth',1), hold on,
    plot(data(id).t,data(id).i,'c--','LineWidth',.6),
    plot(data(id).t,data(id).is,'r--','LineWidth',1),
    title(sprintf('Dataset %d: %s current',id,label));
    xlabel('[s]'), ylabel('[A]'), legend('Raw data','Trimmed data');
    
    subplot(2,1,2)
    plot(time,voltage,'b','LineWidth',1), hold on,
    plot(data(id).t,data(id).v,'c--','LineWidth',.6),
    line(t0(id)*[1 1],sort([Vi Vf]),'LineStyle','--','Color','k','LineWidth',1)
    line(t1(id)*[1 1],sort([Vi Vf]),'LineStyle','--','Color','g','LineWidth',1)
    line([time(1) time(end)],Vi*[1 1],'LineStyle','--','Color','r','LineWidth',1)
    line([time(1) time(end)],Vf*[1 1],'LineStyle','--','Color','m','LineWidth',1)
    title(sprintf('Dataset %d: %s voltage',id,label));
    xlabel('[s]'), ylabel('[V]'), legend('Raw data','Trimmed data');
end

save('experiments.mat','data');
