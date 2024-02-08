function [errValues,totalSCVoltage] = voltageModels (x, par, modelStr)
% Input:    1. x <- [C, Rs, Rx, rho, alpha]
%           2. par <- struct(V0,Vs,Vk,tk)
%           3. modelStr <- 'Caputo' | 'Caputo-Fabrizio' |
%               'Atangana-Baleanu-Caputo' | 'Conformable' | 'Traditional'*
%           4. N <- 1, 2, 3, ..., 100*, 101, ...
% Output:   1. errValues <- FVU values determined using the chosen model
%           2. vValues <- voltage values determined using the chosen model
%
if nargin < 3
    modelStr = 'Traditional';
end

%% Extract parameters

% Unknown variables
C       = x(1);         % Internal capacitance [F]
Rs      = x(2);         % Internal series resistance [Ohm]
Rx      = x(3);         % Thévenin equivalent resistance [Ohm]
rho     = x(4);         % Voltage divider, rho = 1 + (Rx + Rs)/Rp [1/1]

% Fractional order [1/1]
if numel(x) == 5, alpha = x(5); else, alpha = 1; end     

% Known parameters
V0      = par.Vi;               % Initial internal capacitance voltage [V]
% Vf      = par.Vf;               % Initial internal capacitance voltage [V]
vk      = par.v(:)';            % Experimental voltage data [V]
tk      = par.t(:)';            % Experimental time data [V]
I       = par.Imax;
Vs      = I*Rx;

% Other important parameters

tau     = C*(Rx + Rs);  % Time constant [s] 
td      = rho*tk/tau;   % Dimensionless time series [s/s]

%% Choose the corresponding model (the impulse responses)
% td(1) -> delta(t)

switch modelStr
    case 'Caputo-Dzhrbashyan'
        exponentialPart = ml(-td.^alpha, alpha);
        
    case 'Caputo-Fabrizio'
        exponentialPart = (1/(2-alpha))*exp(-alpha*td/(2-alpha));
        
    case 'Atangana-Baleanu-Caputo'
        exponentialPart = (1/(2-alpha))*ml(-alpha*(td.^alpha)/(2-alpha), alpha);
        
    case 'Conformable'
        exponentialPart = exp(-td.^alpha/alpha);
        
    otherwise % Traditional
        exponentialPart = exp(-td);
end

% Determine the capacitor's internal voltage 
internalCVoltage = V0 + (Vs/rho - V0)*(1 - exponentialPart);

% Determine the total supercapacitor voltage
totalSCVoltage  = I*Rs + internalCVoltage;

% Obtaine the corresponding FVU value
errValues = sum((vk - totalSCVoltage).^2)/sum((vk - mean(vk)).^2);

end
