function [i_funct,err_funct] = models (parameters,modelStr,Napp)
% Input:    1. parameters <- struct(E0,I,t,C,Re,Rs)
%           2. modelStr <- 'Caputo' | 'Caputo-Fabrizio' |
%               'Atangana-Baleanu-Caputo' | 'Conformable' | 'Traditional'*
%           3. Napp <- 1, 2, 3, ..., 100*
% Output:   1. i_funct <- handle function @(t,x), x = (A,alpha)
%           2. err_funct <- handle function @(x) // R-squared
if nargin < 3
    Napp = 100;
    if nargin < 2
        modelStr = 'Traditional';
    end
end

%% Define some important parameters
% Time constant [in seconds]
tau     = @(x) x(1)*(parameters.Re + x(2));

% Magnitude constant [dimensionles]
%A_const = @(x) 1 + (parameters.Re + x(2))/(10^x(3));

% Total Sum of Squares (SStot)
SStot   = sum((parameters.I - mean(parameters.I)).^2);

% Define the Current's functions
% x = [C, Rs, rho, alpha]
switch modelStr
    case 'Caputo-Dzhrbashyan'      
        i_funct = @(t,x) (parameters.E0*x(1)/(tau(x)*x(3))) * ...
            Exp( -(x(3)*t(:)'/tau(x)).^(x(4)) , x(4), Napp);
%         wrong (?) model
%         i_funct = @(t,x) (parameters.E0*x(1)/(tau(x)*x(3))) * ...
%             Exp( -x(3)*(t(:)'/tau(x)).^(x(4)) , x(4), Napp);
        
    case 'Caputo-Fabrizio'
        i_funct = @(t,x) (parameters.E0*x(1)/(tau(x)*x(3))) * ...
            (1/(2 - x(4))) * exp( -(x(3)*x(4)/(2 - x(4)))*t(:)'/tau(x) );
%         wrong (?) model
%         i_funct = @(t,x) (parameters.E0*x(1)/( ...
%             tau(x)*x(3) + (tau(x)^(1 - x(4)))*(x(3)^2)*(1 - x(4)))) * ...
%             exp( -x(3)*x(4)*t(:)'/(tau(x)^x(4) + x(3)*(1 - x(4))));
        
    case 'Atangana-Baleanu-Caputo'
        i_funct = @(t,x) (parameters.E0*x(1)/(tau(x)*x(3)))* ...
            (1/(2 - x(4)))*Exp( - (x(4)/(2 - x(4)))*(x(3)*t(:)'/tau(x)).^(x(4)),x(4), Napp);
%         wrong (?) model
%         i_funct = @(t,x) (parameters.E0*x(1)/( ...
%             tau(x)*x(3) + (tau(x)^(1 - x(4)))*(x(3)^2)*(1 - x(4)))) * ...
%             Exp( -x(3)*x(4)*(t(:)').^x(4)/(tau(x)^x(4) + x(3)*(1 - x(4))), x(4), Napp);
        
    case 'Conformable'
        i_funct = @(t,x) (parameters.E0*x(1)/(tau(x)*x(3)))* ...
            exp( -(1/x(4))*(x(3)*t(:)'/tau(x)).^(x(4)) );
%         wrong (?) model
%         i_funct = @(t,x) (parameters.E0*x(1)/(tau(x)*x(3))) * ...
%             Exp( -(x(3)/x(4))*(t(:)'/tau(x)).^(x(4)) , 1, Napp);
        
    otherwise % Traditional
        i_funct = @(t,x) (parameters.E0*x(1)/(tau(x)*x(3))) * ...
            exp( -x(3)*t(:)'/tau(x) );
%         wrong (?) model
%         i_funct = @(t,x) (parameters.E0*x(1)/(tau(x)*x(3))) * ...
%             exp( -(x(3)*t(:)'/tau(x)));
end

% R-squared's function
err_funct = @(x) sum((parameters.I - i_funct(parameters.t,x)).^2)/SStot;

end

%% Simple Mittag Leffler Function
function E = Exp(z,alpha,N)
k   = repmat((0:N-1)',1,length(z));
E   = sum(repmat(z,N,1).^k./gamma(alpha*k + 1));
end

% Alpha   = x(1);
% Sigma   = x(2);
%
% ik      = nan(size(tk));
% n       = 0 : NN;
%
% for k = 1 : length(tk)
%
%     cumSeries = (((-1).^n)./gamma(n*Alpha + 1)).*...
%         ((Sigma/(Re*C)).^((1 - Alpha)*n)).*((tk(k)/(Re*C)).^(Alpha*n));
%
%     ik(k)   = sum(cumSeries);
% end
%
% ik  = E0*ik/Re;
