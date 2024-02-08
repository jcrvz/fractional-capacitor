
load experiments.mat

time = dataset(1).t;
voltage = dataset(1).v;
current = dataset(1).i;
rs = 22e-3;

plot(time,voltage,time,current)

%%

s = sym('s');
R  = 7.39e-3;
Ca = 130.21;
Cb = 308.64;
Ra = rand*R;
Rb = rand*R;
Rc = rand*R;
a = 0.2848;
b = 0.866;
D = Ca*Cb*Ra*Rb/Rc;
J = 1;
v0 = 0.36;

H = R + 1/(Ca*s^a) + 1/(Cb*s^b) + 1/(D*s^(a+b));

v = v0 + (R + time.^a/(Ca*gamma(1+a)) + time.^b/(Cb*gamma(1+b)) + ...
    time.^(a+b)./(D*gamma(1+a+b)))*J;

plot(time,v)

%%


load experiments.mat

time = data(1).t;
voltage = data(1).v;
current = data(1).i;

% Parameters
C1  = 58.48;
C2  = 53.44;
Rp  = 100e3;
Rs  = 22e-3;
Rx  = 100;

tau1    = C1*Rx;
tau2    = C2*Rx;
rho     = 1 + Rx/Rp;

A = [-rho/tau1 1/tau1; -1/tau2 1/tau2];
B = [Rx/tau1; 0];
C = [1 0];
D = Rs;

sys = ss(A,B,C,D);

y = lsim(sys,current,time,[0;0]);

plot(time, (voltage - min(voltage))/diff(minmax(voltage(:)')),'k'), hold on,
plot(time, (y - min(y))/diff(minmax(y(:)')),'b')



%%

load experiments.mat

dataset = data(1);

time = dataset.t;
voltage = dataset.v;
current = dataset.i;
currentSource = dataset.is;

% Parameters
C   = 70;
Rp  = 10000e3;
Rs  = 22e-3;
Rx  = 1;

V0  = voltage(1);

Vs  = Rx*currentSource;

rho = 1 + (Rx + Rs)/Rp;
tau = C*(Rx + Rs);

A = -rho/tau;
B = 1/tau; 
C = Rx/(Rx + Rs);
D = Rs/(Rx + Rs);

sys = ss(A,B,C,D); hold on,

y = lsim(sys,Vs,time,V0);

alpha = 0.5;
h = (rho*time/tau).^(alpha-1).*ml(-(rho*time/tau).^alpha,alpha,alpha);

vc_raw = conv(h,Vs,'full');

zsVc = (1/tau)*vc_raw(1:numel(time));

ziVc = V0*h;

totalVc = ziVc + zsVc;


vsc = (Rs*Vs + Rx*totalVc)/(Rx + Rs);

plot(time, (voltage - min(voltage))/diff(minmax(voltage(:)')),'k'), hold on,
plot(time, (y - min(y))/diff(minmax(y(:)')),'b')
plot(time, (vsc - min(vsc))/diff(minmax(vsc(:)')),'r')

%% Respuesta al impulso

load Datosn.mat

dataset = test4;

time = dataset.t;
voltage = dataset.v;
current = dataset.i;

I   = 1;
Rp  = 1000e3;
Rs  = 22e-3;
C   = 58;
tau = Rp*C;

V0  = 0;

v_imp = I*exp(-time/tau);

zsVc = conv(current,v_imp,'same') + current*Rs;

plot(time,voltage,time,zsVc)



%%
% Rs  = 22e-3;
% 
% [t,y] = ode45(@modelSC,[0 150],[0 0]);
% plot(t,Rs*1 + y(:,1))
% 
% function dv = modelSC (t,v)
% 
% % Parameters
% C1  = 58.48;
% C2  = 53.44;
% Rp  = 100e3;
% 
% I   = 1;
% 
% Rx      = 13.10*sqrt(t) + 0.001*(t<=0);
% tau1    = C1*Rx;
% tau2    = C2*Rx;
% rho     = 1 + Rx/Rp;
% 
% dv      = [(-rho*v(1,1) + v(2,1) + Rx*I)./tau1;
%             (-v(2,1) + v(1,1))./tau2];
% 
% end
