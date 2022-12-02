clear all;
close all;
warning ('off','all');
clc;

addpath('Datasources');

%% Perdas Volumétricas x Frequência
% Curva para B_max = 0.05T

load B1

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x1 = interp1(x_escalalin',x_escalalog',B1(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y1 = interp1(y_escalalin',y_escalalog',B1(:,2));

% Mostrar gráfico
plot(x1,y1);
hold on;
grid on;
xlabel ('Frequency (kHz)');
ylabel ('Volumetric Losses (kW/m^3)');
title ('Volumetric Losses X Frequency - N97'); %azul


%% Perdas Volumétricas x Frequência
% Curva para B_max = 0.2T

load B2

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x2 = interp1(x_escalalin',x_escalalog',B2(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y2 = interp1(y_escalalin',y_escalalog',B2(:,2));

% Mostrar gráfico
plot(x2,y2);
hold on;
grid on;
xlabel ('Frequency (kHz)');
ylabel ('Volumetric Losses (kW/m^3)');
title ('Volumetric Losses X Frequency - N97'); %vermelha


%% Perdas Volumétricas x Frequência
% Curva para B_max = 0.1T

load B3

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x3 = interp1(x_escalalin',x_escalalog',B3(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y3 = interp1(y_escalalin',y_escalalog',B3(:,2));

% Mostrar gráfico

plot(x3,y3);
hold on;
grid on;
xlabel ('Frequency (kHz)');
ylabel ('Volumetric Losses (kW/m^3)');
title ('Volumetric Losses X Frequency - N97'); %amarela


%% Perdas Volumétricas para cada curva

Target.f_min = 50;	% Min frequency - kHz
Target.f_max = 300;	% Max frequency - kHz
Range_f = Target.f_min:50:Target.f_max;

for i = Target.f_min:Target.f_max
    Pv_B1 = interp1(x1,y1,Range_f);
    Pv_B2 = interp1(x2,y2,Range_f);
    Pv_B3 = interp1(x3,y3,Range_f);
    i = i+1;
end

%% Pontos para cálculos dos coeficientes de Steinmetz

f1 = 200; %kHz
f3 = 100; %kHz
f4 = 300; %kHz

B_1 = 0.5e-3; %T
B_2 = 200e-3; %T
B_3 = 100e-3; %T

P_p1 = interp1(x1,y1,f1);
P_p2 = interp1(x2,y2,f1);
P_p3 = interp1(x3,y3,f3);
P_p4 = interp1(x3,y3,f4);


%% Cálculo dos coeficientes de Steinmetz

% Frequência fixa f1 = 200kHz, varia fluxo
Beta = (log(P_p1)-log(P_p2))/(log(B_1)-log(B_2))

% Fluxo magnético fixo B_3 = 100mT, varia frequência
Alpha = (log(P_p3)-log(P_p4))/(log(f1)-log(f4))

% Coeficiente Cm
Cm = 10^((-Beta*log(B_3))+(-Alpha*log(f3))+log(P_p3))


%% Cálculo das perdas volumétricas para volume, frequência e fluxo fixos

fsw = 100e3;
Bm = 200e-3;
V = 1;
Tsw = 1/fsw;

% tsw = 0:1e-3:Tsw
% 
% Range_B = 0:1e-3:200e-3;
% 
% g = (diff(Range_B)).^2;
% 
% c = integral(g,0,Tsw)
% 
% feq = (2/(((200e-3)^2)*(pi^2)))   
% 
% Pv_B2
% 
% % Simplified Steinmetz Equation - Tese Renan
% Volumetric_Loss_MSE = Cm * (feq^(Alpha)) * (Bm^Beta) * fsw 
% 
% Steinmetz Equation - Tese Maikel
Volumetric_Loss_SE =  (pi/4) * Cm * (fsw^Alpha) * (Bm^Beta) 




























