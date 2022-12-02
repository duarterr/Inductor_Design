clear all;
close all;
warning ('off','all');
clc;

addpath('Datasources');

%% Parâmetros do conversor

Pin = 185.56;                    % W
Vin = 400;                       % V
Pout = 180;                      % W
Vout = 103.68;                   % V
fsw = 100;                       % kHz
D = Vout/Vin;

%% Cálculo de Perdas

% Steinmetz Modificado
% Pv = k * feq^(a-1) * Bmax^b * fsw

% Frequência efetiva para forma de onda de excitação não senoidal
feq = (2/pi^2)*(fsw/(D - D^2));

%% Perdas Volumétricas x Frequência
% Curva para B_max = 0.05T

Bmax1 = 0.05; % T

load B1

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x1 = interp1(x_escalalin',x_escalalog',B1(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y1 = interp1(y_escalalin',y_escalalog',B1(:,2));

% Fit Curve
[xData1, yData1] = prepareCurveData( x1, y1 );

% Set up fittype and options.
ft1 = fittype( 'k*x^a*0.05^b', 'independent', 'x', 'dependent', 'y' );
opts1 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts1.Display = 'Off';
opts1.StartPoint = [0.894702366496388 0.480240315621718 0.957506835434298];

% Fit model to data.
[fitresult1, gof1] = fit( xData1, yData1, ft1, opts1 );

% Plot fit with data.
figure()
h1 = plot( fitresult1, xData1, yData1 );
legend( h1, 'Curve B = 50mT', 'Interpolation', 'Location', 'SouthEast', 'Interpreter', 'none' );
xlabel('Frequency (kHz)');
ylabel('Volumetric Losses (kW/m^3)');
grid on

% Coeficientes de Steinmetz extraídos com Fit Curve Tool
r1 = coeffvalues(fitresult1);
a1 = r1(1);
b1 = r1(2);
k1 = r1(3);

% Perdas no núcleo 
Pcore1 = (k1 * (feq^(a1-1)) * (Bmax1^b1) * fsw)

%% Perdas Volumétricas x Frequência
% Curva para B_max = 0.2T

Bmax2 = 0.2; % T

load B2

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x2 = interp1(x_escalalin',x_escalalog',B2(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y2 = interp1(y_escalalin',y_escalalog',B2(:,2));

% Fit Curve
[xData2, yData2] = prepareCurveData( x2, y2 );

% Set up fittype and options.
ft2 = fittype( 'k*x^a*0.2^b', 'independent', 'x', 'dependent', 'y' );
opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
opts2.StartPoint = [0.0172895593763746 0.211907074257636 0.913375856139019];

% Fit model to data.
[fitresult2, gof2] = fit( xData2, yData2, ft2, opts2 );

% Plot fit with data.
figure()
h2 = plot( fitresult2, xData2, yData2 );
legend( h2, 'Curve B = 200mT', 'Interpolation', 'Location', 'SouthEast', 'Interpreter', 'none' );
xlabel('Frequency (kHz)');
ylabel('Volumetric Losses (kW/m^3)');
grid on

% Coeficientes de Steinmetz extraídos com Fit Curve Tool
r2 = coeffvalues(fitresult2);
a2 = r2(1);
b2 = r2(2);
k2 = r2(3);

% Perdas no núcleo 
Pcore2 = (k2 * (feq^(a2-1)) * (Bmax2^b2) * fsw) 

%% Perdas Volumétricas x Frequência
% Curva para B_max = 0.1T

Bmax3 = 0.1; % T

load B3

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x3 = interp1(x_escalalin',x_escalalog',B3(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y3 = interp1(y_escalalin',y_escalalog',B3(:,2));

% Fit Curve
[xData3, yData3] = prepareCurveData( x3, y3 );

% Set up fittype and options.
ft3 = fittype( 'k*x^a*0.1^b', 'independent', 'x', 'dependent', 'y' );
opts3 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts3.Display = 'Off';
opts3.StartPoint = [0.894702366496388 0.480240315621718 0.957506835434298];

% Fit model to data.
[fitresult3, gof3] = fit( xData3, yData3, ft3, opts3 );

% Plot fit with data.
figure()
h3 = plot( fitresult3, xData3, yData3 );
legend( h3, 'Curve B = 100mT', 'Interpolation', 'Location', 'SouthEast', 'Interpreter', 'none' );
xlabel('Frequency (kHz)');
ylabel('Volumetric Losses (kW/m^3)');
grid on

% Coeficientes de Steinmetz extraídos com Fit Curve Tool
r3 = coeffvalues(fitresult3);
a3 = r3(1);
b3 = r3(2);
k3 = r3(3);

% Perdas no núcleo 
Pcore3 = (k3 * (feq^(a3-1)) * (Bmax3^b3) * fsw) 

%%
% Mostrar gráfico de todas as curvas
figure()
plot(x1,y1);
hold on;
plot(x2,y2);
hold on;
plot(x3,y3);
grid on;
xlabel ('Frequency (kHz)');
ylabel ('Volumetric Losses (kW/m^3)');
title ('Volumetric Losses X Frequency - N97'); 


Flux = linspace(0.05,0.2,130);
Flux = Flux';