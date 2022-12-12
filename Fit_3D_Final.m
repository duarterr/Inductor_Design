clear all;
clc;
close all;
addpath('Datasources');

%% Importação dos dados das curvas do datasheet

load Data_25mT.mat
load Data_50mT.mat
load Data_100mT.mat
load Data_200mT.mat
load Data_300mT.mat

%% Ajuste dos pontos 

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);

% Interpolação dos pontos - *1000 para obter em Hz e W/m3
x25 = interp1(x_escalalin',x_escalalog',Data_25mT(:,1))*1000;
y25 = interp1(y_escalalin',y_escalalog',Data_25mT(:,2))*1000;

x50 = interp1(x_escalalin',x_escalalog',Data_50mT(:,1))*1000;
y50 = interp1(y_escalalin',y_escalalog',Data_50mT(:,2))*1000;

x100 = interp1(x_escalalin',x_escalalog',Data_100mT(:,1))*1000;
y100 = interp1(y_escalalin',y_escalalog',Data_100mT(:,2))*1000;

x200 = interp1(x_escalalin',x_escalalog',Data_200mT(:,1))*1000;
y200 = interp1(y_escalalin',y_escalalog',Data_200mT(:,2))*1000;

x300 = interp1(x_escalalin',x_escalalog',Data_300mT(:,1))*1000;
y300 = interp1(y_escalalin',y_escalalog',Data_300mT(:,2))*1000;

clear x_escalalog x_escalalin y_escalalog y_escalalin;

%%

figure;
hold on;
grid on;
plot (x25,y25);
plot (x50,y50);
plot (x100,y100);
plot (x200,y200);
plot (x300,y300);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
%xlim([10^1 10^3]);
%ylim([10^-1 10^4]); 
legend ({'25mT', '50mT', '100mT', '200mT', '300mT'});

%% Fit: '25mT'.

% Set up fittype and options.
ft = fittype( 'k*f^a*0.025^b', 'independent', 'f', 'dependent', 'P' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
%opts.StartPoint = [0.14399839834894 0.73773725496228 0.694828622975817];

% Fit model to data.
[fitresult, ~] = fit( x25, y25, ft, opts );

% Plot fit with data.
figure( 'Name', '25mT' );
h = plot( fitresult, x25, y25 );
legend( h, 'P vs. f', '25mT', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'f', 'Interpreter', 'none' );
ylabel( 'P', 'Interpreter', 'none' );
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on


%% Parâmetros do conversor

Vin = 400;                       % V
Vout = 103.68;                   % V
D = Vout/Vin;
fsw = 100;                       % kHz
%Ve = 7.46e-6;

% Cálculo de Perdas - Steinmetz Modificado
% Pv = k * feq^(a-1) * Bmax^b * fsw
% Frequência efetiva para forma de onda de excitação não senoidal
feq = (2/pi^2)*(fsw/(D-D^2));

%% Perdas Volumétricas x Frequência

% CURVA PARA B_MAX = 0.025T
Bmax1 = 0.025; % T
load Data_25mT

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x25 = interp1(x_escalalin',x_escalalog',Data_25mT(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y25 = interp1(y_escalalin',y_escalalog',Data_25mT(:,2));

% Fit Curve
[xData1, yData1] = prepareCurveData( x25, y25 );

% Set up fittype and options.
ft1 = fittype( 'k*x^a*0.025^b', 'independent', 'x', 'dependent', 'y' );
opts1 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts1.Display = 'Off';
opts1.StartPoint = [0.14399839834894 0.73773725496228 0.694828622975817];

% Fit model to data.
[fitresult1, gof1] = fit( xData1, yData1, ft1, opts1 );

% Coeficientes de Steinmetz extraídos com Fit Curve Tool
r1 = coeffvalues(fitresult1);
a1 = r1(1);
b1 = r1(2);
k1 = r1(3);


% CURVA PRA B_MAX = 0.05T
Bmax2 = 0.05; % T
load Data_50mT

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x2 = interp1(x_escalalin',x_escalalog',Data_50mT(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y2 = interp1(y_escalalin',y_escalalog',Data_50mT(:,2));

% Fit Curve
[xData2, yData2] = prepareCurveData( x2, y2 );

% Set up fittype and options.
ft2 = fittype( 'k*x^a*0.05^b', 'independent', 'x', 'dependent', 'y' );
opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
opts2.StartPoint = [0.0843065411551075 0.739014598693051 0.438744359656398];

% Fit model to data.
[fitresult2, gof2] = fit( xData2, yData2, ft2, opts2 );

% Coeficientes de Steinmetz extraídos com Fit Curve Tool
r2 = coeffvalues(fitresult2);
a2 = r2(1);
b2 = r2(2);
k2 = r2(3);


% CURVA PARA B_MAX = 0.1T
Bmax3 = 0.1; % T
load Data_100mT

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x3 = interp1(x_escalalin',x_escalalog',Data_100mT(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y3 = interp1(y_escalalin',y_escalalog',Data_100mT(:,2));

% Fit Curve
[xData3, yData3] = prepareCurveData( x3, y3 );

% Set up fittype and options.
ft3 = fittype( 'k*x^a*0.1^b', 'independent', 'x', 'dependent', 'y' );
opts3 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts3.Display = 'Off';
opts3.StartPoint = [0.112145636748411 0.47831407362237 0.186872604554379];

% Fit model to data.
[fitresult3, gof3] = fit( xData3, yData3, ft3, opts3 );

% Coeficientes de Steinmetz extraídos com Fit Curve Tool
r3 = coeffvalues(fitresult3);
a3 = r3(1);
b3 = r3(2);
k3 = r3(3);


% CURVA PARA B_MAX = 0.025T
Bmax4 = 0.2; % T
load Data_200mT

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x4 = interp1(x_escalalin',x_escalalog',Data_200mT(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y4 = interp1(y_escalalin',y_escalalog',Data_200mT(:,2));

% Fit Curve
[xData4, yData4] = prepareCurveData( x4, y4 );

% Set up fittype and options.
ft4 = fittype( 'k*x^a*0.2^b', 'independent', 'x', 'dependent', 'y' );
opts4 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts4.Display = 'Off';
opts4.StartPoint = [0.747489478402464 0.126820712754585 0.709364830858073];

% Fit model to data.
[fitresult4, gof4] = fit( xData4, yData4, ft4, opts4 );

% Coeficientes de Steinmetz extraídos com Fit Curve Tool
r4 = coeffvalues(fitresult4);
a4 = r4(1);
b4 = r4(2);
k4 = r4(3);


% CURVA PARA B_MAX = 0.3T
Bmax5 = 0.3; % T
load Data_300mT

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x5 = interp1(x_escalalin',x_escalalog',Data_300mT(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y5 = interp1(y_escalalin',y_escalalog',Data_300mT(:,2));

% Fit Curve
[xData5, yData5] = prepareCurveData( x5, y5 );

% Set up fittype and options.
ft5 = fittype( 'k*x^a*0.3^b', 'independent', 'x', 'dependent', 'y' );
opts5 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts5.Display = 'Off';
opts5.StartPoint = [0.89883322382211 0.294754399637592 0.505957051665142];

% Fit model to data.
[fitresult5, gof5] = fit( xData5, yData5, ft5, opts5 );

% Coeficientes de Steinmetz extraídos com Fit Curve Tool
r5 = coeffvalues(fitresult5);
a5 = r5(1);
b5 = r5(2);
k5 = r5(3);


% Plot 5 Interpolações
figure()
h1 = plot( fitresult1, 'c', xData1, yData1 );
hold on;
h2 = plot( fitresult2, 'y', xData2, yData2 );
hold on;
h3 = plot( fitresult3, 'r', xData3, yData3 );
hold on;
h4 = plot( fitresult4, 'b', xData4, yData4 );
hold on;
h5 = plot( fitresult5, 'g', xData5, yData5 );
grid on;
xlabel ('Frequency (kHz)');
ylabel ('Volumetric Losses (kW/m^3)');
title ('Volumetric Losses X Frequency - N97');
legend(['Data - 25mT'],['Curve - 25mT'],['Data - 50mT'],['Curve - 50mT'],['Data - 100mT'],['Curve - 100mT'],['Data - 200mT'],['Curve - 200mT'],['Data - 300mT'],['Curve - 300mT'])
%legend('Curve B = 50mT', 'Interpolation', 'Curve B = 200mT', 'Interpolation', 'Curve B = 100mT', 'Interpolation' );

%% Fit Curve 
% Superfície com dimensões x,y,z = frequência, fluxo magnético, perdas
% volumétricas, a partir dos coeficientes de Steinmetz de cada curva
% extraídos anteriormente

% Coeficientes K 
k_25 = fitresult1.k;
k_50 = fitresult2.k;
k_100 = fitresult3.k;
k_200 = fitresult4.k;
k_300 = fitresult5.k;

% Coeficientes Alpha
a_25 = fitresult1.a;
a_50 = fitresult2.a;
a_100 = fitresult3.a;
a_200 = fitresult4.a;
a_300 = fitresult5.a;

% Coeficientes Beta
b_25 = fitresult1.b;
b_50 = fitresult2.b;
b_100 = fitresult3.b;
b_200 = fitresult4.b;
b_300 = fitresult5.b;

% Vetor de frequências
f_range = 20:1:300;
% Vetor dos fluxo magnético de cada curva
B_range = [0.025 0.05 0.1 0.2 0.3];

% Matriz de cada coeficiente
k = [k_25 k_50 k_100 k_200 k_300];
a = [a_25 a_50 a_100 a_200 a_300];
b = [b_25 b_50 b_100 b_200 b_300];

% Cálculo das perdas por Steinmetz Clássico
for B = 1:numel(B_range)
   for f = 1:numel(f_range)
        P(B,f) = k(B)*f_range(f)^a(B)*B_range(B)^b(B);  
   end
end
%P = k*x^a*y^b
% x = frequência
% y = fluxo magnético

%% Fit: 'Superfície'.
[xData, yData, zData] = prepareSurfaceData( f_range, B_range, P );

% Set up fittype and options.
ft = fittype( 'k*x^a*y^b', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.260428644600092 0.940333711980278 0.547215529963803];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Surface' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 'Surface', 'Volumetric Losses (kW/m^3) X Frequency (kHz) X Magnetic Flux (T)', 'Location', 'NorthEast' );
xlabel( 'Frequency (kHz)' );
ylabel( 'Magnetic Flux (T)' );
zlabel( 'Volumetric Losses (kW/m^3)' );
grid on
view( -33.9, 13.6 );
set(gca, 'XScale', 'log');


