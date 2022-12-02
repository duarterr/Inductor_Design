clc;
clear all;

% Steinmetz Modificado
% Pv = k * feq^(a-1) * Bmax^b * fsw

%% Parâmetros
Bmax = 0.2;                      % T

Pin = 185.56;                    % W
Vin = 400;                       % V
Pout = 180;                      % W
Vout = 103.68;                   % V
fsw = 100;                       % kHz

D = Vout/Vin;

Ve = 7.46e-6;                    % m^3

%% Fit: 'B_0.2T'.

load B2

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x2 = interp1(x_escalalin',x_escalalog',B2(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y2 = interp1(y_escalalin',y_escalalog',B2(:,2));

[xData, yData] = prepareCurveData( x2, y2 );

% Set up fittype and options.
ft = fittype( 'k*x^a*0.2^b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.0172895593763746 0.211907074257636 0.913375856139019];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'B = 200mT' );
h = plot( fitresult, xData, yData );
legend( h, 'Curve B = 200mT', 'Interpolation', 'Location', 'SouthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Frequency (kHz)' );
ylabel( 'Volumetric Losses (kW/m^3)');
grid on


fitresult
%% Cálculo Perdas em W

% Frequência efetiva para forma de onda de excitação não senoidal
feq = (2/pi^2)*(fsw/(D - D^2));

% Coeficientes de Steinmetz extraídos com Fit Curve Tool
r1 = coeffvalues(fitresult);
a1 = r1(1);
b1 = r1(2);
k1 = r1(3);

% Perdas no núcleo 
Pcore = (k1 * (feq^(a1-1)) * (Bmax^b1) * fsw) 



%%
%% Parâmetros do conversor

Pin = 185.56;                    % W
Vin = 400;                       % V
Pout = 180;                      % W
Vout = 103.68;                   % V
fsw = 100;                       % kHz
D = Vout/Vin;
Ve = 7.46e-6;

% Steinmetz Modificado
% Pv = k * feq^(a-1) * Bmax^b * fsw

% Frequência efetiva para forma de onda de excitação não senoidal
feq = (2/pi^2)*(fsw/(D - D^2));

Bmax = 0.2;

load B2

% Eixo x = Frequência
x_escalalog = logspace(1,3,100000);
x_escalalin = linspace(10,1000,100000);
x2 = interp1(x_escalalin',x_escalalog',B2(:,1));

% Eixo y = Perdas volumétricas
y_escalalog = logspace(-1,4,100000);
y_escalalin = linspace(0.1,10000,100000);
y2 = interp1(y_escalalin',y_escalalog',B2(:,2));


Flux = linspace(0.05,0.2,158);
Flux = Flux';

[xData, yData, zData] = prepareSurfaceData( x2, y2, Flux );

% Set up fittype and options.
ft = fittype( 'k*x^a*y^b', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.258034668870745 0.486915030471507 0.421761282626275];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 'untitled fit 1', 'Flux vs. x2, y2', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x2', 'Interpreter', 'none' );
ylabel( 'y2', 'Interpreter', 'none' );
zlabel( 'Flux', 'Interpreter', 'none' );
grid on
view( -13.8, 53.9 );

% Coeficientes de Steinmetz extraídos com Fit Curve Tool
r = coeffvalues(fitresult)
a = r(1);
b = r(2);
k = r(3);

% Perdas no núcleo 
Pcore = (k * (feq^(a-1)) * (Bmax^b) * fsw) 









