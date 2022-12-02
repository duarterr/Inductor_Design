clear all;
clc;


%% Parâmetros do conversor

Pin = 185.56;                    % W
Vin = 400;                       % V
Pout = 180;                      % W
Vout = 103.68;                   % V
fsw = 100;                       % kHz
D = Vout/Vin;
Ve = 7.46e-6;

%% Cálculo de Perdas

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


%% Fit: '200mT'.
[xData, yData, zData] = prepareSurfaceData( x2, Flux, y2 );

% Set up fittype and options.
ft = fittype( 'k*x^a*y^b', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.655740699156587 0.0357116785741896 0.849129305868777];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 'untitled fit 1', 'y2 vs. x2, Flux', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x2', 'Interpreter', 'none' );
ylabel( 'Flux', 'Interpreter', 'none' );
zlabel( 'y2', 'Interpreter', 'none' );
grid on
view( -27.1, 19.9 );

%%
% Coeficientes de Steinmetz extraídos com Fit Curve Tool
r = coeffvalues(fitresult)
a = r(1);
b = r(2);
k = r(3);

% Perdas no núcleo 
Pcore = (k * (feq^(a-1)) * (Bmax^b) * fsw) 
