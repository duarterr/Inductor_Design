function VL_Core = Volumetric_Losses(Target)

% Função que retorna as perdas volumétricas para uma dada frequência e um
% fluxo magnético máximo, percorrendo a superfície gerada a partir dos
% dados do datasheet do núcleo N97.

% Matriz de 3 dimensões 
% x = frequência (Hz);
% y = fluxo magnético (T);
% z = perdas volumétricas (W/m^3);

M = [100e3 0.025 5e3    ;
     150e3 0.025 7e3    ;
     200e3 0.025 11e3   ;
     300e3 0.025 20e3   ;
     20e3  0.05  4e3    ;
     50e3  0.05  12e3   ;
     80e3  0.05  11e3   ;
     100e3 0.05  27e3   ;
     200e3 0.05  65e3   ;
     20e3  0.1   11e3   ;
     30e3  0.1   12e3   ;
     50e3  0.1   15e3   ;
     80e3  0.1   100e3  ;
     100e3 0.1   120e3  ;
     150e3 0.1   200e3  ;
     200e3 0.1   300e3  ;
     300e3 0.1   600e3  ;
     20e3  0.2   90e3   ;
     40e3  0.2   200e3  ;
     70e3  0.2   400e3  ;
     100e3 0.2   600e3  ;
     150e3 0.2   900e3  ;
     200e3 0.2   1700e3 ;
     300e3 0.2   3000e3 ;
     20e3  0.3   230e3  ;
     30e3  0.3   350e3  ;
     40e3  0.3   470e3  ;
     50e3  0.3   600e3 ]; 
 
 x = M(:,1);
 y = M(:,2);
 z = M(:,3);
 
% Fit Surface - Steinmetz Original Equation (SOE)
[xData, yData, zData] = prepareSurfaceData( x, y, z );
% Set up fittype and options.
ft = fittype( 'k*x^a*y^b', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.817303220653433 0.86869470536351 0.0844358455109103];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Coeficientes de Steinmetz da Superfície
c = coeffvalues(fitresult);
k = c(3);
b = c(2);
a = c(1);

% Cálculo da perda volumétrica
VL_Core = k*Target.f^a*Target.B_max^b;

%Plot fit with data.
figure( 'Name', 'Relative Core Losses Surface' );
h = plot( fitresult, [xData, yData], zData );
%legend( h, 'Surface', 'Data (z vs. x, y)', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Frequency (Hz)', 'Interpreter', 'none' );
ylabel( 'Magnetic Flux (T)', 'Interpreter', 'none' );
zlabel( 'Relative Core Losses (W/m^3)', 'Interpreter', 'none' );
grid on
hold on
scatter3(Target.f,Target.B_max,VL_Core,'filled', 'MarkerFaceColor','r');

end

