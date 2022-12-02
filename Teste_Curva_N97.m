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
title ('Core Losses X Frequency - N97');


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
title ('Core Losses X Frequency - N97');


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
title ('Core Losses X Frequency - N97');


%% Variar frequência --> Encontrar Pv

% Definir range de frequências
%Target.f_min = 100e3;	% Min frequency - Hz
%Target.f_max = 300e3;	% Max frequency - Hz

Target.f_min = 100;	% Min frequency - kHz
Target.f_max = 300;	% Max frequency - kHz

% Encontrar perdas volumétricas para o range de frequências
Range_f = Target.f_min:50:Target.f_max; %igualmente espaçado de 5 em 5
%Range_f = Target.f_min:Target.f_max;

for i = Target.f_min:Target.f_max
    Range_Pv = interp1(x2,y2,Range_f);
    i = i+1;
end

    
%% Cálculo de Perdas no Núcleo (Pv x Ve)

% CORE DATABASE %

Core_Table = readtable('Cores_New.xlsx', 'Sheet', 1);

% Save values to new variables
Cores.Name = string(Core_Table.Core');          % String
% Cores.AeAw = Core_Table.AeAw_mm_4_' * 1e-12;    % m^4    
% Cores.Ae = Core_Table.Ae_mm_2_' * 1e-6;         % m^2
% Cores.Aw = Core_Table.Aw_mm_2_' * 1e-6;         % m^2;
% Cores.le = Core_Table.le_mm_' * 1e-3;           % m;
% Cores.lt = Core_Table.lt_mm_' * 1e-3;           % m
Cores.Ve = Core_Table.Ve_mm_3_' * 1e-9;           % m^3

% Clear data
clear Core_Table;

%% Cálculo Perdas Volumétricas para frequência fixa

f_fixa = 100;
Pv = interp1(x2,y2,f_fixa)*1e3; % kW/m^3

% Criar matriz para armazenar valores
Loss.Values = zeros(1, size(Cores.Ve, 2));

% Perda Volumétrica fixa para uma frequência
% Variando apenas o volume do núcleo de acordo com tabela
Loss.Values = Pv*Cores.Ve;


%% Cálculo Perdas Volumétricas para frequência variando

% for i = Target.f_min:Target.f_max
%     Range_Pv = interp1(x,y,Range_f)*1e3;
%     Core.Loss = Range_Pv'*Cores.Ve;
%     i = i+1;
% end

% Criar matriz para armazenar valores
Core.Loss = zeros(1, size(Cores.Ve, 2));

for idx_f = 1:numel(Range_f)
 f_now = Range_f(idx_f);
 Range_Pv = interp1(x2,y2,Range_f);

 for idx_size = 1:numel(Cores.Ve)
     %Range_Pv = interp1(x,y,f_now);
     Core.Loss = Range_Pv*Cores.Ve(idx_size);
     
 end
end












