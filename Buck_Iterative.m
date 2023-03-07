clear all;
close all;
warning ('off','all');
clc;

%% Add subfolders to search path

addpath('Functions', 'Datasources');

%% Design parameters

% Frequency range
Param.Fsw_Range = 50e3:5e3:350e3;

% Switching period range
Param.Tsw_Range = Param.Fsw_Range.^-1;

% Time resolution to integrate waveforms and calculate power loss (seconds)
Param.Time_Step = 1e-10;

%% Converter parameters

Param.Po = 180;             % Nominal output power

Param.Vi_min = 400;         % Min input voltage
Param.Vi_max = 400;         % Max input voltage
Param.Delta_Vi = 0.05;      % Max high frequency input voltage ripple (%)

Param.Vo_min = 100;         % Min output voltage
Param.Vo_max = 100;         % Max output voltage
Param.Delta_Vo = 0.01;      % Max output voltage ripple (%)

Param.Delta_IL = 0.8;       % Max inductor current ripple (%)

%% Switch parameters - GS-065-018-2-L

S.RDS = 78e-3;     % RDS on
S.Rg = 1.3;        % Gate resistance
S.Coss = 34e-12;   % Output capacitance
S.Crss = 0.4e-12;  % Reverse capacitance
S.Qg = 4e-9;       % Total gate charge
S.Qth = 0.8e-9;    % Gate charge at threshold
S.Qgs = 1.2e-9;    % Gate to source charge
S.Qgd = 1.2e-9;    % Gate to drain charge
S.Vpl = 3.5;       % Plateau voltage
S.Vsd = 2.0;       % Reverse conduction voltage drop
S.Qrr = 0;         % Reverse recovery charge
S.kt = 0.0105;     % RDSon temperature coefficient (1/oC)

%% Gate driver parameters - LM5113

Drv.Vgs = 5;       % Driver voltage
Drv.Td = 150e-9;   % Dead time
Drv.Ron = 3.3;     % Gate on resistance
Drv.Roff = 0;      % Gate off resistance
Drv.Rcon = 2.1;    % IC internal on resistance
Drv.Rcoff = 0.6;   % IC internal off resistance
Drv.Iq = 2.2e-3;   % IC operating current (IDDo + IDD + IHB)

%% Inductor parameters

Inductor.kw_min = 0.20;	% Max window utilization factor - Unitless
Inductor.kw_max = 0.70;	% Max window utilization factor - Unitless
Inductor.B_min = 0.00;	% Min peak flux density - T
Inductor.B_max = 0.30;	% Max peak flux density - T
Inductor.J_max = 450e4;	% Max current density - A/m^2

%% Fixed parameters

% Duty cycle
Param.D_min = Param.Vo_min/Param.Vi_min;
Param.D_max = Param.Vo_max/Param.Vi_max;

% Nominal inductor currents
Param.IL_avg_min = Param.Po/Param.Vo_max;
Param.IL_avg_max = Param.Po/Param.Vo_min;

% Inductor currents for component design -> Safety margin
Inductor.I_rms = Param.IL_avg_max*1.25;       
Inductor.I_pk = (Param.IL_avg_max + Param.IL_avg_max*Param.Delta_IL/2)*1.25;

% Gate currents
Drv.Ion = (Drv.Vgs - S.Vpl)/(Drv.Ron + Drv.Rcon + S.Rg);
Drv.Ioff = (Drv.Vgs - S.Vpl)/(Drv.Roff + Drv.Rcoff + S.Rg);

% Switching times
Drv.Ton = (S.Qgs - S.Qth + S.Qgd)/Drv.Ion;
Drv.Toff = (S.Qgs - S.Qth + S.Qgd)/Drv.Ioff;

%% Converter design

% Prealocate variables to store results
Res_Ci = ones(1, numel(Param.Fsw_Range))*NaN;     % Input capacitance value
Res_Co = ones(1, numel(Param.Fsw_Range))*NaN;     % Output capacitance value
Res_L_Ind = ones(1, numel(Param.Fsw_Range))*NaN;  % Inductor value
Res_L_Core = strings(1, numel(Param.Fsw_Range));  % Inductor core
Res_L_Turns = ones(1, numel(Param.Fsw_Range))*NaN;% Inductor winding turns
Res_L_Cond = ones(1, numel(Param.Fsw_Range))*NaN; % Inductor winding conductors
Res_L_AWG = ones(1, numel(Param.Fsw_Range))*NaN;  % Inductor winding AWG
Res_L_kw = ones(1, numel(Param.Fsw_Range))*NaN;   % Inductor window utilization factor
Res_P_Lw = ones(1, numel(Param.Fsw_Range))*NaN;   % Inductor wire losses
Res_P_Lc = ones(1, numel(Param.Fsw_Range))*NaN;   % Inductor core losses
Res_P_Lt = ones(1, numel(Param.Fsw_Range))*NaN;   % Inductor total losses
    
% Define converter parameters for each frequency
for Idx_Freq = 1:numel(Param.Fsw_Range)
    % Calculate required capacitances
    Res_Ci(Idx_Freq) = (Param.IL_avg_max*(1-Param.D_min)*Param.D_min)/(Param.Vi_min*Param.Delta_Vi*Param.Fsw_Range(Idx_Freq));
    Res_Co(Idx_Freq) = (Param.IL_avg_max*Param.Delta_IL)/(8*Param.Vo_min*Param.Delta_Vo*Param.Fsw_Range(Idx_Freq));

    % Define remaining inductor parameters
    Inductor.f = Param.Fsw_Range(Idx_Freq);
    Inductor.L =(Param.D_max*(Param.Vi_max-Param.Vo_max))/(Param.IL_avg_min*Param.Delta_IL*Param.Fsw_Range(Idx_Freq));
    
    % Design inductor
    Temp_Result = Inductor_Design (Inductor);

    % Save results
    Res_L_Ind(Idx_Freq) = Inductor.L;
    Res_L_Core(Idx_Freq) = Temp_Result.Core;
    Res_L_Turns(Idx_Freq) = Temp_Result.Turns;
    Res_L_Cond(Idx_Freq) =  Temp_Result.Cond;
    Res_L_AWG(Idx_Freq) = Temp_Result.AWG;
    Res_L_kw(Idx_Freq) = Temp_Result.kw;
    Res_P_Lw(Idx_Freq) = Temp_Result.Pw; 
    Res_P_Lc(Idx_Freq) = Temp_Result.Pc;
    Res_P_Lt(Idx_Freq) = Temp_Result.Pt;
end

% Clear data
clear -regexp ^Temp_;

%% Calculate power losses

% Prealocate variables to store power losses
Res_P_Smc = ones(1, numel(Param.Fsw_Range))*NaN;        % Main switch conduction losses
Res_P_Smsw = ones(1, numel(Param.Fsw_Range))*NaN;       % Main switch switching losses
Res_P_Smt = ones(1, numel(Param.Fsw_Range))*NaN;        % Main switch total losses
Res_P_Sfwc = ones(1, numel(Param.Fsw_Range))*NaN;       % Free-wheeling switch conduction losses
Res_P_Sfwsw = ones(1, numel(Param.Fsw_Range))*NaN;      % Free-wheeling switch switching losses
Res_P_Sfwt = ones(1, numel(Param.Fsw_Range))*NaN;       % Free-wheeling switch total losses
Res_P_Aux = ones(1, numel(Param.Fsw_Range))*NaN;        % Auxiliary circuits losses
Res_Pt = ones(1, numel(Param.Fsw_Range))*NaN;           % Total losses
Res_Efficiency = ones(1, numel(Param.Fsw_Range))*NaN;   % Efficiency
    
% Define converter losses for each frequency
for Idx_Freq = 1:numel(Param.Fsw_Range)
    % Time array
    Temp_t = 0:Param.Time_Step:Param.Tsw_Range(Idx_Freq);
  
    % Inductor currents - Limits
    Temp_Delta_I_L = (Param.Vi_min-Param.Vo_min)*Param.D_min*Param.Tsw_Range(Idx_Freq)/Res_L_Ind(Idx_Freq);
    Temp_IL_min = (Param.IL_avg_max - Temp_Delta_I_L/2);
    Temp_IL_max = (Param.IL_avg_max + Temp_Delta_I_L/2);
    
    % Inductor currents - Array
    Temp_I_L = zeros (size(Temp_t));
    Temp_Interval = logical(Temp_t<=Param.D_min*Param.Tsw_Range(Idx_Freq));
    Temp_I_L(Temp_Interval) = (Param.Vi_min-Param.Vo_min)*Temp_t(Temp_Interval)/Res_L_Ind(Idx_Freq) + Temp_IL_min;
    Temp_Interval = logical(Temp_t>Param.D_min*Param.Tsw_Range(Idx_Freq));
    Temp_I_L(Temp_Interval) = (-Param.Vo_min*(Temp_t(Temp_Interval) - Param.D_min*Param.Tsw_Range(Idx_Freq)))/Res_L_Ind(Idx_Freq) + Temp_IL_max;
    Temp_I_L_rms = rms(Temp_I_L);
    
    % Main switch currents
    Temp_I_Sm = zeros (size(Temp_t));
    Temp_Interval = logical(Temp_t<=Param.D_min*Param.Tsw_Range(Idx_Freq));
    Temp_I_Sm(Temp_Interval) = Temp_I_L(Temp_Interval);
    Temp_I_Sm_rms = rms(Temp_I_Sm);

    % Free wheeling switch currents - Channel conduction
    Temp_I_Sfwc = zeros (size(Temp_t));
    Temp_Interval = logical((Temp_t >=(Param.D_min*Param.Tsw_Range(Idx_Freq) + Drv.Td)).*(Temp_t < (Param.Tsw_Range(Idx_Freq) - Drv.Td)));
    Temp_I_Sfwc(Temp_Interval) = Temp_I_L(Temp_Interval);  
    Temp_I_Sfwc_rms = rms(Temp_I_Sfwc);
    
    % Free wheeling switch currents - Diode conduction
    Temp_I_Sfwd = Temp_I_L - Temp_I_Sm - Temp_I_Sfwc;
    Temp_I_Sfwd_avg = mean(Temp_I_Sfwd);   
      
    % Main switch losses

    % Main switch conduction losses
    Res_P_Smc(Idx_Freq) = Temp_I_Sm_rms^2 * S.RDS;

    % Main switch Coss losses
    Temp_PSm_Coss = 1/2 * (S.Coss - S.Crss)*Param.Vi_min^2 * Param.Fsw_Range(Idx_Freq);

    % Main switch turn on losses
    Temp_PSm_on = Param.Fsw_Range(Idx_Freq)*Param.Vi_min*Temp_IL_min*Drv.Ton / 2;

    % Main switch turn off losses
    Temp_PSm_off = Param.Fsw_Range(Idx_Freq)*Param.Vi_min*Temp_IL_max*Drv.Toff / 2;

    % Main switch switching losses
    Res_P_Smsw(Idx_Freq) = Temp_PSm_Coss + Temp_PSm_on + Temp_PSm_off;

    % Main switch total losses
    Res_P_Smt(Idx_Freq) = Res_P_Smc(Idx_Freq) + Res_P_Smsw(Idx_Freq);

    % Free wheeling switch losses

    % Free wheeling switch conduction losses (Channel + Diode)
    Res_P_Sfwc(Idx_Freq) = Temp_I_Sfwc_rms^2 * S.RDS + S.Vsd * Temp_I_Sfwd_avg;

    % Free wheeling switch Coss losses
    Temp_PSfw_Coss = 1/2 * (S.Coss - S.Crss)*Param.Vi_min^2 * Param.Fsw_Range(Idx_Freq);

    % Free wheeling switch reverse recovery losses
    Temp_PSfw_Qrr = S.Qrr * Param.Vi_min * Param.Fsw_Range(Idx_Freq);

    % Free wheeling switch switching losses
    Res_P_Sfwsw(Idx_Freq) = Temp_PSfw_Coss + Temp_PSfw_Qrr;

    % Free wheeling switch total losses
    Res_P_Sfwt(Idx_Freq) = Res_P_Sfwc(Idx_Freq) + Res_P_Sfwsw(Idx_Freq);       

    % Auxiliary circuits losses

    % Main switch driving current
    Temp_PAux_Sm_drv = Drv.Vgs*S.Qg*Param.Fsw_Range(Idx_Freq);

    % Free wheeling switch driving current
    Temp_PAux_Sfw_drv = Drv.Vgs*S.Qg*Param.Fsw_Range(Idx_Freq);
    
    % Total auxiliary losses
    Res_P_Aux(Idx_Freq) = Temp_PAux_Sm_drv + Temp_PAux_Sfw_drv;
    
    % Total losses
    Res_Pt(Idx_Freq) = Res_P_Lt(Idx_Freq) + Res_P_Smt(Idx_Freq) + Res_P_Sfwt(Idx_Freq) + Res_P_Aux(Idx_Freq);
    
    % Efficiency
    Res_Efficiency(Idx_Freq) = Param.Po/(Param.Po + Res_Pt(Idx_Freq));
end

% Clear data
clear -regexp ^Temp_;

%% Plot component values

Temp_hFig = figure ();
set(Temp_hFig, 'units', 'normalized', 'InnerPosition',[0 0 1 1]);

subplot 121;
hold on;
grid on;
plot(Param.Fsw_Range/1e3, Res_L_Ind/1e-6);
xlabel('Frequency (kHz)');
ylabel('Inductance (uH)'); 

subplot 122;
hold on;
grid on;
plot(Param.Fsw_Range/1e3, Res_Ci/1e-6);
plot(Param.Fsw_Range/1e3, Res_Co/1e-6);
xlabel('Frequency (kHz)');
ylabel('Capacitance (uF)'); 
legend(["Input capacitor", "Output capacitor"]);

% Clear data
clear -regexp ^Temp_;

%% Plot results

Temp_hFig = figure ();
set(Temp_hFig, 'units', 'normalized', 'InnerPosition',[0 0 1 1]);
hold on;
Temp_ymax = ceil(max(Res_Pt)) + 0.5;

% Patches for core labels
Temp_Unique_Cores = unique(Res_L_Core);
for Idx = numel(Temp_Unique_Cores):-1:1
    Temp_Begin = Param.Fsw_Range(find(Res_L_Core == Temp_Unique_Cores(Idx), 1, 'first'))/1e3;
    Temp_End = Param.Fsw_Range(find(Res_L_Core == Temp_Unique_Cores(Idx), 1, 'last'))/1e3;
    Temp_Middle = (Temp_Begin + Temp_End)/2;

    Temp_Color = [0.9 0.9 0.9] - 0.1*mod(Idx,2);
    patch([Temp_Begin Temp_Begin Temp_End Temp_End],[0 Temp_ymax Temp_ymax 0], 'k', 'FaceColor', Temp_Color, 'EdgeColor', 'none', 'HandleVisibility','on');
    text(Temp_Middle, Temp_ymax*0.99, Temp_Unique_Cores(Idx), 'Rotation', 90, 'VerticalAlignment','middle','HorizontalAlignment','right');
end

% All core labels
% for Idx_Freq = 1:numel(Param.Fsw_Range)
%     text(Param.Fsw_Range(Idx_Freq)/1e3, Temp_ymax*0.99, Res_L_Core(Idx_Freq), 'Rotation', 90, 'VerticalAlignment','middle','HorizontalAlignment','right');
% end

% Losses
Temp_b = bar(Param.Fsw_Range/1e3, [Res_P_Aux; Res_P_Sfwc; Res_P_Sfwsw; Res_P_Smc; Res_P_Smsw; Res_P_Lw; Res_P_Lc], 'stacked');
Temp_Legend = ["Other", "Sfw cond", "Sfw sw", "Sm cond", "Sm sw", "L windings", "L core"];
uistack(Temp_b, 'top');

% Efficiency
yyaxis right;
plot(Param.Fsw_Range/1e3, Res_Efficiency, 'k');

% Plot formatting
yyaxis left;
ylabel('Losses (W)'); 
ylim([0, Temp_ymax]);
legend (flip(Temp_b), flip(Temp_Legend), 'Location', 'south', 'Orientation', 'horizontal')

% Get left Tick Values
Temp_ytl = get(gca, 'YTick'); 

yyaxis right;
ylabel('Efficiency (%)');

% Get right Tick Values
Temp_ytr = get(gca, 'YTick');                                    

% Create New Right Tick Values Matching Number Of Left Ticks
Temp_ytrv = linspace(min(Temp_ytr), max(Temp_ytr), numel(Temp_ytl));
% Tick Label Cell Array
Temp_ytrc = compose('%.1f',Temp_ytrv*100);

% Set New Right Tick Labels
set(gca, 'YTick', Temp_ytrv, 'YTickLabel', Temp_ytrc);

grid on;
set(gca, 'Layer', 'top');
Temp_ax = gca;
Temp_ax.YAxis(1).Color = 'k';
Temp_ax.YAxis(2).Color = 'k';
xlabel('Frequency (kHz)');

% Clear data
clear -regexp ^Temp_;