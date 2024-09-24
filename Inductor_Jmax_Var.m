clear all;
close all;
warning ('off','all');
clc;

%% Add subfolders to search path

addpath('Datasources');

%% Inductor parameters

Target.L = 550e-6;          % Inductance - H
Target.fmin = 21e3;         % Minimum frequency - Hz (For losses calculation)
Target.fmax = 100e3;        % Maximum frequency - Hz (For wire calculation)
Target.I_rms = 2;         % RMS current - A
Target.I_pk = 3;          % Max current - A
Target.kw_min = 0.50;       % Max window utilization factor - Unitless
Target.kw_max = 0.70;       % Max window utilization factor - Unitless
Target.B_min = 0.20;        % Min peak flux density - T
Target.B_max = 0.35;        % Max peak flux density - T
Target.J_min = 100e4;       % Min current density - A/m^2
Target.J_max = 450e4;       % Max current density - A/m^2

Param.Standard_Wire = 0;    % Use only wire data from catalog (no custom strands)

Param.kh = 40;              % Hysteresis losses constant - W*s*(T^-2.4)/m^3
Param.kf = 4e-5;            % Foucault losses constant - W*s*(T^-2.4)/m^3

% Print message
fprintf("Inductor parameters: \n");
fprintf("\t L: %.2fuH @ %.2fkHz \n", Target.L/1e-6, Target.fmin/1e3);
fprintf("\t Irms: %.2fA, Ipk: %.2fA \n", Target.I_rms, Target.I_pk);
fprintf("\t kW between %d%% and %d%% \n", Target.kw_min*100, Target.kw_max*100);
fprintf("\t Bmax between %.2fT and %.2fT @ Ipk \n", Target.B_min, Target.B_max);
fprintf("\t Jmax between %dA/cm^2 and %dA/cm^2 @ Irms \n", Target.J_min/1e4, Target.J_max/1e4);

%% Core database

% Print message
fprintf("Step 1 - Loading cores database... ");

% Load core data
try
    Table_Cores = readtable('Cores.xlsx', 'Sheet', 1, 'TreatAsEmpty', {'' '.'});
    Table_Cores.Properties.VariableNames = {'Name', 'AeAw', 'Ae', 'Aw', 'le', 'lt', 'Ve', 'Datasheet'};
    Table_Cores.Name = string(Table_Cores.Name);    % String
    Table_Cores.AeAw = Table_Cores.AeAw * 1e-12;    % m^4    
    Table_Cores.Ae = Table_Cores.Ae * 1e-6;         % m^2
    Table_Cores.Aw = Table_Cores.Aw * 1e-6;         % m^2;
    Table_Cores.le = Table_Cores.le * 1e-3;         % m;
    Table_Cores.lt = Table_Cores.lt * 1e-3;         % m
    Table_Cores.Ve = Table_Cores.Ve * 1e-9;         % m^3
catch
    fprintf("\n\t Error loading file. Aborting");
    return;
end

% Sort table by AwAw value - Ascending order
Temp_Idx = find(strcmp(Table_Cores.Properties.VariableNames(:), 'AeAw'));
Table_Cores = sortrows(Table_Cores, Temp_Idx, 'ascend');

% Remove entries with no AeAw info
Temp_Idx = (Table_Cores.AeAw(:) == 0);
Table_Cores(Temp_Idx, :) = [];

% Print message
fprintf("Done.\n\t %d cores will be considered. \n", numel(Table_Cores.Name));

% Clear data
clear -regexp ^Temp_ ^Idx_;

%% Wire database

% Print message
fprintf("Step 2 - Loading wires database... ");

% Load wire data
try
    Table_Wires = readtable('Wires.xlsx', 'Sheet', 1, 'TreatAsEmpty', {'' '.'});
    Table_Wires.Properties.VariableNames = {'Name', 'Strands', 'AWG', 'S_Cu', 'S_Total', 'Ohm_m', 'Custom'};
    Table_Wires.Name = string(Table_Wires.Name);    % String
catch
    fprintf("\n\t Error loading file. Aborting");
    return;
end

% % Sort table by S_Cu value - Descending order
% Temp_Idx = find(strcmp(Table_Wires.Properties.VariableNames(:), 'S_Cu'));
% Table_Wires = sortrows(Table_Wires, Temp_Idx, 'descend');

% Print message
fprintf("Done.\n\t %d wires will be considered. \n", numel(Table_Wires.AWG));

% Clear data
clear -regexp ^Temp_ ^Idx_;

%% Wire section selection

% Print message
fprintf("Step 3 - Finding valid litz wires... ");

% Remove custom strands if necessary
if (Param.Standard_Wire == 1)
    Temp_Idx = (Table_Wires.Custom == 0);
    Res_Wires = Table_Wires(Temp_Idx, :);
else 
    Res_Wires = Table_Wires;
end

% Skin effect penetration
Temp_S_skin = pi*(7.5e-2^2)/Target.fmax;

% Find all wires that can be used due to skin effect
Temp_Idx = (Res_Wires.S_Cu./Res_Wires.Strands) <= Temp_S_skin;

% Save results
Res_Wires = Res_Wires(Temp_Idx, :);

% Required copper section - Limits
Temp_S_Cu_min = Target.I_rms/Target.J_max;
if (Target.J_min > 0)
    Temp_S_Cu_max = Target.I_rms/Target.J_min;
else
    Temp_S_Cu_max = 1;
end

% Find all wires that can be used due to current density
Temp_Idx = (Res_Wires.S_Cu <= Temp_S_Cu_max) & (Res_Wires.S_Cu >= Temp_S_Cu_min);

% No valid wire configuration was found
if (sum(Temp_Idx) == 0)
    fprintf ("\n\t No valid wire configuration was found. Please check the design constrains.");
    return;
end

% Save results
Res_Wires = Res_Wires(Temp_Idx, :);

% Save current density in wires table
Res_J = Target.I_rms./Res_Wires.S_Cu;

% Print message
fprintf("Done.\n\t %d possible wiring configurations will be tested. \n", sum(Temp_Idx,'all'));

% Clear data
clear -regexp ^Temp_ ^Idx_;

%% Number of turns selection

% Print message
fprintf("Step 5 - Finding adequate cores... ");

% Absolute minimum number of turns - Largest core
Temp_Turns_Abs_min = ceil((Target.L*Target.I_pk)/(Target.B_max*max(Table_Cores.Ae)));

% Absolute maximum number of turns - Smallest core
Temp_Turns_Abs_max = ceil((Target.L*Target.I_pk)/(Target.B_max*min(Table_Cores.Ae)));

% Restrain absolutes if necessary
if (Temp_Turns_Abs_min < 1)
    Temp_Turns_Abs_min = 1;
end

% Save results
Res_Turns = Temp_Turns_Abs_min:Temp_Turns_Abs_max;

% Rows are cores, columns are turns
Temp_Configurations_Dimensions = [numel(Table_Cores.Name), numel(Res_Turns)];
Temp_Idx_Max = prod(Temp_Configurations_Dimensions);

% Preallocate matrices to get results
% Res_Bpk(Cores, Turns)
Res_Bpk = ones(Temp_Configurations_Dimensions)*NaN;

% Check the maximum flux for each core and number of turns
for Idx_Cfg = 1:Temp_Idx_Max
    % Get indexes
    [Idx_Core, Idx_Turns] = ind2sub(Temp_Configurations_Dimensions,Idx_Cfg);
    
    % Calculate peak flux
    Temp_B = (Target.L*Target.I_pk)/(Res_Turns(Idx_Turns)*Table_Cores.Ae(Idx_Core));

    % Save result if Bpk is valid
    if ((Temp_B >= Target.B_min) && (Temp_B <= Target.B_max))
      Res_Bpk (Idx_Core, Idx_Turns) = Temp_B;
    end    
end

% No valid core was found
if (sum(~isnan(Res_Bpk),'all') == 0)
    fprintf ("\n\t No valid core was found. Please check the design constrains.");
    return;
end

% Find turn numbers with no valid configuration
Temp_Idx = all(isnan(Res_Bpk), 1);

% Remove invalids from result
if (sum(Temp_Idx, 'all') ~= 0)
    Res_Bpk (:, Temp_Idx)= [];
    Res_Turns(Temp_Idx) = [];
end

% Find cores with no valid configuration
Temp_Idx = all(isnan(Res_Bpk), 2);

% Remove invalids from result
if (sum(Temp_Idx, 'all') ~= 0)
    Res_Bpk (Temp_Idx, :)= [];
end

% Save valid cores
Res_Cores = Table_Cores(~Temp_Idx, :);

% Print message
fprintf("Done.\n\t %d cores will be tested, resulting in %d Core vs Turns combinations to be analyzed. \n", numel(Res_Cores.Name), sum(~isnan(Res_Bpk),'all'));
fprintf("\t %d different turn counts are valid to at least one core. \n", numel(Res_Turns));

% Clear data
clear -regexp ^Temp_ ^Idx_;

%% Valid combinations

% Print message
fprintf("Step 6 - Finding adequate combinations of cores, wires and parallel conductors... ");

% Rows are cores, columns are turns, 3rd dimension is wires
Temp_Configurations_Dimensions = [numel(Res_Cores.Name), numel(Res_Turns), numel(Res_Wires.Name)];
Temp_Idx_Max = prod(Temp_Configurations_Dimensions);

% Preallocate matrices to receive the results
Res_kw = ones(Temp_Configurations_Dimensions)*NaN;

% Find winding configurations that fit in the core Aw
for Idx_Cfg = 1:Temp_Idx_Max
    % Get indexes
    [Idx_Core, Idx_Turns, Idx_Wire] = ind2sub(Temp_Configurations_Dimensions,Idx_Cfg);

    % Selected number of turns for selected core is valid
    if (~isnan(Res_Bpk(Idx_Core, Idx_Turns)))
        % Calculate required Aw to allocate windind
        Temp_Aw_min = (Res_Turns(Idx_Turns)*Res_Wires.S_Total(Idx_Wire));

        % Calculate window usage factor
        Temp_kw = Temp_Aw_min/Res_Cores.Aw(Idx_Core);

        % Save result if usage factor is valid
        if ((Temp_kw >= Target.kw_min) && (Temp_kw <= Target.kw_max))
            Res_kw(Idx_Core, Idx_Turns, Idx_Wire) = Temp_kw;
        end
    end
end

% No valid combination was found
if (sum(~isnan(Res_kw),'all') == 0)
    fprintf ("\n\t No valid combination was found. Please check the design constrains.");
    return;
end

% Remove invalid cores from the results
for Idx_Core = numel(Res_Cores.Name):-1:1
    if all(isnan(Res_kw(Idx_Core, :, :)), 'all')
        Res_kw(Idx_Core, :, :) = [];
        Res_Bpk (Idx_Core, :) = [];
        Res_Cores(Idx_Core, :) = [];
    end
end

% Remove invalid turns from the results
for Idx_Turns = size(Res_Turns, 2):-1:1
    if all(isnan(Res_kw(:, Idx_Turns, :)), 'all')
        Res_kw(:, Idx_Turns, :) = [];
        Res_Bpk (:, Idx_Turns) = [];
        Res_Turns(Idx_Turns) = [];
    end
end

% Remove invalid wires from the results
for Idx_Wire = numel(Res_Wires.Name):-1:1
    if all(isnan(Res_kw(:, :, Idx_Wire)), 'all')
        Res_kw(:, :, Idx_Wire) = [];
        Res_Wires(Idx_Wire) = [];
    end
end

% Print message
fprintf("Done.\n\t %d combinations are adequate and will be tested. \n", sum(~isnan(Res_kw),'all'));
fprintf("\t %d different cores, %d different turn counts and %d wirings are valid. \n", numel(Res_Cores.Name), numel(Res_Turns), numel(Res_Wires.Name));

% Clear data
clear -regexp ^Temp_ ^Idx_;

%% Losses calculation

% Print message
fprintf("Step 7 - Calculating losses... ");

% Rows are cores, columns are turns, 3rd dimension is wires
Temp_Configurations_Dimensions = [numel(Res_Cores.Name), numel(Res_Turns), numel(Res_Wires.Name)];
Temp_Idx_Max = prod(Temp_Configurations_Dimensions);

% Preallocate matrices to receive the results
% Res_P (Cores, Turns, Wires)
Res_P_Wire = ones(Temp_Configurations_Dimensions)*NaN;
Res_P_Core = ones(Temp_Configurations_Dimensions)*NaN;

% Calculate losses for each configuration
for Idx_Cfg = 1:Temp_Idx_Max
    % Get indexes
    [Idx_Core, Idx_Turns, Idx_Wire] = ind2sub(Temp_Configurations_Dimensions,Idx_Cfg);

    % Selected configuration is valid
    if (~isnan(Res_kw(Idx_Core, Idx_Turns, Idx_Wire)))    
        % Calculate losses
        Res_P_Wire(Idx_Core, Idx_Turns, Idx_Wire) = (Res_Turns(Idx_Turns) * Res_Cores.lt(Idx_Core) * Res_Wires.Ohm_m(Idx_Wire))*Target.I_rms^2;
        Res_P_Core(Idx_Core, Idx_Turns, Idx_Wire) = (Res_Bpk(Idx_Core, Idx_Turns)^2.4)*(Param.kh * Target.fmin + Param.kf *Target.fmin)*Res_Cores.Ve(Idx_Core);
    end
end

% Total losses
Res_P_Total = Res_P_Wire + Res_P_Core;

% Choose best projet based on losses
[Temp_Min_Losses, Temp_Idx] = min(Res_P_Total(:));
[Idx_Core, Idx_Turns, Idx_Wire] = ind2sub(size(Res_P_Total), Temp_Idx);

% Print message
fprintf("Done.\n\" );
fprintf ("\t Best project: %s - %d turns - %s \n", Res_Cores.Name(Idx_Core), Res_Turns(Idx_Turns), Res_Wires.Name(Idx_Wire));
fprintf ("\t J %.1fA/cm^2 - Bmax %.3fT - kw %.1f%% - Losses %.3fW \n", Res_J(Idx_Wire)/1e4, Res_Bpk(Idx_Core, Idx_Turns), Res_kw(Idx_Core, Idx_Turns, Idx_Wire)*100, Temp_Min_Losses);

% Clear data
clear -regexp ^Temp_ ^Idx_;

%% Plot all results - J

Temp_hFig = figure ();
set(Temp_hFig, 'units', 'normalized', 'InnerPosition',[0 0 1 1]);

bar(Res_J*1e-4);
set(gca, 'XTickLabel', Res_Wires.Name);
title('Current density');
ylabel('Current density (A/cm^2)');

% Clear data
clear -regexp ^Temp_ ^Idx_;

%% Plot all results - Bpk

Temp_hFig = figure ();
set(Temp_hFig, 'units', 'normalized', 'InnerPosition',[0 0 1 1]);

% Res_Bpk(Cores, Turns)

Temp_h = heatmap(Res_Turns, Res_Cores.Name, Res_Bpk, 'Colormap', jet);
Temp_h.MissingDataColor = [1 1 1];
Temp_h.MissingDataLabel = 'Inv'; 
Temp_h.GridVisible = 'off';

Temp_h.Title = 'Max flux density (T)';
Temp_h.XLabel = 'Turns';
Temp_h.YLabel = 'Core';

% Clear data
clear -regexp ^Temp_ ^Idx_;

%% Plot all results - Window utilization

Temp_hFig = figure ();
set(Temp_hFig, 'units', 'normalized', 'InnerPosition',[0 0 1 1]);
sgtitle ("Window utilization");

% Plot only in columns 1-3. Column 4 if for legend
% Define the right number of rows and columns for subplots
if numel(Res_Cores.Name) < 4
    Temp_Rows = 1;
    Temp_Columns = numel(Res_Cores.Name) + 1;
    
    % Subplots to receive legend
    Temp_Legend_Space = Temp_Columns;
else
    Temp_Rows = ceil(numel(Res_Cores.Name)/3);
    Temp_Columns = 4;
    
    % Subplots to receive legend
    Temp_Legend_Space = (1:Temp_Rows)*4;
end

% Subplots to receive results
Temp_Plot_Space = setdiff(1:(numel(Res_Cores.Name) + numel(Temp_Legend_Space)),Temp_Legend_Space);

% Plot each core in one subplot
for Idx_Core = 1:numel(Res_Cores.Name)
    subplot(Temp_Rows, Temp_Columns, Temp_Plot_Space(Idx_Core));
    hold on;
    grid on;
    xlabel('Turns');
    ylabel('kw (%)'); 
    title (Res_Cores.Name(Idx_Core));
    ylim([Target.kw_min*100, Target.kw_max*100]);

    for Idx_Wire = 1:numel(Res_Wires.Name)
        if (sum(~isnan(Res_kw(Idx_Core, :, Idx_Wire)),'all') ~= 0)    
            % Configuration has only one valid turn count - Put a marker
            if (sum(~isnan(Res_kw(Idx_Core, :, Idx_Wire)), 'all') == 1)
                Temp_Marker = '-o';
            else
                Temp_Marker = '-';
            end
            
            % Plot line
            Temp_h = plot (Res_Turns, Res_kw(Idx_Core, :, Idx_Wire)*100, Temp_Marker, 'DisplayName', Res_Wires.Name(Idx_Wire));

            % Add data tip to line
            Temp_h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Trace',repmat({Temp_h.DisplayName},size(Temp_h.XData)));
        end
    end
end

% Merge all rows of 4th column and plot only the legend
subplot(Temp_Rows, Temp_Columns, Temp_Legend_Space);  
plot(1:sum(~isnan(Res_J), 'all'), nan);
axis off;
legend (Res_Wires.Name, 'Location', 'east', 'NumColumns', 2);

% Clear data
clear -regexp ^Temp_ ^Idx_;

% %% Plot all results - Peak flux
% 
% Temp_hFig = figure ();
% set(Temp_hFig, 'units', 'normalized', 'InnerPosition',[0 0 1 1]);
% sgtitle ("Magnetic flux");
% 
% % Plot only in columns 1-3. Column 4 if for legend
% % Define the right number of rows and columns for subplots
% if numel(Res_Cores.Name) < 4
%     Temp_Rows = 1;
%     Temp_Columns = numel(Res_Cores.Name);
% else
%     Temp_Rows = ceil(numel(Res_Cores.Name)/3);
%     Temp_Columns = 3;
% end
% 
% % Subplots to receive results
% Temp_Plot_Space = 1:(numel(Res_Cores.Name));
% 
% % Plot each core in one subplot
% for Idx_Core = 1:numel(Res_Cores.Name)
%     subplot(Temp_Rows, Temp_Columns, Temp_Plot_Space(Idx_Core));
%     hold on;
%     grid on;
%     xlabel('Turns');
%     ylabel('B (T)'); 
%     title (Res_Cores.Name(Idx_Core));
%     ylim([Target.B_min, Target.B_max]);
% 
%     if (sum(~isnan(Res_Bpk(Idx_Core, :)),'all') ~= 0)    
%         % Configuration has only one valid turn count - Put a marker
%         if (sum(~isnan(Res_Bpk(Idx_Core, :)), 'all') == 1)
%             Temp_Marker = '-o';
%         else
%             Temp_Marker = '-';
%         end
% 
%         % Plot line
%         Temp_h = plot (Res_Turns, Res_Bpk(Idx_Core, :), Temp_Marker);
%     end
% end
% 
% % Clear data
% clear -regexp ^Temp_ ^Idx_;

%% Plot all results - Losses

Temp_hFig = figure ();
set(Temp_hFig, 'units', 'normalized', 'InnerPosition',[0 0 1 1]);
sgtitle ("Total losses");

% Plot only in columns 1-3. Column 4 if for legend
% Define the right number of rows and columns for subplots
if numel(Res_Cores.Name) < 4
    Temp_Rows = 1;
    Temp_Columns = numel(Res_Cores.Name) + 1;
    
    % Subplots to receive legend
    Temp_Legend_Space = Temp_Columns;
else
    Temp_Rows = ceil(numel(Res_Cores.Name)/3);
    Temp_Columns = 4;
    
    % Subplots to receive legend
    Temp_Legend_Space = (1:Temp_Rows)*4;
end

% Subplots to receive results
Temp_Plot_Space = setdiff(1:(numel(Res_Cores.Name) + numel(Temp_Legend_Space)), Temp_Legend_Space);

% Plot each core in one subplot
for Idx_Core = 1:numel(Res_Cores.Name)
    subplot(Temp_Rows, Temp_Columns, Temp_Plot_Space(Idx_Core));
    hold on;
    grid on;
    xlabel('Turns');
    ylabel('Losses (W)'); 
    title (Res_Cores.Name(Idx_Core));
    ylim([min(Res_P_Total(:)) max(Res_P_Total(:))]); % Normalized Y axis

    for Idx_Wire = 1:numel(Res_Wires.Name)
        if (sum(~isnan(Res_kw(Idx_Core, :, Idx_Wire)),'all') ~= 0)    
            % Configuration has only one valid turn count - Put a marker
            if (sum(~isnan(Res_kw(Idx_Core, :, Idx_Wire)), 'all') == 1)
                Temp_Marker = '-o';
            else
                Temp_Marker = '-';
            end
            
            % Plot line
            Temp_h = plot (Res_Turns, Res_P_Total(Idx_Core, :, Idx_Wire), Temp_Marker, 'DisplayName', Res_Wires.Name(Idx_Wire));            
            
            % Add data tip to line
            Temp_h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Trace',repmat({Temp_h.DisplayName},size(Temp_h.XData)));
        end
    end
end

% Merge all rows of 4th column and plot only the legend
subplot(Temp_Rows, Temp_Columns, Temp_Legend_Space);  
plot(1:sum(~isnan(Res_J), 'all'), nan);
axis off;
legend (Res_Wires.Name, 'Location', 'east', 'NumColumns', 2);

% Clear data
clear -regexp ^Temp_ ^Idx_;

%% Plot best result

% Choose best projet based on losses
[Temp_Min_Losses, Temp_Idx] = min(Res_P_Total(:));
[Idx_Core, Idx_Turns, Idx_Wire] = ind2sub(size(Res_P_Total), Temp_Idx);

Temp_hFig = figure ();
set(Temp_hFig, 'units', 'normalized', 'InnerPosition',[0 0 1 1]);
sgtitle (sprintf ("Best core: %s", Res_Cores.Name(Idx_Core)));
hold on;
grid on;
xlabel('Turns');
ylabel('Total losses (W)');

% Plot best point
Temp_h = plot (Res_Turns(Idx_Turns), Res_P_Total(Idx_Core, Idx_Turns, Idx_Wire), '-o', 'HandleVisibility','off');     

% Array to receive wiring labels
Temp_Label = strings(numel(Res_Wires.Name), 1);

for Idx_Wire = 1:numel(Res_Wires.Name)
    if (sum(~isnan(Res_kw(Idx_Core, :, Idx_Wire)),'all') ~= 0)    
        % Configuration has only one valid turn count - Put a marker
        if (sum(~isnan(Res_kw(Idx_Core, :, Idx_Wire)), 'all') == 1)
            Temp_Marker = '-o';
        else
            Temp_Marker = '-';
        end

        % Create line label
        Temp_Label(Idx_Wire) = Res_Wires.Name(Idx_Wire);

        % Plot line
        Temp_h = plot (Res_Turns, Res_P_Total(Idx_Core, :, Idx_Wire), Temp_Marker, 'DisplayName', Res_Wires.Name(Idx_Wire));            

        % Add data tip to line
        Temp_h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Trace',repmat({Temp_h.DisplayName},size(Temp_h.XData)));
    end
end

% Remove empty labels
Temp_Label = Temp_Label(~strcmp(Temp_Label(:), ""));

legend (Temp_Label, 'Location', 'best', 'NumColumns', 2);

% Clear data
clear -regexp ^Temp_ ^Idx_;