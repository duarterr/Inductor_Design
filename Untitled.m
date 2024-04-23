clear all;
close all;
warning ('off','all');
clc;

%% Add subfolders to search path

addpath('Datasources');

% Load wire data
try
    Table_Wires = readtable('Wires_New.xlsx', 'Sheet', 1, 'TreatAsEmpty', {'' '.'});
    Table_Wires.Properties.VariableNames = {'AWG', 'S_Cu', 'S_Total'};
catch
    fprintf("\n\t Error loading file. Aborting");
    return;
end

% Sort table by AWG value - Descending order
Temp_Idx = find(strcmp(Table_Wires.Properties.VariableNames(:), 'AWG'));
Table_Wires = sortrows(Table_Wires, Temp_Idx, 'ascend');

Idx_total = 1;
for Idx_Awg = 1:numel(Table_Wires.AWG)
    for Strand = 1:10
        Strands_New(Idx_total) = Strand;
        AWG_New(Idx_total) = Table_Wires.AWG(Idx_Awg);
        S_Cu_New(Idx_total) = Table_Wires.S_Cu(Idx_Awg)*Strand;
        S_Total_New(Idx_total) = Table_Wires.S_Total(Idx_Awg)*Strand;
        Ohm_m_New(Idx_total) = 1.72E-8/(S_Cu_New(Idx_total));
        
        
        Idx_total = Idx_total + 1;
    end
end

T = table (Strands_New', AWG_New', S_Cu_New', S_Total_New', Ohm_m_New');
T.Properties.VariableNames = {'Strands','AWG', 'S_Cu', 'S_Total', 'Ohm_m'};


writetable(T,'Wires.xlsx','Sheet',1)