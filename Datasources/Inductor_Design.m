function Result = Inductor_Design(Target)
    
    % Datasources - Persistent to speedup calculations
    persistent Cores;
    persistent Wires; 
    
    % Check if core database was loaded
    if isempty(Cores)
        % Load core data
        Core_Table = readtable('Cores_New.xlsx', 'Sheet', 1);

        % Save values to new variables
        Cores.Name = string(Core_Table.Core');          % String
        Cores.AeAw = Core_Table.AeAw_mm_4_' * 1e-12;    % m^4    
        Cores.Ae = Core_Table.Ae_mm_2_' * 1e-6;         % m^2
        Cores.Aw = Core_Table.Aw_mm_2_' * 1e-6;         % m^2;
        Cores.le = Core_Table.le_mm_' * 1e-3;           % m;
        Cores.lt = Core_Table.lt_mm_' * 1e-3;           % m
        Cores.Ve = Core_Table.Ve_mm_3_' * 1e-9;         % m^3

        % Crear data
        clear Core_Table;        
    end
    
    % Check if wire database was loaded
    if isempty(Wires)
        % Load wire data
        Wire_Table = readtable('Wires_new.xlsx', 'Sheet', 1);

        % Save values to new variables
        Wires.AWG = Wire_Table.AWG';                    % AWG  
        Wires.S_Cu = Wire_Table.S_Cu_m_2_';             % m^2
        Wires.S_Total = Wire_Table.S_Total_m_2_';       % m^2

        % Crear data
        clear Wire_Table;
    end

    %% Constant parameters
    
    Param.RCu = 1.72e-8;	% Annealed copper resistance - Ohm/m
    %Param.kh = 40;          % Hysteresis losses constant - W*s*(T^-2.4)/m^3
    %Param.kf = 4e-5;        % Foucault losses constant - W*s*(T^-2.4)/m^3
    Param.Turns = 5:500;    % Number of turns to be considered

    %% Wire section selection

    % Required conductor section
    S_Cu_min = Target.I_rms/Target.J_max;

    % Skin effect penetration
    S_skin = pi*(7.5e-2^2)/Target.f;

    % Find all wires that can be used
    Valid_S = Wires.S_Cu<=S_skin;
    Valid_Wires.AWG = Wires.AWG(Valid_S);
    Valid_Wires.S_Cu = Wires.S_Cu(Valid_S);
    Valid_Wires.S_Total = Wires.S_Total(Valid_S);

    % Preallocate matrices to get results
    Valid_Wires.Cond = zeros(1, size(Valid_Wires.AWG, 2));
    Label_Windind = strings(1, size(Valid_Wires.AWG, 2));

    % Find the required number of parallel wires for each valid AWG
    for Idx_Wire = 1:size(Valid_Wires.AWG, 2)
        % Calculate number of required parallel wires for each valid AWG
        Valid_Wires.Cond(Idx_Wire) = ceil(S_Cu_min/Valid_Wires.S_Cu(Idx_Wire));

        % Create labels for each combination
        Label_Windind(Idx_Wire) = sprintf("%d x AWG %d", Valid_Wires.Cond(Idx_Wire), Valid_Wires.AWG(Idx_Wire));
    end

    % Clear data
    clear S_Cu_min S_skin Valid_S Idx_Wire;

    %% Number of turns selection

    Configurations_Dimensions = [size(Cores.Name, 2), size(Param.Turns, 2)];
    Configurations_Max = prod(Configurations_Dimensions);

    % Preallocate matrices to get results
    Bpk = ones(Configurations_Dimensions)*NaN;

    % Check the maximum flux for each core and number of turns
    for Idx_Cfg = 1:Configurations_Max
        % Get indexes
        [Idx_Core, Idx_Turns] = ind2sub(Configurations_Dimensions,Idx_Cfg);

        % Calculate peak flux
        B = (Target.L*Target.I_pk)/(Param.Turns(Idx_Turns)*Cores.Ae(Idx_Core));

        % Save result if Bpk is valid
        if (B >= Target.B_min && B <= Target.B_max)
          Bpk (Idx_Core, Idx_Turns) = B;
        end    
    end

    % Find turn numbers with no valid configuration
    Invalid_Turns = all(isnan(Bpk), 1);

    % Remove invalids from result
    Bpk (:, Invalid_Turns)= [];

    % Save valid turns
    Valid_Turns = Param.Turns(~Invalid_Turns);

    % Find cores numbers with no valid configuration
    Invalid_Cores = all(isnan(Bpk), 2);

    % Remove invalids from result
    Bpk (Invalid_Cores, :)= [];

    % Save valid cores
    Valid_Cores.Name = Cores.Name(~Invalid_Cores);
    Valid_Cores.AeAw = Cores.AeAw(~Invalid_Cores);
    Valid_Cores.Ae = Cores.Ae(~Invalid_Cores);
    Valid_Cores.Aw = Cores.Aw(~Invalid_Cores);
    Valid_Cores.le = Cores.le(~Invalid_Cores);
    Valid_Cores.lt = Cores.lt(~Invalid_Cores);
    Valid_Cores.Ve = Cores.Ve(~Invalid_Cores);

    % Clear data
    clear Configurations_Dimensions Configurations_Max Idx_Cfg Idx_Core Idx_Turns B;
    clear Invalid_Turns Invalid_Cores;

    %% Valid configurations

    Configurations_Dimensions = [size(Valid_Cores.Name, 2), size(Valid_Turns, 2), size(Valid_Wires.AWG, 2)];
    Configurations_Max = prod(Configurations_Dimensions);

    % Preallocate matrices to receive the results
    kw_req = ones(Configurations_Dimensions)*NaN;

    % Find winding configurations that fit in the core Aw
    for Idx_Cfg = 1:Configurations_Max
        % Get indexes
        [Idx_Core, Idx_Turns, Idx_Wire] = ind2sub(Configurations_Dimensions,Idx_Cfg);

        % Selected number of turns for selected core is valid
        if (~isnan(Bpk(Idx_Core, Idx_Turns)))
            % Calculate required Aw to allocate windind
            Aw_min = (Valid_Turns(Idx_Turns)*Valid_Wires.Cond(Idx_Wire)*Valid_Wires.S_Total(Idx_Wire));

            % Calculate window usage factor
            kw = Aw_min/Valid_Cores.Aw(Idx_Core);

            % Save result if usage factor is valid
            if (kw >= Target.kw_min) && (kw <= Target.kw_max)
                kw_req(Idx_Core, Idx_Turns, Idx_Wire) = kw;
            end
        end
    end

    % Remove invalid cores from the results
    for Idx_Core = size(Valid_Cores.Name, 2):-1:1
        if all(isnan(kw_req(Idx_Core, :, :)), 'all')
            kw_req(Idx_Core, :, :) = [];
            Bpk (Idx_Core, :) = [];
            Valid_Cores.Name(Idx_Core) = [];
            Valid_Cores.AeAw(Idx_Core) = [];
            Valid_Cores.Ae(Idx_Core) = [];
            Valid_Cores.Aw(Idx_Core) = [];
            Valid_Cores.le(Idx_Core) = [];
            Valid_Cores.lt(Idx_Core) = [];
            Valid_Cores.Ve(Idx_Core) = [];
        end
    end

    % Remove invalid turns from the results
    for Idx_Turns = size(Valid_Turns, 2):-1:1
        if all(isnan(kw_req(:, Idx_Turns, :)), 'all')
            kw_req(:, Idx_Turns, :) = [];
            Bpk (:, Idx_Turns) = [];
            Valid_Turns(Idx_Turns) = [];
        end
    end

    % Clear data
    clear Configurations_Dimensions Configurations_Max Aw_min kw;

    %% Losses

    Configurations_Dimensions = [size(Valid_Cores.Name, 2), size(Valid_Turns, 2), size(Valid_Wires.AWG, 2)];
    Configurations_Max = prod(Configurations_Dimensions);    
   
    % Preallocate matrices to receive the results
    P_Wire = ones(Configurations_Dimensions)*NaN;
    P_Core = ones(Configurations_Dimensions)*NaN;
    VL_Core = vpa(Volumetric_Losses(Target));
    P_Total = ones(Configurations_Dimensions)*NaN;
    
    % Calculate losses for each configuration
    for Idx_Cfg = 1:Configurations_Max
        % Get indexes
        [Idx_Core, Idx_Turns, Idx_Wire] = ind2sub(Configurations_Dimensions,Idx_Cfg);

        % Selected configuration is valid
        if (~isnan(kw_req(Idx_Core, Idx_Turns, Idx_Wire)))    
            % Calculate losses
            P_Wire(Idx_Core, Idx_Turns, Idx_Wire) = (Param.RCu * Valid_Turns(Idx_Turns) * Valid_Cores.lt(Idx_Core))/(Valid_Wires.Cond(Idx_Wire)*Valid_Wires.S_Cu(Idx_Wire))*Target.I_rms^2;
            P_Core(Idx_Core) = VL_Core*Valid_Cores.Ve(Idx_Core);
            
        end
        % Total losses
        % P_Total(Idx_Cfg) = P_Wire(Idx_Core, Idx_Turns, Idx_Wire) + P_Core(Idx_Core);
    end

    % Total losses
    %P_Total = P_Wire + P_Core;

    % Clear data
    clear Configurations_Dimensions Configurations_Max Idx_Core Idx_Turns Idx_Wire;
    
    %%
    [Min_Losses, Idx] = min(P_Total(:));
    [Idx_Core, Idx_Turns, Idx_Wire] = ind2sub(size(P_Total),Idx);

    Result.Core = Valid_Cores.Name(Idx_Core);
    Result.Turns = Valid_Turns(Idx_Turns);
    Result.Cond =  Valid_Wires.Cond(Idx_Wire);
    Result.AWG = Valid_Wires.AWG(Idx_Wire);
    Result.kw = kw_req(Idx_Core, Idx_Turns, Idx_Wire);
    Result.Pw = P_Wire(Idx_Core, Idx_Turns, Idx_Wire);
    Result.Pcore = P_Core(Idx_Core);
    Result.Pt = Min_Losses;

    % Clear data
    clear Min_Losses Idx Idx_Core Idx_Turns Idx_Wire;
end