% MAIN  Master script of the GHAME 6-DOF simulation.
%
%   Initializes the model, trims the vehicle and runs either a nominal
%   simulation or a Monte Carlo study with aerodynamic uncertainty.
%   Settings are made in Initialize_Simulation_Public.

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

clc; close all; clear

addpath('Helper Functions', 'Helper Functions/Miscellaneous', 'Simulink', 'TrimMod')

mdl      = 'GHAME_Public';
SAVE_EPS = true;                        % Export figures as EPS

% Data and helper functions
Data = Load_Aero_Data();

F = struct( ...
    'mat3tr',        @mat3tr, ...
    'cad_tdi84',     @cad_tdi84, ...
    'cad_in_geo84',  @cad_in_geo84, ...
    'init_ins',      @init_ins, ...
    'pol_from_cart', @pol_from_cart, ...
    'absolute',      @absolute);

% Initialization and trim
[Sim, IC, MC] = Initialize_Simulation_Public(F);
Print_Mission_Description(Sim, IC, MC);
%%
% Arguments: console output, force retrim
[Sim, IC] = Trim_Run_jj(true, true, Sim, IC);

% Initial conditions from the trimmed alpha and theta
IC = Reinitialize_Simulation(IC, Sim, F);

%% Nominal run, also the reference for Monte Carlo
DataPlant = Data;
out_nom   = sim(mdl);

%% Monte Carlo, one perturbed aero database per run
if Sim.mode.muncertainty > 0 && Sim.uncertainty.N > 0
    for k = 1:Sim.uncertainty.N
        UNC       = Create_Uncertainty(Data, MC.scope, MC.level, Sim.uncertainty.seed + k);
        DataPlant = UNC.Perturbed;      %#ok
        out(k)    = sim(mdl);           %#ok
    end
    DataPlant = Data;                   % Back to nominal after MC
end