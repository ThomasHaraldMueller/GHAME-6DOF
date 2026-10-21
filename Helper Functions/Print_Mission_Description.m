function Print_Mission_Description(Sim, IC, MC)
% PRINT_MISSION_DESCRIPTION  Prints the flight condition, active simulation
%   modes and the uncertainty scenario to the command window.
%
%   Print_Mission_Description(Sim, IC, MC)
%
%   Inputs
%     Sim   struct Simulation configuration and mode switches
%     IC    struct Initial conditions
%     MC    struct Monte Carlo settings (scope and levels), optional

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

line = repmat('=', 1, 60);
fprintf('\n%s\n', line);

%% Flight condition
fprintf('Simulation mission setup:\n');
fprintf('  - Airspeed:           %.1f m/s\n', IC.dvba);
fprintf('  - Altitude:           %.1f m\n',   IC.alt);
fprintf('  - Fuel fraction:      %.2f\n',     Sim.vehicle.fuel.fraction);
fprintf('  - Duration:           %.1f s\n\n', Sim.simulation.duration);

%% Models and modes
if Sim.mode.mEarth > 0
    earth = 'WGS-84 ellipsoid with J2 gravity';
else
    earth = 'Spherical, point-mass gravity';
end

fprintf('Simulation modes:\n');
fprintf('  - Earth model:        %s\n', earth);
fprintf('  - Wind:               %s\n', on_off(Sim.mode.mwind, 'HWM14, smoothed'));
fprintf('  - Turbulence:         %s\n', on_off(Sim.mode.mturb, 'Dryden, normal gust'));
fprintf('  - INS errors:         %s\n', on_off(Sim.mode.mINS));
fprintf('  - Uncertainty:        %s\n', on_off(Sim.mode.muncertainty));

%% Uncertainty scenario
if Sim.mode.muncertainty > 0
    fprintf('  - Monte Carlo:        %d runs, seed %d\n', ...
            Sim.uncertainty.N, Sim.uncertainty.seed);

    if nargin > 2
        fprintf('  - Uncertain tables:   %s\n', scope_text(MC));
    end
end

fprintf('%s\n\n', line);

end

%% Local functions

function str = on_off(flag, detail)
% ON_OFF  'Enabled' or 'Disabled', with an optional detail when enabled
    if flag > 0
        str = 'Enabled';
        if nargin > 1
            str = sprintf('Enabled (%s)', detail);
        end
    else
        str = 'Disabled';
    end
end

function str = scope_text(MC)
% SCOPE_TEXT  Scope entries with their levels, e.g. 'all +-10 %, control +-30 %'
    scope = cellstr(MC.scope);
    level = MC.level(:);
    if isscalar(level)
        level = repmat(level, numel(scope), 1);
    end

    parts = cell(1, numel(scope));
    for i = 1:numel(scope)
        parts{i} = sprintf('%s +-%g %%', scope{i}, 100 * level(i));
    end
    str = strjoin(parts, ', ');
end