function [Sim, IC] = Trim_Run_jj(doPrint, forceTrim, Sim, IC)
% TRIM_RUN_JJ  Longitudinal steady-state trim of the GHAME vehicle with jj_trim.
%
%   [Sim, IC] = Trim_Run_jj(doPrint, forceTrim, Sim, IC)
%
%   Trims the Simulink model 'GHAME_Trim' for wings-level flight at the
%   initial speed and altitude. Free variables are angle of attack, pitch
%   angle, elevator and throttle, set so that alpha_dot, q_dot, the vertical
%   acceleration and V_dot vanish. The trim model integrates position and
%   velocity relative to the rotating Earth in Earth-fixed axes, which
%   coincide with the inertial axes at t = 0.
%
%   The vertical acceleration is taken as the inertial axis closest to the
%   local vertical. This is exact when the launch point lies on an inertial
%   axis, e.g. at 0 deg latitude and longitude, and an approximation
%   otherwise.
%
%   The trim is only recomputed when one of the trim-defining conditions or
%   the saved trim model has changed since the last run (or when forceTrim
%   is true). Otherwise the cached results are loaded. Changes to the aero
%   or engine data are not detected and need forceTrim = true. Call Reinitialize_Simulation afterwards to
%   rebuild the initial conditions from the trim results.
%
%   Inputs
%     doPrint    [1x1]  Console output on (default true)
%     forceTrim  [1x1]  Run the trim even if the conditions are unchanged
%                       (default false)
%     Sim        struct Simulation configuration
%     IC         struct Initial conditions from Initialize_Simulation
%
%   Outputs
%     Sim        struct Unchanged, returned for a uniform call interface
%     IC         struct Initial conditions with the trim results
%                       alpha_trim, thtbd_trim (rad)
%                       dp_trim, dq_trim, dr_trim (rad), dT_trim

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

if nargin < 1 || isempty(doPrint),   doPrint   = true;  end
if nargin < 2 || isempty(forceTrim), forceTrim = false; end

RAD   = Sim.constants.RAD;
model = 'GHAME_Trim';

trimFileState = 'last_trim_state.mat';
trimFileCond  = 'last_trim_conditions.mat';

%% Trim-defining conditions
% A change in any of these invalidates the cached trim.
current.alt      = Sim.initial.altitude.alt;
current.dvbe     = Sim.initial.velocity.dvbe;
current.fuelFrac = Sim.vehicle.fuel.fraction;
current.mEarth   = Sim.mode.mEarth;          % Gravity model
current.latx     = Sim.launch.latx;          % Gravity and Coriolis vary with latitude
current.psibdx   = Sim.initial.attitude.psibdx;  % Vertical Coriolis term depends on heading

% Any saved change to the trim model also invalidates the cache
modelFile         = dir(which(model));
current.modelDate = modelFile.datenum;

%% Cached trim
runTrim = true;

if ~forceTrim && exist(trimFileCond, 'file') && exist(trimFileState, 'file')
    load(trimFileCond, 'lastCond');
    if isequal(current, lastCond)
        tmp = load(trimFileState, 'IC');

        % Only the trim results are copied, the rest of IC stays current
        IC.alpha_trim = tmp.IC.alpha_trim;
        IC.thtbd_trim = tmp.IC.thtbd_trim;
        IC.dp_trim    = tmp.IC.dp_trim;
        IC.dq_trim    = tmp.IC.dq_trim;
        IC.dr_trim    = tmp.IC.dr_trim;
        IC.dT_trim    = tmp.IC.dT_trim;

        runTrim = false;
        if doPrint
            fprintf('Trim skipped, using cached trim results.\n');
        end
    end
end

if ~runTrim
    return
end

if doPrint
    fprintf('Running trim routine...\n');
end

%% State indices from the compiled model
sys        = linmod(model);
stateNames = sys.StateName;

idx_rot = find_states(stateNames, 'Inertial Body Rotation', 3);
idx_pos = find_states(stateNames, 'Inertial Position',      3);
idx_vel = find_states(stateNames, 'Inertial Velocity',      3);
idx_vab = find_states(stateNames, 'V_alpha_beta',           3);
idx_eul = find_states(stateNames, 'Euler Angles',           3);

idx_q     = idx_rot(2);
% Velocity component closest to the local vertical (down row of TDI)
[~, i_rad] = max(abs(IC.TDI(3, :)));
idx_vrad   = idx_vel(i_rad);
idx_V     = idx_vab(1);
idx_alpha = idx_vab(2);
idx_theta = idx_eul(2);

%% Initial state
% Position and velocity relative to the Earth, in Earth-fixed axes that
% coincide with the inertial axes at t = 0.
nx = numel(stateNames);
x0 = zeros(nx, 1);

x0(idx_rot) = IC.WBIB;
x0(idx_pos) = IC.SBII;
x0(idx_vel) = IC.VBEI;            % Earth-relative, without the Earth surface speed
x0(idx_vab) = [IC.dvba; IC.alpha; IC.beta];
x0(idx_eul) = [IC.phibd; IC.thtbd; IC.psibd];

%% Initial input guess [aileron; elevator; rudder; throttle]
u0 = [IC.dp_trim; IC.dq_trim; IC.dr_trim; IC.dT_trim];
nu = numel(u0);

%% Trim problem
d0 = zeros(nx, 1);              % Desired derivatives, steady flight

ny      = 1;                    % Outputs not used, required by jj_trim
y0      = zeros(ny, 1);
iy_trim = [];

ix = [idx_alpha; idx_theta];    % Free states
iu = [2; 4];                    % Free inputs, elevator and throttle
id = [idx_alpha;                % alpha_dot = 0
      idx_q;                    % q_dot     = 0
      idx_vrad;                 % Vertical acceleration = 0
      idx_V];                   % V_dot     = 0

if numel(ix) + numel(iu) ~= numel(id)
    error('Trim_Run_jj:illPosed', ...
          'Trim problem ill-posed, %d free variables for %d conditions.', ...
          numel(ix) + numel(iu), numel(id));
end

%% Names for diagnostics
x_nam = cellstr(stateNames);
u_nam = {'Aileron_Com'; 'Elevator_Com'; 'Rudder_Com'; 'Throttle_Com'};
d_nam = strcat(x_nam, '_dot');

%% Bounds
% Surfaces from the actuator limits, throttle as in thrust.m
ctrlLim = Sim.actuator.deflection_limit_rad;
u_min   = [-ctrlLim; -ctrlLim; -ctrlLim; 0.05];
u_max   = [ ctrlLim;  ctrlLim;  ctrlLim; 2   ];

% Alpha within the range of the aerodynamic tables
x_min = -inf(nx, 1);
x_max =  inf(nx, 1);
x_min(idx_alpha) = -3 * RAD;
x_max(idx_alpha) = 21 * RAD;

%% jj_trim
options.n_iter_max     = 50;
options.cost_tbg       = 1e-9;
options.CompileFlag    = 1;
options.EnableMessages = doPrint;

if doPrint
    fprintf('jj_trim: %s, nx = %d, nu = %d, ny = %d\n', model, nx, nu, ny);
end

[x_tr, u_tr, ~, ~, err, info] = jj_trim( ...
    model, ...
    x0, u0, d0, y0, ...
    ix, iu, id, iy_trim, ...
    x_nam, u_nam, d_nam, {}, ...
    [], [], [], [], ...
    options, ...
    x_min, x_max, u_min, u_max);

if err ~= 0
    error('Trim_Run_jj:failed', 'jj_trim failed, exit code %d.', info.ExitCode);
end

%% Results
IC.alpha_trim = x_tr(idx_alpha);
IC.thtbd_trim = x_tr(idx_theta);

IC.dp_trim = u_tr(1);
IC.dq_trim = u_tr(2);
IC.dr_trim = u_tr(3);
IC.dT_trim = u_tr(4);

if doPrint
    fprintf('\nTrim successful\n');
    fprintf('  alpha    : %7.3f deg\n', IC.alpha_trim / RAD);
    fprintf('  theta    : %7.3f deg\n', IC.thtbd_trim / RAD);
    fprintf('  elevator : %7.3f deg\n', IC.dq_trim    / RAD);
    fprintf('  throttle : %7.3f\n\n',   IC.dT_trim);
end

%% Cache
lastCond = current; 
save(trimFileCond,  'lastCond');
save(trimFileState, 'Sim', 'IC');

end

%% Local functions

function idx = find_states(names, pattern, n)
% FIND_STATES  Indices of the states whose name contains pattern, with a
%   clear error if the block was renamed or removed from the trim model.
    idx = find(contains(names, pattern));
    if numel(idx) ~= n
        error('Trim_Run_jj:states', ...
              'Expected %d states containing ''%s'' in the trim model, found %d.', ...
              n, pattern, numel(idx));
    end
end