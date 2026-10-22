function [Sim, IC, MC] = Initialize_Simulation_Public(F)
% INITIALIZE_SIMULATION  Simulation configuration and initial conditions
%   for the GHAME 6-DOF model.
%
%   [Sim, IC, MC] = Initialize_Simulation_Public(F)
%
%   Input
%     F     struct Helper function container
%
%   Outputs
%     Sim   struct Simulation configuration, constants and vehicle data
%     IC    struct Initial conditions, consistent with the Newton, Euler
%                  and kinematics blocks
%     MC    struct Monte Carlo settings (scope and levels of the aero
%                  uncertainty), used outside Simulink only

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

%% ========================================================================
% USER-DEFINED SIMULATION INPUTS
%
% These parameters define the simulation scenario, configuration,
% and initial conditions. Users are expected to modify ONLY this section.
% ========================================================================

% ---------------- Simulation and integration ----------------
Sim.sample_time             = 0.01;                     % Integration time step [s]
Sim.simulation.duration     = 90;                       % Simulation duration

% ---------------- Model toggles ----------------
Sim.mode.mEarth             = 1;            % 0: spherical Earth, 1: oblate Earth
Sim.mode.mwind              = 0;            % Enable wind model
Sim.mode.mturb              = 0;            % Enable turbulence
Sim.mode.mINS               = 0;            % INS errors, 0: ideal INS, 1: with INS errors
Sim.mode.mgps               = 0;            % GPS aiding of the INS, 0: without, 1: with
Sim.mode.mstar              = 0;            % Star tracker aiding of the INS, 0: without, 1: with
Sim.mode.muncertainty       = 1;            % Enable uncertain simulations

% ---------------- Actuator model ----------------
Sim.actuator.mact           = 2;            % Actuator model selector
Sim.actuator.zeta           = 0.70;         % Actuator damping ratio
Sim.actuator.omega          = 50;           % Actuator natural frequency [rad/s]
Sim.actuator.defl_lim_deg   = 20;           % Actuator deflection limit  [deg]
Sim.actuator.rate_lim_deg   = 50;           % Actuator rate limit        [deg/s]

% ------------------ Wind model ------------------
Sim.wind.turb_sigma         = 5.0;          % Wind turbulence intensity
Sim.wind.turb_length        = 150;          % Wind turbulence length    [m]
Sim.wind.tau_smooth         = 10;           % Mean wind smoothing time constant [s]

% ---------------- Filter models ----------------
Sim.filters.omega_AA        = 150;          % Anti-aliasing cutoff frequency [rad/s]
Sim.filters.omega_N         = 10;           % Noise filter natural frequency [rad/s]
Sim.filters.zeta_N          = 0.707;        % Noise filter damping ratio 
Sim.filters.delay_step      = 0;            % Effective actuator delay [s]

% --------- Aerodynamic uncertainty (Monte Carlo) ----------
Sim.uncertainty.N           = 5;            % Number of uncertainty runs
Sim.uncertainty.seed        = 67;           % Seed for reproducibility
MC.scope                    = {'control'};  % Groups and/or tables
MC.level                    = 0.3;          % One level, or one per scope entry

% Examples
%   {'control'}            0.3          All control derivatives +-30 %
%   {'static'}             0.2          All static stability derivatives +-20 %
%   {'damping'}            0.3          All rate damping derivatives +-30 %
%   {'longitudinal'}       0.2          Lift, drag and pitch tables +-20 %
%   {'lateral'}            0.2          Side force, roll and yaw tables +-20 %
%   {'lift', 'drag'}       0.1          Lift and drag tables +-10 %
%   {'pitch'}              0.2          Pitch moment tables +-20 %
%   {'CMA'}                0.3          Pitch stability derivative only +-30 %
%   {'all'}                0.1          All tables +-10 %
%   {'all', 'control'}     [0.1, 0.3]   All tables +-10 %, control derivatives +-30 %
%   {'all', 'CMQ'}         [0.2, 0]     All tables +-20 % except pitch damping

% ---------------- INS model ----------------
ins_seed                    = 1;            % Seed for the INS errors of this run
ins_frax                    = 0;            % Initial error scaling, XX0*(1+frax)

% ---------------- GPS ----------------
gps_seed                    = 1;            % Seed for the GPS errors of this run

% ---------------- Star tracker ----------------
star_seed                   = 1;            % Seed for the star tracker errors of this run

% ---------------- Launch location ----------------
Sim.launch.lonx             = 0;            % Geodetic longitude [deg]
Sim.launch.latx             = 15;            % Geodetic latitude  [deg]

% ---------------- Initial flight condition ----------------
Sim.initial.velocity.dvbe   = 1000;         % Initial airspeed [m/s]
Sim.initial.altitude.alt    = 10000;        % Initial altitude [m]

Sim.initial.attitude.psibdx = 0;            % Initial yaw   [deg]
Sim.initial.attitude.thtbdx = 3;            % Initial pitch [deg]
Sim.initial.attitude.phibdx = 0;            % Initial roll  [deg]

Sim.initial.aero.alpha_deg  = 3;            % Initial angle of attack [deg]
Sim.initial.aero.beta_deg   = 0;            % Initial sideslip angle [deg]

Sim.initial.rates.ppx       = 0;            % Initial roll rate  [deg/s]
Sim.initial.rates.qqx       = 0;            % Initial pitch rate [deg/s]
Sim.initial.rates.rrx       = 0;            % Initial yaw rate   [deg/s]

% ---------------- Vehicle mass state ----------------
Sim.vehicle.fuel.fraction   = 0.4;          % Fuel fraction (0 = empty, 1 = full)

% ----------------- Uncertainty Consistency Check ------------------

if Sim.uncertainty.N > 0 && Sim.mode.muncertainty == 0
    warning(['Uncertainty runs requested (N = %d) but ', ...
             'Sim.mode.muncertainty is disabled. ', ...
             'No Monte Carlo simulations will be executed.'], ...
             Sim.uncertainty.N);
end

if Sim.uncertainty.N == 0 && Sim.mode.muncertainty == 1
    warning(['Sim.mode.muncertainty is enabled but N = 0. ', ...
             'Only nominal simulation will be executed.']);
end

%%% Warning DO NOT change values from here onwards %%%

%% ========================================================================
% Constants
Sim.constants.PI    = 3.1415926536;
Sim.constants.DEG   = 57.2957795130823;
Sim.constants.RAD   = pi / 180;
Sim.constants.EPS   = 1e-10;
Sim.constants.SMALL = 1e-7;

RAD = Sim.constants.RAD;

Sim.simulation.Sim_Extender = 1;

%% Earth and environment model
Sim.earth.AGRAV    = 9.80675445;     % Standard gravity WGS-84
Sim.earth.GM       = 3.9860044e14;   % Earth gravitational parameter [m^3/s^2]
Sim.earth.WEII3    = 7.292115e-5;    % Earth rotation rate [rad/s];
Sim.earth.GW_CLONG = 0;              % Celestial longitude 

if Sim.mode.mEarth == 0
    Sim.earth.SMAJOR_AXIS = 6370987.308;   % Approximated to REARTH [m]
    Sim.earth.REARTH      = 6370987.308;   % Mean Earth radius [m]
    Sim.earth.FLATTENING  = 0;             % No flattening for spherical Earth model
    Sim.earth.C20         = 0;             % Point-mass gravity, consistent with the sphere
else
    Sim.earth.SMAJOR_AXIS = 6378137.0;     % WGS-84 semi-major axis (equatorial radius) [m]
    Sim.earth.REARTH      = 6370987.308;   % Mean Earth radius [m]
    Sim.earth.FLATTENING  = 3.33528106e-3; % WGS-84 flattening factor (f = (a-b)/a)
    Sim.earth.C20         = -4.8416685e-4; % Normalized second zonal harmonic (J2 = -sqrt(5)*C20)
end

%% ========================================================================
% Vehicle properties

% Vehicle mass properties
Sim.vehicle.mass.total_full = 136077;  % Maximum total vehicle mass
Sim.vehicle.mass.fuel_full  = 81646;   % Maximum fuel mass

Sim.vehicle.mass.structural = ...
    Sim.vehicle.mass.total_full - Sim.vehicle.mass.fuel_full; % Structural vehicle mass

% Current fuel and vehicle mass
Sim.vehicle.mass.fuel = ...
    Sim.vehicle.fuel.fraction * Sim.vehicle.mass.fuel_full;  % Actual fuel mass

Sim.vehicle.mass.current = ...
    Sim.vehicle.mass.structural + Sim.vehicle.mass.fuel;     % Actual vehicle mass

Sim.vehicle.inertia.IBBB0 = [ ...
    1.5730e+06  0            3.8000e+05;
    0           3.1600e+07   0;
    3.8000e+05  0            3.2540e+07 ]; % MOI at full fuel

Sim.vehicle.inertia.IBBB1 = [ ...
    1.1800e+06  0            2.4000e+05;
    0           1.9625e+07   0;
    2.4000e+05  0            2.0200e+07 ]; % MOI at empty fuel

Sim.vehicle.inertia.IBBBr = ...
    Sim.vehicle.fuel.fraction     * Sim.vehicle.inertia.IBBB0 + ...
    (1-Sim.vehicle.fuel.fraction) * Sim.vehicle.inertia.IBBB1; % Actual MOI

Sim.vehicle.geometry.refa = 557.42; % Reference surface area [m^2]
Sim.vehicle.geometry.refb = 24.38;  % Reference span         [m]
Sim.vehicle.geometry.refc = 22.86;  % Reference chord        [m]

Sim.vehicle.acowl = 27.27; % Engine cowl area [m^2]


%% ========================================================================
% ===================== INITIAL CONDITIONS ===============================
% ========================================================================

% Environment
IC.dvba = Sim.initial.velocity.dvbe;
IC.alt  = Sim.initial.altitude.alt;
IC.GW_CLONG = Sim.earth.GW_CLONG;

% Kinematics
IC.psibd = Sim.initial.attitude.psibdx * RAD;
IC.thtbd = Sim.initial.attitude.thtbdx * RAD;
IC.phibd = Sim.initial.attitude.phibdx * RAD;

IC.lon = Sim.launch.lonx*RAD;
IC.lat = Sim.launch.latx*RAD;

IC.TBD = F.mat3tr(IC.psibd, IC.thtbd, IC.phibd);
IC.TDI = F.cad_tdi84(IC.GW_CLONG + IC.lon, IC.lat);
IC.TBI = IC.TBD * IC.TDI;

% Mass and inertia
IC.vmassr = Sim.vehicle.mass.current;
IC.IBBBr = Sim.vehicle.inertia.IBBBr;

Sim.actuator.deflection_limit_rad      = Sim.actuator.defl_lim_deg * RAD;
Sim.actuator.deflection_rate_limit_rad = Sim.actuator.rate_lim_deg * RAD;

% INS sensor errors and initial error states. The ins block initializes
% its own error states from Sim.ins, IC.ESBI, IC.EVBI, IC.RICI are the same
% values for logging and plotting.
Sim.ins = ins_setup(ins_seed, ins_frax);
[IC.ESBI, IC.EVBI, IC.RICI] = F.init_ins(Sim.mode.mINS, Sim.ins.frax, Sim.ins.GAUSS0);

% GPS receiver, error model and filter data
Sim.gps = gps_setup(gps_seed);
if Sim.mode.mgps > 0 && Sim.mode.mINS == 0
    warning('Sim.mode.mgps is enabled but Sim.mode.mINS = 0. GPS has no effect on an ideal INS.');
end

% Star tracker data
Sim.star = startrack_setup(star_seed);
if Sim.mode.mstar > 0 && Sim.mode.mINS == 0
    warning('Sim.mode.mstar is enabled but Sim.mode.mINS = 0. The star tracker has no effect on an ideal INS.');
end

% Position and velocity
IC.WEII3 = Sim.earth.WEII3;

% Skew matrix for the cross products in the Newton block,
% vector IC.WEII (below) for TBI*WEII in the Euler block
IC.WEII_skew = [ ...
    0 -IC.WEII3 0;
    IC.WEII3 0 0;
    0 0 0 ];

IC.SBII = F.cad_in_geo84( ...
    IC.lon, IC.lat, IC.alt, ...
    Sim.earth.SMAJOR_AXIS, Sim.earth.FLATTENING,...
    IC.GW_CLONG);

IC.alpha = Sim.initial.aero.alpha_deg * RAD;
IC.beta  = Sim.initial.aero.beta_deg  * RAD;
IC.dvbe  = Sim.initial.velocity.dvbe;

IC.VBEB = [ ...
    cos(IC.alpha)*cos(IC.beta)*IC.dvbe;
    sin(IC.beta)*IC.dvbe;
    sin(IC.alpha)*cos(IC.beta)*IC.dvbe];

IC.VBED = IC.TBD' * IC.VBEB;

% Earth-relative velocity in inertial axes, initial state of the trim model
IC.VBEI = IC.TDI' * IC.VBED;

% Inertial velocity includes the transport by Earth rotation
IC.VBII = IC.VBEI + IC.WEII_skew * IC.SBII;
IC.dvbi = F.absolute(IC.VBII);

% Flight-path angles
POLAR = F.pol_from_cart(IC.VBED);
IC.psivd = POLAR(2);
IC.thtvd = POLAR(3);

% Angular rates
IC.WBEB = [ ...
    Sim.initial.rates.ppx;
    Sim.initial.rates.qqx;
    Sim.initial.rates.rrx ] * RAD;

IC.WEII = [0; 0; IC.WEII3];

IC.WBIB = IC.WBEB + IC.TBI * IC.WEII;

% Initialize trim values
IC.dp_trim = 0;
IC.dq_trim = 0;
IC.dr_trim = 0;
IC.dT_trim = 0.05;


end