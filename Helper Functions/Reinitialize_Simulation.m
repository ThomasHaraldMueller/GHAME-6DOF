function IC = Reinitialize_Simulation(IC, Sim, F)
% REINITIALIZE_SIMULATION  Updates the initial conditions after trim.
%
%   IC = Reinitialize_Simulation(IC, Sim, F)
%
%   Takes the trimmed angle of attack and pitch angle, rebuilds attitude,
%   body, Earth-relative and inertial velocity and the flight path angles,
%   and activates the Earth rotation for the simulation. Sideslip and body
%   rates are set to zero, position, speed, heading and roll are kept.
%
%   Inputs
%     IC    struct Initial conditions including the trim results
%                  IC.alpha_trim and IC.thtbd_trim (rad)
%     Sim   struct Earth rotation rate
%     F     struct Helper functions (mat3tr, absolute, pol_from_cart)
%
%   Output
%     IC    struct Updated initial conditions

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

%% Trim results
IC.alpha = IC.alpha_trim;
IC.thtbd = IC.thtbd_trim;
IC.beta  = 0;                       % Symmetric trim, no sideslip

%% Attitude
IC.TBD = F.mat3tr(IC.psibd, IC.thtbd, IC.phibd);
IC.TBI = IC.TBD * IC.TDI;

%% Earth rotation, activated for the simulation
IC.WEII3     = Sim.earth.WEII3;
IC.WEII      = [0; 0; IC.WEII3];
IC.WEII_skew = [ 0,        -IC.WEII3, 0;
                 IC.WEII3,  0,        0;
                 0,         0,        0 ];

%% Velocities
IC.VBEB = IC.dvbe * [ cos(IC.alpha) * cos(IC.beta);
                      sin(IC.beta);
                      sin(IC.alpha) * cos(IC.beta) ];
IC.VBED = IC.TBD' * IC.VBEB;

% Earth-relative velocity in inertial axes, initial state of the trim model
IC.VBEI = IC.TDI' * IC.VBED;

% Inertial velocity includes the transport by Earth rotation
IC.VBII = IC.VBEI + IC.WEII_skew * IC.SBII;
IC.dvbi = F.absolute(IC.VBII);

%% Flight path angles
POLAR    = F.pol_from_cart(IC.VBED);
IC.psivd = POLAR(2);
IC.thtvd = POLAR(3);

%% Body rates, zero relative to Earth in steady trim
IC.WBEB = [0; 0; 0];
IC.WBIB = IC.WBEB + IC.TBI * IC.WEII;

end