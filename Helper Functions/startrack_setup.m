function STAR = startrack_setup(seed)
% STARTRACK_SETUP  Star tracker data for one simulation run.
%
%   STAR = startrack_setup(seed)
%
%   Fixed tracker settings plus the per run random values (azimuth and
%   elevation biases, initial Markov states), drawn with a fixed seed. The
%   result goes into Sim.star before the run. The Markov noises are driven
%   during the run by the ZST input of ins from a Random Number block, whose
%   seed can be set from STAR.seed_noise.
%
%   The numbers below are placeholders. Check them against your CADAC input
%   deck or the tracker you want to represent.
%
%   Inputs
%     seed   [1x1]  Random seed of this run (-)
%
%   Outputs
%     STAR   struct Tracker and error model data

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

% Own random stream, the global generator state is restored at the end
rng_state = rng;
rng(seed, 'twister');
% ---------------- Tracker ----------------
STAR.star_el_min   = 10;                % Minimum star elevation over horizon (deg)
STAR.startrack_alt = 5000;              % Altitude above which tracking works (m)
STAR.star_acqtime  = 10;                % Initial acquisition time (s)
STAR.star_step     = 0.25;              % Fix interval (s)
% ---------------- Measurement errors ----------------
sig_bias  = 5e-6;                       % Azimuth and elevation bias, 1 sigma, about 1 arcsec (rad)
sig_noise = 5e-6;                       % Azimuth and elevation noise, 1 sigma (rad)
tau_noise = 1;                          % Noise correlation time (s)
% ---------------- Per run random values ----------------
STAR.az_bias    = sig_bias * randn(3,1);
STAR.el_bias    = sig_bias * randn(3,1);
STAR.sig_markov = sig_noise * ones(6,1);
STAR.tau_markov = tau_noise * ones(6,1);
STAR.markov0    = STAR.sig_markov .* randn(6,1);
STAR.seed_noise = randi(2^31 - 2);
rng(rng_state);
end