function GPS = gps_setup(seed)
% GPS_SETUP  GPS receiver and filter data for one simulation run.
%
%   GPS = gps_setup(seed)
%
%   Fixed receiver and filter settings plus the per run random values
%   (pseudo-range biases, initial clock bias, initial Markov states), drawn
%   with a fixed seed. The result goes into Sim.gps before the run. The
%   Markov noises are driven during the run by the ZGPS input of ins from a
%   Random Number block, whose seed can be set from GPS.seed_noise.
%
%   The numbers below follow the CADAC++ hypersonic vehicle example Hyper6
%
%   Inputs
%     seed   [1x1]  Random seed of this run (-)
%
%   Outputs
%     GPS    struct Receiver, error model and filter data

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

% Own random stream, the global generator state is restored at the end
rng_state = rng;
rng(seed, 'twister');
% ---------------- Receiver ----------------
GPS.almanac_time = 0;                   % Time since almanac epoch at sim start (s)
GPS.del_rearth   = 50e3;                % Added to Earth radius for line of sight (m)
GPS.gps_acqtime  = 10;                  % Signal acquisition time (s)
GPS.gps_step     = 0.25;                % Update interval, 4 Hz (s)
% ---------------- Measurement errors ----------------
sig_pr_bias = 0.842;                    % Pseudo-range bias, 1 sigma (m)
sig_ucbias  = 3;                        % Initial user clock bias, 1 sigma (m)
sig_pr      = 0.25;                     % Pseudo-range noise, 1 sigma (m)
sig_dr      = 0.03;                     % Range rate noise, 1 sigma (m/s)
sig_ucfreq  = 0.1;                      % User clock frequency error, 1 sigma (m/s)
tau_pr      = 500;                      % Pseudo-range noise correlation time (s)
tau_dr      = 500;                      % Range rate noise correlation time (s)
tau_ucfreq  = 100;                      % Clock frequency correlation time (s)
% ---------------- Filter ----------------
GPS.uctime_cor = 100;                   % User clock correlation time (s)
GPS.ppos    = 5;                        % Initial P, position 1 sigma (m)
GPS.pvel    = 0.2;                      % Initial P, velocity 1 sigma (m/s)
GPS.pclockb = 3;                        % Initial P, clock bias 1 sigma (m)
GPS.pclockf = 1;                        % Initial P, clock frequency 1 sigma (m/s)
GPS.qpos    = 0.1;                      % Q, position 1 sigma (m)
GPS.qvel    = 0.01;                     % Q, velocity 1 sigma (m/s)
GPS.qclockb = 0.5;                      % Q, clock bias 1 sigma (m)
GPS.qclockf = 0.1;                      % Q, clock frequency 1 sigma (m/s)
GPS.rpos    = 1;                        % R, pseudo-range 1 sigma (m)
GPS.rvel    = 0.1;                      % R, range rate 1 sigma (m/s)
GPS.factp   = 0;                        % Scaling P*(1+factp)^2 (-)
GPS.factq   = 0;                        % Scaling Q*(1+factq)^2 (-)
GPS.factr   = 0;                        % Scaling R*(1+factr)^2 (-)
% ---------------- Per run random values ----------------
GPS.pr_bias    = sig_pr_bias * randn(4,1);
GPS.ucbias0    = sig_ucbias  * randn;
GPS.sig_markov = [sig_pr * ones(4,1); sig_dr * ones(4,1); sig_ucfreq];
GPS.tau_markov = [tau_pr * ones(4,1); tau_dr * ones(4,1); tau_ucfreq];
GPS.markov0    = GPS.sig_markov .* randn(9,1);
GPS.seed_noise = randi(2^31 - 2);
rng(rng_state);
end