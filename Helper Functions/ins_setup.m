function INS = ins_setup(seed, frax)
% INS_SETUP  Random INS errors for one simulation run.
%
%   INS = ins_setup(seed, frax)
%
%   Draws all per run random values of the INS model with a fixed seed, so
%   a Monte Carlo run can be reproduced. The result goes into Sim.ins before
%   the run. The gyro random walk is not drawn here, it enters ins as the
%   input ZWG from a Random Number block, whose seed can be set from
%   INS.seed_walk.
%
%   Inputs
%     seed   [1x1]  Random seed of this run (-)
%     frax   [1x1]  Initial error scaling, XX0*(1+frax), default 0 (-)
%
%   Outputs
%     INS    struct Sensor errors, initial error sample and seeds
%       .EMISG      [3x1]  Gyro misalignment (rad)
%       .ESCALG     [3x1]  Gyro scale factor (-)
%       .EBIASG     [3x1]  Gyro bias (rad/s)
%       .EUNBG      [3x1]  Gyro spin axis mass unbalance (rad/s per m/s^2)
%       .EWALKG     [1x1]  Gyro random walk intensity (rad/sqrt(s))
%       .EMISA      [3x1]  Accel misalignment (rad)
%       .ESCALA     [3x1]  Accel scale factor (-)
%       .EBIASA     [3x1]  Accel bias (m/s^2)
%       .EWALKA     [3x1]  Accel random bias (m/s^2)
%       .GAUSS0     [9x1]  Unit Gaussian sample for init_ins (-)
%       .frax       [1x1]  Initial error scaling (-)
%       .seed_walk  [1x1]  Seed for the random walk Random Number block (-)

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

if nargin < 2
    frax = 0;
end
rng(seed, 'twister');
% Gyro errors, drawn from their 1-sigma values
INS.EMISG  = 1.1e-4  * randn(3,1);      % 1 sigma 1.1e-4 rad
INS.ESCALG = 2.0e-5  * randn(3,1);      % 1 sigma 20 ppm
INS.EBIASG = 1.0e-6  * randn(3,1);      % 1 sigma 0.2 deg/hr
INS.EUNBG  = [0; 0; 0];
INS.EWALKG = 0.10*pi/180/60;            % 0.1 deg/sqrt(hr)
% Accelerometer errors, drawn from their 1-sigma values
INS.EMISA  = 1.1e-4  * randn(3,1);      % 1 sigma 1.1e-4 rad
INS.ESCALA = 5.0e-4  * randn(3,1);      % 1 sigma 500 ppm
INS.EBIASA = 3.56e-3 * randn(3,1);      % 1 sigma 3.56e-3 m/s^2
INS.EWALKA = [0; 0; 0];
% Initial error states
INS.GAUSS0 = randn(9,1);
INS.frax   = frax;
% Seed for the random walk noise block
INS.seed_walk = randi(2^31 - 2);
end