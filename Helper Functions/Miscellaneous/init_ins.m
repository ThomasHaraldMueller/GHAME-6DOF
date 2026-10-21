function [ESBI, EVBI, RICI] = init_ins(mINS, frax, GAUSS)
% INIT_INS  Initial INS errors from the GPS quality covariance.
%
%   [ESBI, EVBI, RICI] = init_ins(mINS, frax, GAUSS)
%
%   A correlated Gaussian error vector XX0 = L * GAUSS is formed, with
%   PP0 = L * L' the initial covariance of a GPS equipped aircraft and
%   GAUSS ~ N(0, I) drawn outside, so the function is deterministic and
%   can run in a Simulink MATLAB Function block. Its parts are the
%   position, velocity and tilt errors, scaled by (1 + frax). Without INS
%   errors (perfect transfer alignment) all outputs are zero.
%
%   Inputs
%     mINS   [1x1]  INS mode, 0 = no INS errors
%     frax   [1x1]  Uncertainty factor applied to all errors (-)
%     GAUSS  [9x1]  Unit Gaussian sample, e.g. randn(9,1) (-)
%
%   Outputs
%     ESBI   [3x1]  Initial position error (m)
%     EVBI   [3x1]  Initial velocity error (m/s)
%     RICI   [3x1]  Initial tilt error (rad)

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

persistent L
ESBI = zeros(3,1);
EVBI = zeros(3,1);
RICI = zeros(3,1);
% Perfect transfer alignment, no INS errors
if mINS == 0
    return;
end
% Cholesky factor of the initial covariance, computed once
if isempty(L)
    % Initial covariance [pos (m); vel (m/s); tilt (mrad)]
    PP0 = [
          20.701,     0.12317,    0.10541,    6.3213E-02, 2.2055E-03, 1.7234E-03, 1.0633E-03, 3.4941E-02,-3.5179E-02
          0.12317,    20.696,    -0.27174,    4.8366E-03, 5.9463E-02,-1.3367E-03,-3.4903E-02, 2.6112E-03,-4.2663E-02
          0.10541,   -0.27174,    114.12,     5.6373E-04,-8.3147E-03, 5.4059E-02, 1.5496E-02, 7.6463E-02,-3.5302E-03
          6.3213E-02, 4.8366E-03, 5.6373E-04, 1.9106E-03, 8.0945E-05, 1.9810E-06, 2.5755E-04, 2.8346E-03,-5.6482E-04
          2.2055E-03, 5.9463E-02,-8.3147E-03, 8.0945E-05, 1.7201E-03,-1.5760E-05,-2.8341E-03, 2.6478E-04,-1.0781E-03
          1.7234E-03,-1.3367E-03, 5.4059E-02, 1.9810E-06,-1.5760E-05, 3.0070E-03, 4.1963E-04,-1.3297E-04, 4.1190E-05
          1.0638E-03,-3.4903E-02, 1.5496E-02, 2.5755E-04,-2.8341E-03, 4.1963E-04, 5.4490E-02,-1.8695E-03, 8.9868E-04
          3.4941E-02, 2.6112E-03, 7.6463E-02, 2.8346E-03, 2.6478E-04,-1.3297E-04,-1.8695E-03, 5.2819E-02, 1.0990E-02
         -3.5179E-02,-4.2663E-02,-3.5302E-03,-5.6482E-04,-1.0781E-03, 4.1190E-05, 8.9868E-04, 1.0990E-02, 0.1291
    ];
    PP0 = (PP0 + PP0.') / 2;            % Remove the small (1,7)/(7,1) asymmetry
    L   = chol(PP0, 'lower');
end
% Correlated error vector [pos (m); vel (m/s); tilt (mrad)]
XX0 = L * GAUSS(:);
s   = 1 + frax;
ESBI = s * XX0(1:3);
EVBI = s * XX0(4:6);
RICI = s * XX0(7:9) * 1e-3;             % mrad to rad
end