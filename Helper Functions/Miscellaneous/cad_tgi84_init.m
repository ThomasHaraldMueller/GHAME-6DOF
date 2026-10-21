function TGI = cad_tgi84_init(lon, lat, alt, GW_CLONG, SMAJOR_AXIS, FLATTENING)
% CAD_TGI84_INIT  Transformation matrix from inertial to geocentric
%   coordinates at the initial epoch.
%
%   TGI = cad_tgi84_init(lon, lat, alt, GW_CLONG, SMAJOR_AXIS, FLATTENING)
%
%   Uses the celestial longitude GW_CLONG + lon, so it is valid at t = 0 only
%   and is intended for initialization. The geocentric frame G is the
%   geodetic frame rotated about its east axis by the deflection dd between
%   the ellipsoid normal and the radius vector. dd and the ellipsoid radius
%   use series expansions in the flattening, accurate to about 1e-5 rad and
%   0.1 m below 100 km altitude.
%
%   Inputs
%     lon          [1x1]  Geodetic longitude (rad)
%     lat          [1x1]  Geodetic latitude (rad)
%     alt          [1x1]  Altitude above the ellipsoid (m)
%     GW_CLONG     [1x1]  Greenwich longitude at the initial epoch (rad)
%     SMAJOR_AXIS  [1x1]  Ellipsoid semi-major axis (m)
%     FLATTENING   [1x1]  Ellipsoid flattening
%
%   Output
%     TGI          [3x3]  Transformation matrix, inertial to geocentric

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

a = SMAJOR_AXIS;
f = FLATTENING;

% Inertial to geodetic at the initial epoch
TDI = cad_tdi84(GW_CLONG + lon, lat);

% Ellipsoid radius at the geodetic latitude
r0 = a * (1 - f*(1 - cos(2*lat))/2 + 5*f^2*(1 - cos(4*lat))/16);

% Deflection between geodetic and geocentric vertical
dd = f * sin(2*lat) * (1 - f/2 - alt/r0);

% Rotation about the east axis from geodetic to geocentric
TGD = [cos(dd), 0, -sin(dd);
       0,       1,  0;
       sin(dd), 0,  cos(dd)];

TGI = TGD * TDI;

end