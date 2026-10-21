function SBII = cad_in_geo84(lon, lat, alt, SMAJOR_AXIS, FLATTENING, GW_CLONG)
% CAD_IN_GEO84  Inertial position from geodetic longitude, latitude and
%   altitude on the WGS-84 ellipsoid at the initial epoch.
%
%   SBII = cad_in_geo84(lon, lat, alt, SMAJOR_AXIS, FLATTENING, GW_CLONG)
%
%   Exact conversion via the prime vertical radius of curvature. 
%   The celestial longitude is GW_CLONG + lon, so the
%   result is valid at t = 0 and is intended for initialization.
%
%   Inputs
%     lon          [1x1]  Geodetic longitude (rad)
%     lat          [1x1]  Geodetic latitude (rad)
%     alt          [1x1]  Altitude above the ellipsoid (m)
%     SMAJOR_AXIS  [1x1]  Ellipsoid semi-major axis (m)
%     FLATTENING   [1x1]  Ellipsoid flattening
%     GW_CLONG     [1x1]  Greenwich longitude at the initial epoch (rad)
%
%   Output
%     SBII         [3x1]  Position in inertial coordinates (m)

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

a  = SMAJOR_AXIS;
f  = FLATTENING;
e2 = 2*f - f^2;                         % First eccentricity squared

lon_cel = GW_CLONG + lon;               % Celestial longitude at t = 0

N = a / sqrt(1 - e2 * sin(lat)^2);      % Prime vertical radius of curvature

SBII = [ (N + alt)          * cos(lat) * cos(lon_cel);
         (N + alt)          * cos(lat) * sin(lon_cel);
         (N * (1 - e2) + alt) * sin(lat) ];

end