function TDI = cad_tdi84(lon_cel, lat)
% CAD_TDI84  Transformation matrix from inertial to geodetic coordinates.
%
%   TDI = cad_tdi84(lon_cel, lat)
%
%   The geodetic frame D is north, east, down at the vehicle position.
%   Using the celestial longitude accounts for Earth rotation directly.
%   At the initial epoch lon_cel = GW_CLONG + lon, during the simulation
%   lon_cel = GW_CLONG + WEII3*t + lon.
%
%   Inputs
%     lon_cel  [1x1]  Celestial longitude, measured in the inertial frame (rad)
%     lat      [1x1]  Geodetic latitude (rad)
%
%   Output
%     TDI      [3x3]  Transformation matrix, v_D = TDI * v_I

%   Author:    Thomas Harald Müller
%   GitHub:    https://github.com/ThomasHaraldMueller
%   Copyright: (c) 2026 Thomas Harald Müller, Hamburg, Germany
%   Free to use and modify, provided this notice is kept and the author is credited.

slat = sin(lat);
clat = cos(lat);
slon = sin(lon_cel);
clon = cos(lon_cel);

TDI = [ -slat*clon,  -slat*slon,   clat;     % North
        -slon,        clon,        0;        % East
        -clat*clon,  -clat*slon,  -slat ];   % Down

end